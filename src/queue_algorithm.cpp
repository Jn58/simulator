#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>
#include "error.h"


namespace ClusterSimulator {
	void GeneAlgorithm::enqueJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for(auto job : jobs)
		{
			if(job->state==JobState::WAIT)
				if(!run_job(job))
				{ 
					for (auto& p : population)
					{
						p.enqueJob(job);
					}
					length += 1;
				}
		}
	}
	void GeneAlgorithm::deleteJobs(std::vector<std::shared_ptr<Job>>& executed_jobs)
	{
		{
			for (auto& p : population)
			{
				p.chromosomeDeleteJobs(executed_jobs);
			}
			length -= executed_jobs.size();
		}
	}
	GeneAlgorithm::Chromosome& GeneAlgorithm::getBestChromosome()
	{
		return population[0];
	}
	bool GeneAlgorithm::run_job(std::shared_ptr<Job> job)
{
		std::vector<Host*> hosts{ job->get_eligible_hosts() };
		if (hosts.empty())
		{
			return false;
		}

		Host* best_host = *std::min_element(hosts.begin(), hosts.end(), 
			[=](const Host* a, const Host* b)
			{
				return a->get_expected_run_time(*job) < b->get_expected_run_time(*job);
			});
		best_host->execute_job(*job);
		return true;
	}
	void GeneAlgorithm::exec()
	{
		if (length != 0)
		{
			Chromosome& best_chromosome{ getBestChromosome() };
			std::map<Host*, Chromosome::HostInfo>& best_hosts{ best_chromosome.hosts };
			std::vector<std::shared_ptr<Job>> excuted_jobs;
			for (auto it = best_hosts.begin(); it != best_hosts.end(); ++it)
			{
				Host* host = it->first;
				Chromosome::HostInfo& host_info{ it->second };
				if(host_info.count>0)
				{
					if (host != host_info.first_job_gene->host_)
					{
						error("host and first job's host must be same");
					}
					std::shared_ptr<Job> job{ host_info.first_job_gene->job_ };
					if (host->is_executable(*job))
					{
						host->execute_job(*job);
						excuted_jobs.push_back(job);
					}
				}
				else
				{
					error("host info count must not be 0");
				}
			}
			deleteJobs(excuted_jobs);
		}
	}
	bool GeneAlgorithm::check(std::vector<std::shared_ptr<Job>>& jobs)
	{
		size_t count = 0;
		size_t pedding_n = 0;
		for (auto job : jobs)
		{
			if (job->state == JobState::PEND) ++pedding_n;
			else if (job->state == JobState::WAIT) ++count;
		}
		if (pedding_n != length)
		{
			error("pedding_n and length must be same");
		}

		if (pedding_n + count != jobs.size())
		{
			for (auto job : jobs)
			{
				if (job->state != JobState::PEND && job->state != JobState::WAIT)
					printf("%d\n", job->state);
			}
			error("wrong state");
		}

		return true;
	}
	void GeneAlgorithm::run(std::vector<std::shared_ptr<Job>>& jobs) {
		
		check(jobs);
		exec();
		enqueJobs(jobs);
		step();
	}
	void GeneAlgorithm::Chromosome::enqueJob(std::shared_ptr<Job> job)
	{
		gens.emplace_back(job);
		auto gene_it = std::prev(gens.end());
		auto host_it = hosts.find(gene_it->host_);
		if (host_it == hosts.end())
		{
			host_it = hosts.emplace(gene_it->host_, gene_it).first;
		}
		auto& host_info = host_it->second;

		bool flag_min = host_info.make_span == min_span;
		bool flag_max = host_info.make_span == max_span;

		host_info.make_span += gene_it->expected_runtime;
		host_info.count += 1;

		if(flag_min)
			min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		if (flag_max)
			max_span = host_info.make_span;
	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for(auto job : jobs)
		{
			auto gene_it = gens.begin();
			while (gene_it != gens.end() && gene_it->job_ != job) ++gene_it;
			if (gene_it == gens.end())
			{
				error("fail to find gene");
			}

			auto host_it = hosts.find(gene_it->host_);

			{
				auto& host_info = host_it->second;
				bool max_flag = (max_span == host_info.make_span);
				bool min_flag = (min_span == host_info.make_span);

				host_info.make_span -= gene_it->expected_runtime;
				host_info.count -= 1;

				if (max_flag)
					max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
				if (min_flag)
					min_span = host_info.make_span;
			}

			if (host_it->second.count == 0)
			{
				hosts.erase(host_it);
			}
			else
			{
				auto next_gene = std::next(gene_it);
				while (next_gene != gens.end() && next_gene->host_!= host_it->first) ++next_gene;
				if (next_gene == gens.end())
				{
					error("fail to find next gen");
				}
				host_it->second.first_job_gene = next_gene;
			}
			gens.erase(gene_it);
		}

	}
	GeneAlgorithm::Chromosome::Gene::Gene(std::shared_ptr<Job> job): job_(job)
	{	
		const std::vector<Host>& all_hosts{ job->queue_managing_this_job->simulation_->get_cluster().vector() };
		//auto queue = job->queue_managing_this_job;
		//auto simulation = queue->simulation_;
		//auto& cluster = simulation->get_cluster();
		//auto& all_hosts = cluster.vector();
		size_t n = all_hosts.size();
		int i = rand() % n;
		while (all_hosts[i].max_slot < job->slot_required || all_hosts[i].max_mem < job->mem_required)
		{
			i = rand() % n;
		}
		host_ = (Host*)&all_hosts[i];
		expected_runtime = host_->get_expected_run_time(*job_);
	}
	GeneAlgorithm::Chromosome::HostInfo::HostInfo(std::list<Gene>::iterator& first_job_gene) :
		first_job_gene(first_job_gene),
		count(0) {};

}
