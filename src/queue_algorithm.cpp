#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>


namespace ClusterSimulator {
	void GeneAlgorithm::enqueJobs(std::vector<std::shared_ptr<Job>>* jobs)
	{
		for (int i = 0; i < jobs->size(); ++i)
		{
			auto job = &(*jobs->at(i));
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
	void GeneAlgorithm::deleteJobs(std::vector<Job*>* jobs)
	{
		if (jobs->size() != 0)
		{
			for (auto& p : population)
			{
				p.chromosomeDeleteJobs(jobs);
			}
			length -= jobs->size();
		}
	}
	GeneAlgorithm::Chromosome& GeneAlgorithm::getBestChromosome()
	{
		return population[0];
	}
	bool GeneAlgorithm::run_job(Job* job)
{
		auto hosts = job->get_eligible_hosts();
		if (hosts.empty())
		{
			return false;
		}

		auto best_host = *std::min_element(hosts.begin(), hosts.end(), 
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
			auto best_hosts = &(getBestChromosome().hosts);
			std::vector<Job*> excuted_jobs;
			for (auto it = best_hosts->begin(); it != best_hosts->end(); ++it)
			{
				auto host = it->first;
				auto host_info = &it->second;
				if(host_info->count>0)
				{
					if (host != host_info->first_job_gene->host_)
					{
						printf(";klasjdf\n");
					}
					auto job = host_info->first_job_gene->job_;
					if (host->is_executable(*job))
					{
						host->execute_job(*job);
						excuted_jobs.push_back(job);
					}
				}
			}
			deleteJobs(&excuted_jobs);
		}
	}
	bool GeneAlgorithm::check(std::vector<Job*>* jobs)
	{
		auto i = 0;
		for (; i < jobs->size() && jobs->at(i)->state == JobState::WAIT; ++i);
		auto pedd_count = jobs->size() - i;
		if (pedd_count != length)
		{
			printf("%d\t%d\n", pedd_count, length);
		}
		if (length == 0) return true;

		if(jobs->at(jobs->size() - length)->state != JobState::PEND  ) return false;

		return true;
	}
	void GeneAlgorithm::run(std::vector<std::shared_ptr<Job>>& jobs) {
		
		exec();
		enqueJobs(&jobs);
		step();
	}
	void GeneAlgorithm::Chromosome::enqueJob(Job* job)
	{
		gens.emplace_back(job);
		auto gene = std::prev(gens.end());
		auto host_it = hosts.find(gene->host_);
		if (host_it == hosts.end())
		{
			host_it = hosts.insert(std::pair(gene->host_, Chromosome::HostInfo(gene))).first;
		}
		auto host_info = &host_it->second;

		bool flag_min = host_info->make_span == min_span;
		bool flag_max = host_info->make_span == max_span;

		host_info->make_span += gene->expected_runtime;
		host_info->count += 1;

		if (host_info->count == 1)
		{
			host_info->first_job_gene = gene;
		}

		if(flag_min)
			min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		if (flag_max)
			max_span = host_info->make_span;
	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<Job*>* jobs)
	{
		for(int i = 0 ; i < jobs->size(); ++i)
		{
			auto gene_it = gens.begin();
			auto job = jobs->at(i);
			while (gene_it != gens.end() && gene_it->job_ != job) ++gene_it;
			if (gene_it == gens.end())
			{
				printf("%d\n", job->id);
				printf("asdfasdf\n");
			}
			auto host_it = hosts.find(gene_it->host_);
			auto host_info = &host_it->second;
			bool max_flag = (max_span == host_info->make_span);
			bool min_flag = (min_span == host_info->make_span);
			host_info->make_span -= gene_it->expected_runtime;
			host_info->count -= 1;

			if (max_flag)
				max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
			if (min_flag)
				min_span = host_info->make_span;

			if (host_info->count == 0)
			{
				//hosts.erase(host_it);
			}
			else
			{
				auto next_gene = gens.begin();
				while (next_gene != gens.end() && next_gene->job_ != job) ++next_gene;
				if (next_gene == gens.end())
				{
					printf("asdfasdf\n");
				}
				host_info->first_job_gene = next_gene;

			}
			gene_it = gens.erase(gene_it);
		}

	}
	GeneAlgorithm::Chromosome::Gene::Gene(Job* job): job_(job)
	{	
		std::vector<Host>* all_hosts{ &(job->queue_managing_this_job->simulation_->get_cluster().nodes_) };
		//auto queue = job->queue_managing_this_job;
		//auto simulation = queue->simulation_;
		//auto& cluster = simulation->get_cluster();
		//auto& all_hosts = cluster.vector();
		int n = all_hosts->size();
		int i = rand() % n;
		while (all_hosts->at(i).max_slot < job->slot_required || all_hosts->at(i).max_mem < job->mem_required)
		{
			i = rand() % n;
		}
		host_ = &all_hosts->at(i);
		expected_runtime = host_->get_expected_run_time(*job_);
	}
	GeneAlgorithm::Chromosome::HostInfo::HostInfo(std::list<Gene>::iterator first_job_gene) : first_job_gene(first_job_gene)
	{
	}
}
