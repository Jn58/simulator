#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>


namespace ClusterSimulator {
	void GeneAlgorithm::enqueJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for (auto job: jobs)
		{
			for (auto p : population)
			{
				p.enqueJob(job);
			}

		}
		length += jobs.size();
	}
	void GeneAlgorithm::deleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for (auto& p : population)
		{
			p.deleteJobs(jobs);
		}
		length -= jobs.size();
	}
	GeneAlgorithm::Chromosome& GeneAlgorithm::getBestChromosome()
	{
		return population[0];
	}
	void GeneAlgorithm::exec()
	{
		if (length != 0)
		{
			auto& best_hosts = getBestChromosome().hosts;
			std::vector<std::shared_ptr<Job>> excuted_jobs;
			for (auto it = best_hosts.begin(); it != best_hosts.end(); ++it)
			{
				auto& host = it->first;
				auto& job = it->second.first_job->job_;
				if (host->is_executable(*job))
				{
					host->execute_job(*job);
					excuted_jobs.push_back(job);
				}
			}
			deleteJobs(excuted_jobs);
		}
	}
	bool GeneAlgorithm::check(std::vector<std::shared_ptr<Job>>& jobs)
	{
		if(jobs.size() < length ) return false;
		if(jobs[jobs.size() - length]->state != JobState::PEND  ) return false;

		return true;
	}
	void GeneAlgorithm::run(std::vector<std::shared_ptr<Job>>& jobs) {
		assert(check());
		exec();
		enqueJobs(jobs);
		step();
	}
	void GeneAlgorithm::Chromosome::enqueJob(std::shared_ptr<Job>& job)
	{
		gens.emplace_back(job);
		auto gene = std::prev(gens.end());
		auto& h = gene->host_;
		hosts[h].count += 1;
		if (hosts[h].count == 1)
		{
			hosts[h].first_job = gene;
			hosts[h].make_span = gene->expected_runtime;
			if (hosts.size() == 1)
			{
				min_span = max_span = hosts[h].make_span;
			}
		}
		else
		{
			if (hosts[h].make_span == max_span)
			{
				if (hosts[h].make_span == min_span)
					min_span = max_span = hosts[h].make_span += gene->expected_runtime;
				else
					max_span = hosts[h].make_span += gene->expected_runtime;
			}
			else if (hosts[h].make_span == min_span)
			{
				hosts[h].make_span += gene->expected_runtime;
				min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;

			}
		}
		
		
	}
	void GeneAlgorithm::Chromosome::deleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		auto gene_it = gens.begin();
		for (auto job_it = jobs.begin(); job_it != jobs.end(); ++job_it)
		{
			while (gene_it->job_ != *job_it) ++gene_it;
			gene_it = gens.erase(gene_it);
		}

	}
	GeneAlgorithm::Chromosome::Gene::Gene(std::shared_ptr<Job>& job): job_(job)
	{
		auto all_hosts = job->queue_managing_this_job->simulation_->get_cluster().vector();
		int n = all_hosts.size();
		int i = rand() % n;
		while (all_hosts[i].max_slot < job->slot_required || all_hosts[i].max_mem < job->mem_required)
		{
			i = rand() % n;
		}
		host_ = std::shared_ptr<Host>(&all_hosts[i]);
		expected_runtime = host_->get_expected_run_time(*job_);
	}
}
