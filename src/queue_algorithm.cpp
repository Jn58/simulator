#include "queue_algorithm.h"
#include <cassert>


namespace ClusterSimulator {
	void GeneAlgorithm::deleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for (auto& p : population)
		{
			p.deleteJobs(jobs);
		}
	}
	GeneAlgorithm::Chromosome& GeneAlgorithm::getBestChromosome()
	{
		return population[0];
	}
	void GeneAlgorithm::exec()
	{
		auto &best_hosts = getBestChromosome().hosts;
		std::vector<std::shared_ptr<Job>> excuted_jobs;
		for (auto it = best_hosts.begin(); it != best_hosts.end() ; ++it)
		{
			auto & host = it->first;
			auto& job = it->second.first_job->job_;
			if (host->is_executable(*job))
			{
				host->execute_job(*job);
				excuted_jobs.push_back(job);
			}
		}
		deleteJobs(excuted_jobs);
		
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
}
