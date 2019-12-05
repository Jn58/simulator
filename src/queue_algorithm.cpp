#include "queue_algorithm.h"
#include <cassert>


namespace ClusterSimulator {
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
