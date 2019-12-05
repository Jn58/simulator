#include "queue_algorithm.h"
#include <cassert>


namespace ClusterSimulator {
	void GeneAlgorithm::run(std::vector<std::shared_ptr<Job>>& jobs) {
		assert(check());
		exec();
		enqueJobs(jobs);
		step();
	}
}
