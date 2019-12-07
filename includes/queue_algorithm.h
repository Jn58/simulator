#pragma once
#include <utility>
#include <functional>
#include <algorithm>
#include "queue.h"
#include "host.h"
#include <list>

#define POPULATION_SIZE 10
namespace ClusterSimulator
{
	// TODO: As Template
	class QueueAlgorithm
	{
	public:
		virtual const std::string& get_name() const noexcept = 0;
		
		virtual void run(std::vector<std::shared_ptr<Job>>& jobs) = 0;
	};

	/**
	 * Implementation of OLB Algorithm
	 */
	class OLBAlgorithm : public QueueAlgorithm
	{
		inline static const std::string name{"OLB"};
	public:
		const std::string& get_name() const noexcept override { return name; }

		void run(std::vector<std::shared_ptr<Job>>& jobs) override
		{
			for (auto& job : jobs)
			{
				auto hosts = job->get_eligible_hosts();
				if (hosts.empty()) continue;
				auto best_host = *std::min_element(hosts.begin(), hosts.end(), 
					[](const Host* a, const Host* b)
					{
						return a->remaining_slots() > b->remaining_slots();
					});
				best_host->execute_job(*job);
			}
		}
	};

	class GeneAlgorithm : public QueueAlgorithm
	{
		private:
		inline static const std::string name{"Genetic Algorithm"};
	public:
		const std::string& get_name() const noexcept override { return name; }
		class Chromosome
		{
		public:
			class Gene
			{
			public:
				std::shared_ptr<Job>& job_;
				std::shared_ptr<Host> host_;
				std::chrono::milliseconds expected_runtime;
				
				Gene(std::shared_ptr<Job>& job) : job_(job) {};
				~Gene() = default;

				std::chrono::milliseconds setHost(std::shared_ptr<Host>& host) {};
			};

			class HostInfo
			{
			public:
				int count = 0;
				std::chrono::milliseconds make_span = std::chrono::milliseconds(0);
				std::list<Gene>::iterator first_job;
			};

			std::list<Gene> gens;
			std::map<std::shared_ptr<Host>, HostInfo> hosts;

			std::chrono::milliseconds max_span = std::chrono::milliseconds(0);
			std::chrono::milliseconds min_span = std::chrono::milliseconds(0);

			size_t len(){ return gens.size(); };

			void enqueJob(std::shared_ptr<Job>& job) {};
			void deleteJobs(std::shared_ptr<Job>& job) {};

			Chromosome mutation() const{
				return Chromosome();
			};

			Chromosome crossOver(const Chromosome& other) const {
				return Chromosome();
			};

		};

		size_t length = 0;

		std::vector<Chromosome> population{POPULATION_SIZE};

		void enqueJobs(std::vector<std::shared_ptr<Job>>& jobs) {};
		void deleteJobs(std::vector<std::shared_ptr<Job>>& jobs) {};
		void update_pedding_job() {};
		GeneAlgorithm::Chromosome& getBestChromosome();
		void mutation() {};
		void crossOver() {};
		void sort() {};
		void dropout() {};

		void step() {};

		void exec();
		bool check(std::vector<std::shared_ptr<Job>>& jobs);

	
		void run(std::vector<std::shared_ptr<Job>>& jobs) override;

	};

	class QueueAlgorithms
	{
	public:
		inline static QueueAlgorithm* const OLB = new OLBAlgorithm();
		inline static QueueAlgorithm* const Genetic = new GeneAlgorithm();
	};
}



