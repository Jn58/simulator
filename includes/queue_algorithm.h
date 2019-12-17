#pragma once
#include <utility>
#include <functional>
#include <algorithm>
#include "queue.h"
#include "host.h"
#include <list>
#include <cstdlib>
#include <ctime>

#define POPULATION_SIZE 2
#define MUTAION_COUNT 2
#define MUTATION_GENE 2

namespace ClusterSimulator
{
	using namespace std;
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
				std::shared_ptr<Job> job_ = nullptr;
				Host* host_ = nullptr;
				std::chrono::milliseconds expected_runtime = std::chrono::milliseconds(0);
				shared_ptr<Gene> pre=nullptr;
				shared_ptr<Gene> next=nullptr;
				shared_ptr<Gene> host_pre=nullptr;
				shared_ptr<Gene> host_next=nullptr;
				
				Gene(std::shared_ptr<Job> job);
				Gene(const Gene& ref);
				Gene() {};
				~Gene();

				std::chrono::milliseconds setHost(Host* host) {};
				Host* setRandomHost();

				void insert_front(shared_ptr<Gene> ptr);
				void insert_back(shared_ptr<Gene> ptr);
				void insert_host_front(shared_ptr<Gene> ptr);
				void insert_host_back(shared_ptr<Gene> ptr);

			};

			class HostInfo
			{
			public:
				std::chrono::milliseconds make_span = std::chrono::milliseconds(0);
				shared_ptr<Gene> host_head=nullptr;
				shared_ptr<Gene> host_tail=nullptr;
				size_t size = 0;
				HostInfo();
				void insert(shared_ptr<Gene> gene);
				void detach(shared_ptr<Gene> gene);
			};

			void push_front(shared_ptr<Gene> ptr);
			void push_back(shared_ptr<Gene> ptr);

			

			shared_ptr<Gene> head=nullptr;
			shared_ptr<Gene> tail=nullptr;
			size_t size = 0;

			//std::map<shared_ptr<Job>, shared_ptr<Gene>> job_map;
			std::map<Host*, HostInfo> hosts;
			void detach(shared_ptr<Gene> gene);

			std::chrono::milliseconds max_span = std::chrono::milliseconds(0);
			std::chrono::milliseconds min_span = std::chrono::milliseconds(0);

			size_t len(){ return size; };

			void enqueJob(std::shared_ptr<Job> job);
			void enqueJobs(std::vector<std::shared_ptr<Job>>& jobs);
			void chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs);

			Chromosome mutation() const;

			void mutate();
			void update_span();

			Chromosome crossOver(const Chromosome& other) const {
				return Chromosome();
			};

			Chromosome();
			Chromosome(const Chromosome& ref);
			void insert(shared_ptr<Gene> gene);
			shared_ptr<Gene> find(shared_ptr<Job> job);

		};

		size_t length = 0;

		std::vector<Chromosome> population{POPULATION_SIZE};

		//GeneAlgorithm() { srand(time(NULL)); }
		GeneAlgorithm() { srand(0); }

		
		bool run_job(std::shared_ptr<Job> job);
		void enqueJobs(std::vector<std::shared_ptr<Job>>& jobs);
		void deleteJobs(std::vector<std::shared_ptr<Job>>& jobs);
		void update_pedding_job() {};
		GeneAlgorithm::Chromosome& getBestChromosome();
		void mutation();
		void crossOver() {};
		void sort();
		void dropout() {
			population.erase(population.begin() + POPULATION_SIZE, population.end());
		};

		void step();

		void exec();
		bool check(const std::vector<std::shared_ptr<Job>>& jobs) const;

	
		void run(std::vector<std::shared_ptr<Job>>& jobs) override;

	};

	class QueueAlgorithms
	{
	public:
		inline static QueueAlgorithm* const OLB = new OLBAlgorithm();
		inline static QueueAlgorithm* const Genetic = new GeneAlgorithm();
	};
}



