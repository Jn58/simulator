#pragma once
#include <utility>
#include <functional>
#include <algorithm>
#include "queue.h"
#include "host.h"
#include <list>
#include <cstdlib>
#include <ctime>
#include <unordered_map>

#define POPULATION_SIZE 2
#define MUTAION_COUNT 1
#define MUTATION_GENE 0.1

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
				Gene* pre=nullptr;
				Gene* next=nullptr;
				Gene* host_pre=nullptr;
				Gene* host_next=nullptr;
				
				Gene(std::shared_ptr<Job> job);
				Gene(const Gene& ref);
				Gene() {};
				~Gene();

				Gene& operator=(const Gene&) = delete;

				//std::chrono::milliseconds setHost(Host* host) {};
				Host* setRandomHost();

				void insert_front(Gene* ptr);
				void insert_back(Gene* ptr);
				void insert_host_front(Gene* ptr);
				void insert_host_back(Gene* ptr);

			};

			class HostInfo
			{
			public:
				std::chrono::milliseconds make_span = std::chrono::milliseconds(0);
				Gene* const host_head;
				Gene* const host_tail;
				size_t size = 0;
				HostInfo(const HostInfo&) = delete;
				HostInfo& operator=(const HostInfo&) = delete;
				HostInfo();
				~HostInfo();
				void insert(Gene* gene);
				void detach(Gene* gene);
			};

			void push_front(Gene* ptr);
			void push_back(Gene* ptr);

			

			Gene* const head=nullptr;
			Gene* const tail=nullptr;
			size_t size = 0;

			//std::map<shared_ptr<Job>, Gene*> job_map;
			unordered_map<Host*, HostInfo> hosts;
			void detach(Gene* gene);

			std::chrono::milliseconds max_span = std::chrono::milliseconds(0);
			std::chrono::milliseconds min_span = std::chrono::milliseconds(0);

			size_t len(){ return size; };

			void enqueJob(std::shared_ptr<Job> job);
			void enqueJobs(std::vector<std::shared_ptr<Job>>& jobs);
			void chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs);

			Chromosome* mutation() const;

			void mutate();
			void update_span();

			Chromosome crossOver(const Chromosome& other) const {
				return Chromosome();
			};

			Chromosome();
			Chromosome(const Chromosome& ref);
			Chromosome& operator=(const Chromosome&) = delete;
			~Chromosome();
			void insert(Gene* gene);
			Gene* find(shared_ptr<Job>& job);

		};

		GeneAlgorithm();
		~GeneAlgorithm();

		size_t length = 0;

		vector<Chromosome*> population;

		//GeneAlgorithm() { srand(time(NULL)); }

		
		bool run_job(std::shared_ptr<Job> job);
		void enqueJobs(std::vector<std::shared_ptr<Job>>& jobs);
		void deleteJobs(std::vector<std::shared_ptr<Job>>& jobs);
		void update_pedding_job() {};
		GeneAlgorithm::Chromosome& getBestChromosome();
		void mutation();
		void crossOver() {};
		void sort();
		void dropout();
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



