#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>
#include "error_.h"
#include <omp.h>

size_t pseudo_random(size_t x)
{
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return x;
}


namespace ClusterSimulator {
	using namespace std;
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
	void GeneAlgorithm::mutation()
	{
		if(length > 0)
			for (int i = 0; i < POPULATION_SIZE; ++i)
			{
				for (int ii = 0; ii < MUTAION_COUNT; ++ii)
				{
					Chromosome c = population[i].mutation();
					population.push_back(c);
				}
			}
	}

	void GeneAlgorithm::sort()
	{
		std::sort(population.begin(), population.end(), [](const Chromosome& left, const Chromosome& right) {
			if (left.max_span < right.max_span) return true;
			if (left.max_span != right.max_span) return false;
			if (left.hosts.size() > right.hosts.size()) return true;
			if (left.hosts.size() != right.hosts.size()) return false;
			if (left.min_span > right.min_span) return true;
			return false;
			});
	}
	bool GeneAlgorithm::run_job(std::shared_ptr<Job> job)
{
		std::vector<Host*> eligible_hosts{ job->get_eligible_hosts() };
		std::vector<Host*> idel_hosts;
		auto& hosts = getBestChromosome().hosts;

		if (eligible_hosts.empty())
		{
			return false;
		}

		for (Host* h_ptr : eligible_hosts)
		{
			if( hosts.find(h_ptr) == hosts.end())
			{
				idel_hosts.push_back(h_ptr);
			}
		}
		if (idel_hosts.empty())
		{
			return false;
		}

		Host* best_host = *std::min_element(idel_hosts.begin(), idel_hosts.end(), 
			[&job](const Host* a, const Host* b)
			{
				return a->get_expected_run_time(*job) < b->get_expected_run_time(*job);
			});
		best_host->execute_job(*job);
		return true;
	}
	void GeneAlgorithm::step()
	{
		mutation();
		crossOver();
		sort();
		dropout();
	}
	void GeneAlgorithm::exec()
	{

		if (length != 0)
		{
			Chromosome& best_chromosome{ getBestChromosome() };
			list<Chromosome::Gene>& gens{ best_chromosome.gens };
			map<Host*, Chromosome::HostInfo>& hosts{ best_chromosome.hosts };
			vector<shared_ptr<Job>> excuted_job;

			for (auto host_it : hosts)
			{
				Host* host{ host_it.first };
				Chromosome::HostInfo& host_info{ host_it.second };
				if (host_info.count > 0)
				{
					shared_ptr<Job> job = find_if(gens.begin(), gens.end(), [&host](const Chromosome::Gene& gene) {return gene.host_ == host; })->job_;

					if (host->is_executable(*job))
					{
						host->execute_job(*job);
						excuted_job.push_back(job);
					}
				}
			}
			deleteJobs(excuted_job);
		}
	}
	bool GeneAlgorithm::check(const std::vector<std::shared_ptr<Job>>& jobs) const
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
			printf("%d\n", this->length);
			error_("pedding_n and length must be same");
		}

		if (pedding_n + count != jobs.size())
		{
			for (auto job : jobs)
			{
				if (job->state != JobState::PEND && job->state != JobState::WAIT)
					printf("%d\n", job->state);
			}
			error_("wrong state");
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
		Gene gene(job);
		this->gens.push_back(gene);
		auto gene_it = std::prev(this->gens.end());
		Host* host = gene.host_;
		if (gene_it->host_ != host || gene_it->job_ != job)
		{
			error_("wrong iterator");
		}
		auto& host_info = hosts[host];


		host_info.make_span += gene_it->expected_runtime;
		++host_info.count;

		min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		max_span = std::max_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;

	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for(auto job : jobs)
		{
			list<Gene>::iterator gene_it = find_if(gens.begin(), gens.end(), [&job](const Gene& gene) {return gene.job_ == job; });
			if (gene_it == gens.end())
			{
				error_("fail to find gene");
			}
			Host* host{ gene_it->host_ };
			auto host_it = hosts.find(host);
			if (hosts.end() == host_it)
			{
				error_("fail to find host info");
			}

			auto& host_info = host_it->second;

			host_info.make_span -= gene_it->expected_runtime;
			--host_info.count;

			if (host_info.count == 0)
			{
				hosts.erase(host);
			}

			gens.erase(gene_it);

		}
	}
	GeneAlgorithm::Chromosome GeneAlgorithm::Chromosome::mutation() const
	{
		Chromosome c(*this);
		c.mutate();
		return c;
	}
	void GeneAlgorithm::Chromosome::mutate()
	{
		size_t seed = pseudo_random(size_t(time(NULL)) ^ size_t(this)^size_t(clock()));
		vector<size_t> index;
		while (index.size() < std::min(size_t(MUTATION_GENE),gens.size()))
		{
			seed = pseudo_random(seed);
			size_t i = seed % gens.size();
			if (find(index.begin(), index.end(), i) == index.end())
				index.push_back(i);
		}
		size_t i = 0;
		list<Gene>::iterator gene_it = gens.begin();
		std::sort(index.begin(), index.end());
		for (size_t ii : index)
		{
			while (i < ii)
			{
				++i;
				++gene_it;
			}
			Host* host{ gene_it->host_ };
			HostInfo& host_info{ hosts[host] };
			--host_info.count;
			host_info.make_span -= gene_it->expected_runtime;
			
			host = gene_it->setRandomHost();

			HostInfo& new_host_info{ hosts[host] };
			++new_host_info.count;
			new_host_info.make_span += gene_it->expected_runtime;
			++i;
			++gene_it;
		}

		max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
	}
	GeneAlgorithm::Chromosome::Chromosome(const Chromosome& ref)
	{
		max_span = ref.max_span;
		min_span = ref.min_span;
		gens = ref.gens;
		for(auto& gene : gens)
		{
			Host* host{ gene.host_ };
			shared_ptr<Job> job{ gene.job_ };
			auto it = hosts.find(host);
			if (it == hosts.end()) hosts[host] = HostInfo();
			auto& host_info{ hosts[host] };

			++host_info.count;
			host_info.make_span += gene.expected_runtime;
		}
	}
	GeneAlgorithm::Chromosome::Gene::Gene(std::shared_ptr<Job> job): job_(job)
	{	
		setRandomHost();
	}
	GeneAlgorithm::Chromosome::Gene::Gene(const Gene& ref) :
		job_(ref.job_),
		host_(ref.host_),
		expected_runtime(ref.expected_runtime)
	{}
	Host* GeneAlgorithm::Chromosome::Gene::setRandomHost()
	{
		std::vector<Host>& all_hosts{ job_->queue_managing_this_job->simulation_->get_cluster().vector() };
		size_t n = all_hosts.size();
		size_t seed = pseudo_random((size_t)(time(NULL)) ^ (size_t)(this)^(size_t)(clock())^(size_t)(&*job_));
		seed = pseudo_random(seed);
		int i = seed % n;
		while (all_hosts[i].max_slot < job_->slot_required || all_hosts[i].max_mem < job_->mem_required)
		{
			seed = pseudo_random(seed);
			i = seed % n;
		}
		host_ = &all_hosts[i]; 
		expected_runtime = host_->get_expected_run_time(*job_);
		return host_;
	}
}
