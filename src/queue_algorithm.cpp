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
		for (auto job : jobs)
		{
			if (job->state == JobState::WAIT)
				if (!run_job(job))
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
		if (length > 0)
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
			if (hosts.find(h_ptr) == hosts.end())
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
			map<Host*, Chromosome::HostInfo>& hosts{ best_chromosome.hosts };
			vector<shared_ptr<Job>> excuted_job;

			for (auto host_it : hosts)
			{
				Host* host{ host_it.first };
				Chromosome::HostInfo& host_info{ host_it.second };
				if (host_info.size > 0)
				{
					shared_ptr<Job> job = host_info.head->host_next->job_;

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
	void GeneAlgorithm::Chromosome::push_front(shared_ptr<Gene> ptr)
	{
		head->insert_back(ptr);
		++size;
	}
	void GeneAlgorithm::Chromosome::push_back(shared_ptr<Gene> ptr)
	{
		tail->insert_front(ptr);
		++size;
	}
	void GeneAlgorithm::Chromosome::detach(shared_ptr<Gene> gene)
	{
		
		job_map.erase(gene->job_);
		gene->pre->next = gene->next;
		gene->next->pre = gene->pre;
		gene->next = gene->pre = nullptr;
		hosts[gene->host_].detach(gene);
		--size;
		
	}
	void GeneAlgorithm::Chromosome::enqueJob(std::shared_ptr<Job> job)
	{
		shared_ptr<Gene> gene(new Gene(job));
		insert(gene);
		Host* host = gene->host_;
		auto& host_info = hosts[host];
		host_info.insert(gene);
		update_span();
	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for (auto job : jobs)
		{
			shared_ptr<Gene> gene = this->find(job);
			detach(gene);
		}
		update_span();
	}
	GeneAlgorithm::Chromosome GeneAlgorithm::Chromosome::mutation() const
	{
		Chromosome c(*this);
		c.mutate();
		return c;
	}
	void GeneAlgorithm::Chromosome::mutate()
	{
		size_t seed = pseudo_random(size_t(time(NULL)) ^ size_t(this) ^ size_t(clock()));
		vector<size_t> index;
		while (index.size() < std::min(size_t(MUTATION_GENE), size))
		{
			seed = pseudo_random(seed);
			size_t i = seed % size;
			if (std::find(index.begin(), index.end(), i) == index.end())
				index.push_back(i);
		}
		size_t i = 0;
		shared_ptr<Gene> cur = head->next;
		std::sort(index.begin(), index.end());
		for (size_t ii : index)
		{
			while (i < ii)
			{
				++i;
				cur = cur->next;
			}
			HostInfo& host_info{ hosts[cur->host_] };
			host_info.detach(cur);

			Host* host = cur->setRandomHost();

			auto it = hosts.insert({ host,HostInfo() }).first;
			HostInfo& new_host_info{ it->second };
			host_info.insert(cur);
			++i;
			cur = cur->next;
		}

		update_span();

	}
	void GeneAlgorithm::Chromosome::update_span()
	{
		max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
	}
	GeneAlgorithm::Chromosome::Chromosome(const Chromosome& ref)
	{
		head = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		tail = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		head->next = tail;
		tail->pre = head;

		max_span = ref.max_span;
		min_span = ref.min_span;
		auto cur = ref.head->next;
		while (cur != ref.tail)
		{
			auto gene = shared_ptr<Gene>(new Gene(*cur));
			insert(gene);
			auto it = hosts.insert({ gene->host_, HostInfo() }).first;
			it->second.insert(gene);
			cur = cur->next;
		}

	}
		GeneAlgorithm::Chromosome::Gene::Gene(std::shared_ptr<Job> job) : job_(job)
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
		size_t seed = pseudo_random((size_t)(time(NULL)) ^ (size_t)(this) ^ (size_t)(clock()) ^ (size_t)(&*job_));
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
	void GeneAlgorithm::Chromosome::Gene::insert_front(shared_ptr<Gene> ptr)
	{
		ptr->pre = pre;
		ptr->next = pre->next;
		pre->next = ptr;
		pre = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_back(shared_ptr<Gene> ptr)
	{
		ptr->next = next;
		ptr->pre = next->pre;
		next->pre = ptr;
		next = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_host_front(shared_ptr<Gene> ptr)
	{
		ptr->host_pre = host_pre;
		ptr->host_next = host_pre->host_next;
		host_pre->next = ptr;
		host_pre = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_host_back(shared_ptr<Gene> ptr)
	{
		ptr->host_next = host_next;
		ptr->host_pre = host_next->host_pre;
		host_next->host_pre = ptr;
		host_next = ptr;
	}


	void GeneAlgorithm::Chromosome::insert(shared_ptr<GeneAlgorithm::Chromosome::Gene> gene)
	{
		job_map[gene->job_] = gene;
		auto cur = tail->pre;
		while (cur != head && cur->job_->submit_time < gene->job_->submit_time) cur = cur->pre;
		cur->insert_back(gene);
		++size;
	}

	shared_ptr<GeneAlgorithm::Chromosome::Gene> GeneAlgorithm::Chromosome::find(shared_ptr<Job> job)
	{
		return job_map[job];
	}

	void GeneAlgorithm::Chromosome::HostInfo::insert(shared_ptr<Gene> gene)
	{
		auto cur = tail->host_pre;
		while (cur != head && cur->job_->submit_time < gene->job_->submit_time) cur = cur->host_pre;
		cur->insert_host_back(gene);
		++size;
		make_span += gene->expected_runtime;
	}

	void GeneAlgorithm::Chromosome::HostInfo::detach(shared_ptr<Gene> gene)
	{
		gene->host_pre->host_next = gene->host_next;
		gene->host_next->host_pre = gene->host_pre;
		gene->host_next = gene->host_pre = nullptr;
		make_span -= gene->expected_runtime;
		--size;
	}

	GeneAlgorithm::Chromosome::Chromosome()
	{
		head = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		tail = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		head->next = tail;
		tail->pre = head;
	}

	GeneAlgorithm::Chromosome::HostInfo::HostInfo()
	{
		head = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		tail = shared_ptr<GeneAlgorithm::Chromosome::Gene>(new GeneAlgorithm::Chromosome::Gene());
		head->host_next = tail;
		tail->host_pre = head;
	}
}
