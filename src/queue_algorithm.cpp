#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>
#include "error_.h"
#include <omp.h>
#include <set>
#include <mutex>




size_t pseudo_random(size_t x)
{
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 17;
	return x;
}


namespace ClusterSimulator {
	using namespace std;
	unique_ptr<mutex> GeneAlgorithm::Chromosome::pool_mutex;

	vector<GeneAlgorithm::Chromosome::Gene*> GeneAlgorithm::Chromosome::genePool;

	void GeneAlgorithm::enqueJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		if (jobs.size()-length >  0)
		{
			{
				//int num_threads = std::min(jobs.size(), size_t(omp_get_max_threads()));
				int num_threads = omp_get_max_threads();
				//int num_threads = 1;
				vector<vector<shared_ptr<Job>>> queued_jobs_vec(num_threads);
				vector<shared_ptr<Job>> queued_jobs;
#pragma omp parallel num_threads(num_threads)
				{
#pragma omp for
					for (int i = 0; i < jobs.size()-length; ++i)
					{
						if (jobs[i]->state == JobState::WAIT)
						{
							int tid = omp_get_thread_num();
							shared_ptr<Job>& job = jobs[i];
							std::vector<Host*> eligible_hosts{ job->get_eligible_hosts() };
							std::vector<Host*> idel_hosts;
							unordered_map<Host*, GeneAlgorithm::Chromosome::HostInfo>& hosts = getBestChromosome().hosts;

							if (eligible_hosts.empty())
							{
								queued_jobs_vec[tid].push_back(job);
								continue;
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
								queued_jobs_vec[tid].push_back(job);
								continue;
							}

							Host* best_host = *std::min_element(idel_hosts.begin(), idel_hosts.end(),
								[&job](const Host* a, const Host* b)
								{
									return a->get_expected_run_time(*job) < b->get_expected_run_time(*job);
								});
#pragma omp critical
							{
								best_host->execute_job(*job);
							}
						}
						
					}
#pragma omp barrier
#pragma omp single
					{
						for (auto& vec : queued_jobs_vec)
						{
							queued_jobs.insert(queued_jobs.end(), vec.begin(), vec.end());
						}
						length += queued_jobs.size();
					}
#pragma omp for
					for (int i = 0; i < POPULATION_SIZE; ++i)
					{
						population[i]->enqueJobs(queued_jobs);
					}
				}
			}
		}
	}
	void GeneAlgorithm::deleteJobs(std::vector<std::shared_ptr<Job>>& executed_jobs)
	{
		{
			for (auto& p : population)
			{
				p->chromosomeDeleteJobs(executed_jobs);
			}
			length -= executed_jobs.size();
		}
	}
	GeneAlgorithm::Chromosome& GeneAlgorithm::getBestChromosome()
	{
		return *population[0];
	}
	void GeneAlgorithm::mutation()
	{
		if (length > 0)
			for (int i = 0; i < POPULATION_SIZE; ++i)
			{
				for (int ii = 0; ii < MUTAION_COUNT; ++ii)
				{
					Chromosome* c = population[i]->mutation();
					population.push_back(c);
				}
			}
	}

	void GeneAlgorithm::crossOver()
	{
		size_t seed = pseudo_random(size_t(rand()) ^ size_t(this) ^ size_t(clock()));
		for (int i = 0; i < CROSS_OVER; ++i)
		{
			seed = pseudo_random(seed);
			size_t l_index = seed % POPULATION_SIZE;
			seed = pseudo_random(seed);
			size_t r_index = seed % POPULATION_SIZE;
			while (l_index == r_index)
			{
				seed = pseudo_random(seed);
				r_index = seed % POPULATION_SIZE;
			}
			population.push_back(population[l_index]->crossOver(population[r_index]));
		}


	}

	void GeneAlgorithm::sort()
	{
		std::sort(population.begin(), population.end(), [](const Chromosome* left, const Chromosome* right) {
			if (left->max_span < right->max_span) return true;
			if (left->max_span != right->max_span) return false;
			if (left->hosts.size() > right->hosts.size()) return true;
			if (left->hosts.size() != right->hosts.size()) return false;
			if (left->min_span > right->min_span) return true;
			return false;
			});
	}
	void GeneAlgorithm::dropout()
	{
		while (population.size() > POPULATION_SIZE)
		{
			delete population.back();
			population.pop_back();
		}
	}
	GeneAlgorithm::GeneAlgorithm()
	{
		srand(time(NULL));
		Chromosome::pool_mutex = make_unique<mutex>();
		for (int i = 0; i < POPULATION_SIZE; ++i)
		{
			population.push_back(new Chromosome);
		}
	}
	GeneAlgorithm::~GeneAlgorithm()
	{
		while (population.size() > 0)
		{
			delete population.back();
			population.pop_back();
		}
	}
	bool GeneAlgorithm::run_job(std::shared_ptr<Job> job)
	{
		std::vector<Host*> eligible_hosts{ job->get_eligible_hosts() };
		std::vector<Host*> idel_hosts;
		unordered_map<Host*, GeneAlgorithm::Chromosome::HostInfo>& hosts = getBestChromosome().hosts;

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
			unordered_map<Host*, Chromosome::HostInfo>& hosts{ best_chromosome.hosts };
			vector<shared_ptr<Job>> excuted_job;

			for (auto& host_it : hosts)
			{
				Host* host{ host_it.first };
				Chromosome::HostInfo& host_info{ host_it.second };
				if (host_info.size > 0)
				{
					shared_ptr<Job>& job = host_info.host_head->host_next->job_;

					if (host->is_executable(*job))
					{
						host->execute_job(*job);
						excuted_job.push_back(job);
					}
				}
				else
				{
					error_("host_info must bigger than 0");
				}
			}
			deleteJobs(excuted_job);
		}
	}
	bool GeneAlgorithm::check(const std::vector<std::shared_ptr<Job>>& jobs) const
	{
		size_t count = 0;
		size_t pedding_n = 0;
		for (auto& job : jobs)
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
			for (auto& job : jobs)
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
	void GeneAlgorithm::Chromosome::push_front(Gene* ptr)
	{
		head->insert_back(ptr);
		++size;
	}
	void GeneAlgorithm::Chromosome::push_back(Gene* ptr)
	{
		tail->insert_front(ptr);
		++size;
	}
	void GeneAlgorithm::Chromosome::detach(Gene* gene)
	{
		
		//job_map.erase(gene->job_);
		gene->pre->next = gene->next;
		gene->next->pre = gene->pre;
		gene->next = gene->pre = nullptr;
		Host* host = gene->host_;
		HostInfo& host_info = hosts[host];
		host_info.detach(gene);
		if (host_info.size == 0)
		{
			hosts.erase(host);
		}
		job_map.erase(gene->job_);
		freeGene(gene);
		--size;
		
	}
	void GeneAlgorithm::Chromosome::enqueJob(std::shared_ptr<Job> job)
	{
		Gene* gene = allocGene(job);
		insert(gene);
		Host* host = gene->host_;
		HostInfo& host_info = hosts[host];
		host_info.insert(gene);
		update_span();
	}
	void GeneAlgorithm::Chromosome::enqueJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		if (jobs.size() > 0 )
		{
			for (auto& job : jobs)
			{
				Gene* gene = allocGene(job);
				insert(gene);
				Host* host = gene->host_;
				HostInfo& host_info = hosts[host];
				host_info.insert(gene);
			}
		}
		update_span();
	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{
		for (auto& job : jobs)
		{
			Gene* gene = this->find(job);
			detach(gene);
		}
		update_span();
	}
	GeneAlgorithm::Chromosome* GeneAlgorithm::Chromosome::mutation() const
	{
		Chromosome* c = new Chromosome(*this);
		c->mutate();
		return c;
	}
	void GeneAlgorithm::Chromosome::mutate()
	{
		size_t seed = pseudo_random(size_t(rand()) ^ size_t(this) ^ size_t(clock()));
		set<size_t> index;
		while (index.size() < size_t(MUTATION_GENE*size) )
		{
			seed = pseudo_random(seed);
			size_t i = seed % size;
			index.insert(i);
		}
		size_t i = 0;
		Gene* cur = head->next;
		for (size_t ii : index)
		{
			while (i < ii)
			{
				++i;
				cur = cur->next;
			}
			Host* host = cur->host_;
			HostInfo& host_info{ hosts[host] };
			host_info.detach(cur);
			if (host_info.size == 0)
			{
				hosts.erase(host);
			}

			host = cur->setRandomHost();

			
			HostInfo& new_host_info = hosts[host];
			new_host_info.insert(cur);
			++i;
			cur = cur->next;
		}

		update_span();

	}
	void GeneAlgorithm::Chromosome::update_span()
	{
		if (size > 0)
		{
			max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
			min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		}
		else
		{
			max_span = min_span = chrono::milliseconds(0);
		}
	}

	GeneAlgorithm::Chromosome::Chromosome(const Chromosome& ref) : head(allocGene()), tail(allocGene())
	{
		head->next = tail;
		tail->pre = head;

		max_span = ref.max_span;
		min_span = ref.min_span;
		Gene* cur = ref.head->next;
		for (int i = 0; i < ref.size; ++i)
		{
			Gene* gene = allocGene(cur);
			insert(gene);
			Host* host = gene->host_;
			HostInfo& host_info = hosts[host];
			host_info.insert(gene);
			cur = cur->next;
		}
		if (cur != ref.tail)
		{
			error_("??");
		}

	}
	GeneAlgorithm::Chromosome::~Chromosome()
	{
		Gene* cur = tail->pre;
		while (cur != head)
		{
			freeGene(cur->next);
			cur = cur->pre;
		}
		freeGene(cur);
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
	GeneAlgorithm::Chromosome::Gene::~Gene()
	{
	}
	GeneAlgorithm::Chromosome::Gene& GeneAlgorithm::Chromosome::Gene::operator=(const Gene& ref)
	{
		host_ = ref.host_;
		job_ = ref.job_;
		expected_runtime = ref.expected_runtime;
		return *this;
	}
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
		expected_runtime = host_->get_expected_run_time(*job_)* job_->slot_required / host_->max_slot;
		return host_;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_front(Gene* ptr)
	{
		ptr->pre = pre;
		ptr->next = this;
		pre->next = ptr;
		pre = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_back(Gene* ptr)
	{
		ptr->next = next;
		ptr->pre = this;
		next->pre = ptr;
		next = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_host_front(Gene* ptr)
	{
		ptr->host_pre = host_pre;
		ptr->host_next = this;
		host_pre->next = ptr;
		host_pre = ptr;
	}
	void GeneAlgorithm::Chromosome::Gene::insert_host_back(Gene* ptr)
	{
		ptr->host_next = host_next;
		ptr->host_pre = this;
		host_next->host_pre = ptr;
		host_next = ptr;
	}


	void GeneAlgorithm::Chromosome::insert(Gene* gene)
	{
		//job_map[gene->job_] = gene;
		Gene* cur = tail->pre;
		while (cur != head && cur->job_->submit_time > gene->job_->submit_time) 
			cur = cur->pre;
		cur->insert_back(gene);
		job_map[gene->job_] = gene;
		++size;
	}

	GeneAlgorithm::Chromosome::Gene* GeneAlgorithm::Chromosome::find(shared_ptr<Job>& job)
	{
		return job_map[job];
		//return job_map[job];
	}

	GeneAlgorithm::Chromosome::Gene* GeneAlgorithm::Chromosome::allocGene()
	{
		
		pool_mutex->lock();
		Gene* gene;
		if (genePool.size() > 0)
		{
			gene = genePool.back();
			genePool.pop_back();
		}
		else
		{
			gene = new Gene;
		}
		pool_mutex->unlock();
		return gene;
	}

	GeneAlgorithm::Chromosome::Gene* GeneAlgorithm::Chromosome::allocGene(shared_ptr<Job>& job)
	{
		Gene* gene = allocGene();
		gene->job_ = job;
		gene->setRandomHost();
		return gene;
	}

	GeneAlgorithm::Chromosome::Gene* GeneAlgorithm::Chromosome::allocGene(const Gene* other)
	{
		Gene* gene = allocGene();
		*gene = *other;
		return gene;
	}

	void GeneAlgorithm::Chromosome::freeGene(Gene* gene)
	{
		gene->pre = gene->next = gene->host_next = gene->host_pre = nullptr;
		gene->host_ = nullptr;
		gene->job_ = nullptr;

		pool_mutex->lock();
		genePool.push_back(gene);
		pool_mutex->unlock();
	}

	void GeneAlgorithm::Chromosome::HostInfo::insert(Gene* gene)
	{
		Gene* cur = host_tail->host_pre;
		while (cur != host_head && cur->job_->submit_time > gene->job_->submit_time) cur = cur->host_pre;
		cur->insert_host_back(gene);
		++size;
		make_span += gene->expected_runtime;
	}

	void GeneAlgorithm::Chromosome::HostInfo::detach(Gene* gene)
	{
		assert(gene->host_pre->host_next == gene);
		assert(gene->host_next->host_pre == gene);
		gene->host_pre->host_next = gene->host_next;
		gene->host_next->host_pre = gene->host_pre;
		gene->host_next = gene->host_pre = nullptr;
		make_span -= gene->expected_runtime;
		--size;
	}

	GeneAlgorithm::Chromosome* GeneAlgorithm::Chromosome::crossOver(const Chromosome* other) const
	{
		size_t seed = pseudo_random(size_t(rand()) ^ size_t(this) ^ size_t(clock()));
		Chromosome* p = new Chromosome;
		Gene* g1 = head->next;
		Gene* g2 = other->head->next;
		while (g1 != tail)
		{
			seed = pseudo_random(seed);
			Gene* g = (seed % 2) ? allocGene(g1) : allocGene(g2);
			p->insert(g);
			p->hosts[g->host_].insert(g);
			g1 = g1->next;
			g2 = g2->next;
		}
		p->update_span();
		return p;
	}

	GeneAlgorithm::Chromosome::Chromosome() : head(allocGene()), tail(allocGene())
	{
		head->next = tail;
		tail->pre = head;
	}

	GeneAlgorithm::Chromosome::HostInfo::HostInfo() : host_head(allocGene()), host_tail(allocGene())
	{
		host_head->host_next = host_tail;
		host_tail->host_pre = host_head;
	}
	GeneAlgorithm::Chromosome::HostInfo::~HostInfo()
	{
		freeGene(host_head);
		freeGene(host_tail);
	}
}
