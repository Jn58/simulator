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
		int i;
		for (i = 0; i < POPULATION_SIZE; ++i)
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
			std::map<Host*, Chromosome::HostInfo>& best_hosts{ best_chromosome.hosts };
			//std::vector<std::shared_ptr<Job>> excuted_jobs;

			std::vector<std::map<Host*, Chromosome::HostInfo>::iterator> hosts_vec;

			for (auto it = best_hosts.begin(); it != best_hosts.end(); ++it)
			{
				hosts_vec.push_back(it);
			}

			int thread_num = std::min(omp_get_max_threads(),int(hosts_vec.size()));
			thread_num = 1;
			//thread_num = 1;
			std::vector <std::vector<std::shared_ptr<Job>>> excuted_jobs_vec(thread_num);
			std::vector<std::shared_ptr<Job>> excuted_jobs;
#pragma omp parallel num_threads(thread_num)
			{
#pragma omp for
				for (int i = 0; i < hosts_vec.size(); ++i)
				{
					auto tid = omp_get_thread_num();

					Host* host{ hosts_vec[i]->first };
					Chromosome::HostInfo& host_info{ hosts_vec[i]->second };
					if (host_info.queue.size() == 0)
					{
						error_("host info queue length must bigger than 0");
					}
					auto& front = host_info.queue.front();
					if (host != front->host_)
					{
						error_("host and first job's host must be same");
					}
					std::shared_ptr<Job>& job{ front->job_ };
					if (host->is_executable(*job))
					{
#pragma omp critical
						{
							 host->execute_job(*job);
						}
						excuted_jobs_vec[tid].push_back(job);
					}
				}
			}
// end of parallel

			for (const auto& vec : excuted_jobs_vec)
			{
				excuted_jobs.insert(excuted_jobs.end(), vec.begin(), vec.end());
			}

		/*	for (i = 0; i < thread_num; ++i)
			{
				for (int j = 0; j < excuted_host[i].size(); ++j)
				{
					auto& job{ excuted_jobs_vec[i][j] };
					excuted_host[i][j]->execute_job(*job);
					excuted_jobs.push_back(job);

				}
			}*/

			if (!excuted_jobs.empty())
			{
				deleteJobs(excuted_jobs);
			}

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
		host_info.queue.push_back(gene_it);
		host_info.sort();

		min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		max_span = std::max_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;

	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{

		int thread_num = std::min(omp_get_max_threads(),int(jobs.size()));
		thread_num = 1;
		std::map<Host*, std::vector<std::list<Gene>::iterator>> assigned_hosts;
		std::vector<std::map<Host*, std::vector<std::list<Gene>::iterator>>::iterator> assigned_hosts_vec;
		int i;
		
		std::vector<std::vector<std::list<Gene>::iterator>> gene_its_vec(thread_num);

#pragma omp parallel num_threads(thread_num)
		{
#pragma omp for
			for (i = 0; i < jobs.size(); ++i)
			{
				auto job{ jobs[i] };
				int tid = omp_get_thread_num();
				auto gene_it = std::find(gens.begin(),gens.end(),job);
				if (gene_it == gens.end())
				{
					error_("fail to find gene");
				}
				gene_its_vec[tid].push_back(gene_it);
			}
		}
// end of parallel

		for (auto& vec : gene_its_vec)
		{
			for (auto it : vec)
			{
				assigned_hosts[it->host_].push_back(it);
			}
		}
		
		for(auto it = assigned_hosts.begin(); it != assigned_hosts.end(); ++it)
		{
			assigned_hosts_vec.push_back(it);
		}

		thread_num = std::min(omp_get_max_threads(),int(assigned_hosts_vec.size()));
		thread_num = 1;

#pragma omp parallel num_threads(thread_num)
		{
#pragma omp for
			for (i = 0; i < assigned_hosts_vec.size(); ++i)
			{
				Host* host{ assigned_hosts_vec[i]->first };
				auto& host_info{ hosts[host] };
				for (auto gene_it : assigned_hosts_vec[i]->second)
				{
					if (gene_it->host_ != host)
					{
						error_("gene host must same with host");
					}
					auto job{ gene_it->job_ };
					host_info.make_span -= gene_it->expected_runtime;
					int n = host_info.queue.size();
					auto it = std::find_if(host_info.queue.begin(), host_info.queue.end(), [&job](const auto& ref) {return job == ref->job_; });
					host_info.queue.erase(it);
					if (n == host_info.queue.size())
					{
						error_("failt to find gene iterator in queue");
					}
				}
				if (host_info.queue.size() == 0)
				{
#pragma omp critical
					{
						int n = hosts.erase(host);
						if (n == 0)
						{
							error_("fail to remove host info");
						}
					}
				}
			}
		}
// end of parallel

		max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		min_span = std::min_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;

		for (auto& vec : gene_its_vec)
		{
			for (auto& it : vec)
			{
				int n = gens.size();
				gens.erase(it);
				if (n == gens.size())
				{
					error_("fail to remove gene");
				}
			}
		}
	}



	GeneAlgorithm::Chromosome GeneAlgorithm::Chromosome::mutation() const
	{
		Chromosome c(*this);
		c.mutate();
		return  c;
	}
	void GeneAlgorithm::Chromosome::mutate()
	{
		size_t seed = pseudo_random((size_t)(time(NULL)) ^ (size_t)(this));
		std::vector<size_t> index;
		int count = 0;
		while (index.size() < std::min(MUTATION_GENE,int(gens.size())))
		{
			seed = pseudo_random(seed);
			size_t i = seed % gens.size();
			if (index.end() == std::find(index.begin(), index.end(), i))
				index.push_back(i);
			++count;
			if (count > 1000)
			{
				error_("something wrong with index");
			}
		}
		std::sort(index.begin(), index.end());
		int i = 0;
		auto gene_it = gens.begin();
		for (auto ii : index)
		{
			while (i < ii)
			{
				++i;
				++gene_it;
			}
			{
				Host* host = gene_it->host_;
				auto job = gene_it->job_;
				auto host_it = hosts.find(host);
				auto& host_info = host_it->second;
				host_info.make_span -= gene_it->expected_runtime;
				max_span = std::max_element(hosts.begin(), hosts.end(), [](const auto& left, const auto& right) {return left.second.make_span < right.second.make_span; })->second.make_span;
				min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& left, const auto& right) {return left.second.make_span < right.second.make_span; })->second.make_span;
				auto it = std::find_if(host_info.queue.begin(), host_info.queue.end(), [&job](const auto& ref) {return job == ref->job_; });
				host_info.queue.erase(it);
				if (host_info.queue.size() == 0)
				{
					hosts.erase(host_it);
				}
			}
			{
				Host* host = gene_it->setRandomHost();
				auto& host_info = hosts[host];
				bool max_flag = host_info.make_span == max_span;
				bool min_flag = host_info.make_span == min_span;
				host_info.make_span += gene_it->expected_runtime;
				max_span = std::max_element(hosts.begin(), hosts.end(), [](const auto& left, const auto& right) {return left.second.make_span < right.second.make_span; })->second.make_span;
				min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& left, const auto& right) {return left.second.make_span < right.second.make_span; })->second.make_span;
				host_info.queue.push_back(gene_it);
				host_info.sort();
			}

		}
	}
	GeneAlgorithm::Chromosome::Chromosome(const Chromosome& ref)
	{
		max_span = ref.max_span;
		min_span = ref.min_span;
		gens = ref.gens;
		for(auto it = gens.begin(); it!= gens.end() ; ++it)
		{
			hosts[it->host_].queue.push_back(it);
			hosts[it->host_].make_span += it->expected_runtime;
		}
		for (auto& host_it : hosts)
		{
			hosts[host_it.first].make_span = ref.hosts.at(host_it.first).make_span;
			host_it.second.sort();
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
	void GeneAlgorithm::Chromosome::HostInfo::sort()
	{
		queue.sort([](const auto& left, const auto& right) {
			return left->job_->submit_time < right->job_->submit_time;
			});
	}
}
