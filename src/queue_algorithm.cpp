#include "queue_algorithm.h"
#include <cassert>
#include <cluster_simulation.h>
#include <algorithm>
#include <iterator>
#include "error.h"
#include <omp.h>

unsigned int pseudo_random(unsigned int x)
{
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
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
	bool GeneAlgorithm::run_job(std::shared_ptr<Job> job)
{
		std::vector<Host*> eligible_hosts{ job->get_eligible_hosts() };
		std::vector<Host*> idel_hosts;
		if (eligible_hosts.empty())
		{
			return false;
		}
		for (auto h_ptr : eligible_hosts)
		{
			auto& hosts = this->getBestChromosome().hosts;
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
			[=](const Host* a, const Host* b)
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

			std::vector<std::map<Host*, Chromosome::HostInfo>::iterator> hosts;

			for (auto it = best_hosts.begin(); it != best_hosts.end(); ++it)
			{
				hosts.push_back(it);
			}

			int thread_num = std::min(omp_get_max_threads(),int(hosts.size()));
			//thread_num = 1;
			std::vector <std::vector<std::shared_ptr<Job>>> excuted_jobs_vec(thread_num);
			std::vector<std::shared_ptr<Job>> excuted_jobs;
			int i;
#pragma omp parallel num_threads(thread_num)
			{
#pragma omp for
				for (i = 0; i < best_hosts.size(); ++i)
				{
					auto tid = omp_get_thread_num();

					Host* host{ hosts[i]->first };
					Chromosome::HostInfo& host_info{ hosts[i]->second };
					if (host_info.queue.size() == 0)
					{
						error("host info queue length must bigger than 0");
					}
					if (host != hosts[i]->second.queue.front()->host_)
					{
						error("host and first job's host must be same");
					}
					std::shared_ptr<Job> job{ host_info.queue.front()->job_ };
					if (host->is_executable(*job))
					{
#pragma omp critical
						host->execute_job(*job);
						excuted_jobs_vec[tid].push_back(job);
					}
				}
			}

			for (auto& vec : excuted_jobs_vec)
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
	bool GeneAlgorithm::check(std::vector<std::shared_ptr<Job>>& jobs)
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
			error("pedding_n and length must be same");
		}

		if (pedding_n + count != jobs.size())
		{
			for (auto job : jobs)
			{
				if (job->state != JobState::PEND && job->state != JobState::WAIT)
					printf("%d\n", job->state);
			}
			error("wrong state");
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
		gens.emplace_back(job);
		auto gene_it = std::prev(gens.end());
		auto host_it = hosts.find(gene_it->host_);
		if (host_it == hosts.end())
		{
			host_it = hosts.emplace(gene_it->host_, HostInfo{}).first;
		}
		auto& host_info = host_it->second;

		bool flag_min = host_info.make_span == min_span;
		bool flag_max = host_info.make_span == max_span;

		host_info.make_span += gene_it->expected_runtime;
		host_info.queue.push_back(gene_it);
		host_info.sort();
		if(flag_min)
			min_span = std::min_element(hosts.begin(), hosts.end(), [](const auto& a, const auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
		if (flag_max)
			max_span = host_info.make_span;
	}
	void GeneAlgorithm::Chromosome::chromosomeDeleteJobs(std::vector<std::shared_ptr<Job>>& jobs)
	{

		int thread_num = std::min(omp_get_max_threads(),int(hosts.size()));
		int i;
#pragma omp parallel num_threads(thread_num)
		{
#pragma omp for
			for(i = 0; i < jobs.size(); ++i)
			{
				auto& job{ jobs[i] };
				auto gene_it = gens.begin();
				while (gene_it != gens.end() && gene_it->job_ != job) gene_it = next(gene_it);
				if (gene_it == gens.end())
				{
					error("fail to find gene");
				}

				auto host{ gene_it->host_ };
				auto host_it = hosts.find(host);
				if (host_it == hosts.end())
				{
					error("can not fild host info for gene");
				}

				
				auto& host_info = host_it->second;
#pragma omp critical
				{
					bool max_flag = (max_span == host_info.make_span);
					bool min_flag = (min_span == host_info.make_span);
					host_info.make_span -= gene_it->expected_runtime;
					if (max_flag)
						max_span = std::max_element(hosts.begin(), hosts.end(), [](auto& a, auto& b) {return a.second.make_span < b.second.make_span; })->second.make_span;
					if (min_flag)
						min_span = host_info.make_span;
				}
				auto it = host_info.queue.begin();
				while (it != host_info.queue.end() && (*it)->job_ != job) ++it;
				if (it == host_info.queue.end())
				{
					error("can not find job in host info queue");
				}
#pragma omp critical
				{
					host_info.queue.erase(it);
					if (host_info.queue.size()== 0)
					{
						hosts.erase(host_it);
					}
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
		unsigned int seed = pseudo_random((unsigned int)(time(NULL)) ^ (unsigned int)(this));
		std::vector<unsigned int> index;
		int count = 0;
		while (index.size() < std::min(3,int(gens.size())))
		{
			seed = pseudo_random(seed);
			unsigned int i = seed % gens.size();
			if (index.end() == std::find(index.begin(), index.end(), i))
				index.push_back(i);
			++count;
			if (count > 1000)
			{
				error("something wrong with index");
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
				auto host_it = hosts.find(host);
				auto& host_info = host_it->second;
				bool max_flag = host_info.make_span == max_span;
				bool min_flag = host_info.make_span == min_span;
				host_info.make_span -= gene_it->expected_runtime;
				if (max_flag)
					max_span = std::max_element(hosts.begin(), hosts.end(), [](const auto& left, const auto& right) {return left.second.make_span < right.second.make_span; })->second.make_span;
				if (min_flag)
					min_span = host_info.make_span;
				host_info.queue.remove(gene_it);
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
				if (max_flag)
					max_span = host_info.make_span;
				if (min_flag)
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
			hosts[host_it.first].make_span = ref.hosts[host_it.first].make_span;
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
		const std::vector<Host>& all_hosts{ job_->queue_managing_this_job->simulation_->get_cluster().vector() };
		size_t n = all_hosts.size();
		int i = rand() % n;
		while (all_hosts[i].max_slot < job_->slot_required || all_hosts[i].max_mem < job_->mem_required)
		{
			i = rand() % n;
		}
		host_ = (Host*)&all_hosts[i]; 
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
