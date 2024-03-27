#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>

#include <set>

#include "util.hpp"
#include "index_set.hpp"
#include "jobs.hpp"
#include "cache.hpp"
#include "statistics.hpp"

namespace NP {

	namespace Global {

	  // RV: Job_index is defined in jobs.hpp in the NP namespace.
	  // typedef std::size_t Job_index;
	  typedef Index_set Job_set;
	  
		typedef std::vector<Job_index> Job_precedence_set;

		template<class Time> class Schedule_node;

		template<class Time> class Schedule_state
		{
			private:

		  // RV: would another data structure improve access to job_finish_times?
		  //     use vector, as in uni/state.hpp
		  //typedef typename std::unordered_map<Job_index, Interval<Time>> JobFinishTimes;
		  typedef typename std::vector<std::pair<Job_index, Interval<Time>>> JobFinishTimes;
			JobFinishTimes job_finish_times;

			public:

			// initial state -- nothing yet has finished, nothing is running
			Schedule_state(unsigned int num_processors)
			: scheduled_jobs()
			, num_jobs_scheduled(0)
			, core_avail{num_processors, Interval<Time>(Time(0), Time(0))}
			, lookup_key{0x9a9a9a9a9a9a9a9aUL}
			{
				assert(core_avail.size() > 0);
			}

			// transition: new state by scheduling a job in an existing state,
			//             by replacing a given running job.
			Schedule_state(
				const Schedule_state& from,
				Job_index j,
				const Job_precedence_set& predecessors,
				Interval<Time> start_times,
				Interval<Time> finish_times,
				hash_value_t key)
			: num_jobs_scheduled(from.num_jobs_scheduled + 1)
			, scheduled_jobs{from.scheduled_jobs, j}
			, lookup_key{from.lookup_key ^ key}
			, job_finish_times{from.job_finish_times}
			{
			  // RV: check what changed in this function.
				auto est = start_times.min();
				auto lst = start_times.max();
				auto eft = finish_times.min();
				auto lft = finish_times.max();

				DM("est: " << est << std::endl
				<< "lst: " << lst << std::endl
				<< "eft: " << eft << std::endl
				<< "lft: " << lft << std::endl);

				std::vector<Time> ca, pa;

				pa.push_back(eft);
				ca.push_back(lft);

				// skip first element in from.core_avail
				for (int i = 1; i < from.core_avail.size(); i++) {
					pa.push_back(std::max(est, from.core_avail[i].min()));
					ca.push_back(std::max(est, from.core_avail[i].max()));
				}

				// update scheduled jobs
				// keep it sorted to make it easier to merge
				bool added_j = false;
				for (const auto& rj : from.certain_jobs) {
					auto x = rj.first;
					auto x_eft = rj.second.min();
					auto x_lft = rj.second.max();
					if (contains(predecessors, x)) {
						if (lst < x_lft) {
							auto pos = std::find(ca.begin(), ca.end(), x_lft);
							if (pos != ca.end())
								*pos = lst;
						}
					} else if (lst <= x_eft) {
						if (!added_j && rj.first > j) {
							// right place to add j
							certain_jobs.emplace_back(j, finish_times);
							added_j = true;
						}
						certain_jobs.emplace_back(rj);
					}
				}
				// if we didn't add it yet, add it at the back
				if (!added_j)
					certain_jobs.emplace_back(j, finish_times);


				// sort in non-decreasing order
				std::sort(pa.begin(), pa.end());
				std::sort(ca.begin(), ca.end());

				for (int i = 0; i < from.core_avail.size(); i++) {
					DM(i << " -> " << pa[i] << ":" << ca[i] << std::endl);
					core_avail.emplace_back(pa[i], ca[i]);
				}

				add_pred(j, finish_times);

				assert(core_avail.size() > 0);
				DM("*** new state: constructed " << *this << std::endl);
			}

			hash_value_t get_key() const
			{
				return lookup_key;
			}

			// RV: this test is used to detect hash key collisions for states.
			//     depending on the data structure of scheduled_jobs, it might be expensive.
			//     to avoid this test, perhaps the key generation for jobs could be improved
			//     or the hash keys for jobs can be checked for collisions.
			bool same_jobs_scheduled(const Schedule_state &other) const
			{
				return scheduled_jobs == other.scheduled_jobs;
			}

			bool can_merge_with(const Schedule_state<Time>& other) const
			{
				assert(core_avail.size() == other.core_avail.size());
				// RV: collect statics:
				//     1: total number of calls
				//     2: #calls with matching hash keys
				//     3: #calls with matching scheduled jobs (<2 means hash collisions)
				//     4: #calls with intersecting core_avail
				//     5: #calls with overlapping job_finish_times
				static StatCollect cmw_stats = StatCollect("can_merge_with",100000);
				cmw_stats.tick(1);
				if (get_key() != other.get_key())
					return false;
				cmw_stats.tick(2);
				for (int i = 0; i < core_avail.size(); i++)
					if (!core_avail[i].intersects(other.core_avail[i]))
						return false;
				cmw_stats.tick(3);
				//check if JobFinishTimes overlap
				bool result = check_overlap(other.job_finish_times);
				if (!result) return false;
				cmw_stats.tick(4);
				if (!same_jobs_scheduled(other))
					return false;
				cmw_stats.tick(5);
				cmw_stats.print();
				return true;
			}

			bool try_to_merge(const Schedule_state<Time>& other)
			{
				if (!can_merge_with(other))
					return false;

				for (int i = 0; i < core_avail.size(); i++)
					core_avail[i] |= other.core_avail[i];

				// vector to collect joint certain jobs
				std::vector<std::pair<Job_index, Interval<Time>>> new_cj;

				// walk both sorted job lists to see if we find matches
				auto it = certain_jobs.begin();
				auto jt = other.certain_jobs.begin();
				while (it != certain_jobs.end() &&
				       jt != other.certain_jobs.end()) {
					if (it->first == jt->first) {
						// same job
						new_cj.emplace_back(it->first, it->second | jt->second);
						it++;
						jt++;
					} else if (it->first < jt->first)
						it++;
					else
						jt++;
				}
				// move new certain jobs into the state
				certain_jobs.swap(new_cj);

				// merge job_finish_times
				widen_overlap(other.job_finish_times);

				DM("+++ merged " << other << " into " << *this << std::endl);

				return true;
			}

			const unsigned int number_of_scheduled_jobs() const
			{
				return num_jobs_scheduled;
			}

			Interval<Time> core_availability() const
			{
				assert(core_avail.size() > 0);
				return core_avail[0];
			}

                        // Find the offset in the JobFinishTimes vector where the index should be located.
			int jft_find(const JobFinishTimes &jft, const Job_index pred_job) const
                        {
                                int start=0;
                                int end=jft.size();
                                while (start<end) {
                                        int mid=(start+end)/2;
                                        if (jft[mid].first == pred_job)
                                                return mid;
                                        else if (jft[mid].first < pred_job)
                                                start = mid+1;  // mid is too small, mid+1 might fit.
                                        else
                                                end = mid;
                                }
                                return start;
                        }


		  
			bool get_finish_times(Job_index j, Interval<Time> &ftimes) const
			{
			  int offset = jft_find(certain_jobs, j);
			  if (offset<certain_jobs.size() && certain_jobs[offset].first == j)
			  {
			    ftimes=certain_jobs[offset].second;
			    return true;
			  } else {
			    return false;
			  }
			}

			const bool job_incomplete(Job_index j) const
			{
				return !scheduled_jobs.contains(j);
			}

			const bool job_ready(const Job_precedence_set& predecessors) const
			{
				for (auto j : predecessors)
					if (!scheduled_jobs.contains(j))
						return false;
				return true;
			}

			friend std::ostream& operator<< (std::ostream& stream,
			                                 const Schedule_state<Time>& s)
			{
				stream << "Global::State(";
				for (const auto& a : s.core_avail)
					stream << "[" << a.from() << ", " << a.until() << "] ";
				stream << "(";
				for (const auto& rj : s.certain_jobs)
					stream << rj.first << "";
				stream << ") " << s.scheduled_jobs << ")";
				stream << " @ " << &s;
				return stream;
			}

			void print_vertex_label(std::ostream& out,
				const typename Job<Time>::Job_set& jobs) const
			{
				for (const auto& a : core_avail)
					out << "[" << a.from() << ", " << a.until() << "] ";
				out << "\\n";
				bool first = true;
				out << "{";
				for (const auto& rj : certain_jobs) {
					if (!first)
						out << ", ";
					out << "T" << jobs[rj.first].get_task_id()
					    << "J" << jobs[rj.first].get_job_id() << ":"
					    << rj.second.min() << "-" << rj.second.max();
					first = false;
				}
				out << "}";
			}

			// Functions for pathwise self-suspending exploration

		        // RV:  note that job_finish_times is an unordered_map<Job_index, Interval<Time>>.
		        //      similar to certain_jobs, it could use a vector<std::pair<Job_index,Interval<Time>>> .
		        //      Check the impact on memory and computation time.
		        //      Used operations:  find(), erase(), end().
		  // Reuse code of uni/state.hpp
			void add_pred(const Job_index pred_job, Interval<Time> ft)
			{
				// keep sorted according to Job_index.
				int it = jft_find(job_finish_times, pred_job);
				job_finish_times.insert(job_finish_times.begin()+it, std::make_pair(pred_job, ft));
			}

                        void add_pred_list(JobFinishTimes jft_list)
                        {
			  // RV: add_pred_list() is probably called when job_finish_times is empty.
                                for(auto e: jft_list)
                                {
                                        add_pred(e.first, e.second);
                                }
                        }

			void del_pred(const Job_index pred_job)
			{
				int it = jft_find(job_finish_times, pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job)
					job_finish_times.erase(job_finish_times.begin()+it);
			}

			// RV:  job_finish_times has Job_index, not JobID. Unclear how this works.
			void widen_pathwise_job(const JobID pred_job, const Interval<Time> ft)
			{
				int it = jft_find(job_finish_times, pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					(job_finish_times[it].second).widen(ft);
				}
			}

			// RV: Job_index instead of JobID?
			bool pathwisejob_exists(const JobID pred_job) const
			{
				int it = jft_find(job_finish_times, pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					return true;
				}
				return false;
			}

			// RV: similar to get_finish_times() in global/state.hpp .
			bool pathwisejob_exists(const Job_index pred_job, Interval<Time> &ft) const
			{
				int it = jft_find(job_finish_times, pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					ft = job_finish_times[it].second;
					return true;
				}
				return false;
                        }

			const Interval<Time>& get_pathwisejob_ft(const Job_index pathwise_job) const
			{
				int it = jft_find(job_finish_times, pathwise_job);
				//if (it < job_finish_times.size() && job_finish_times[it].first == pathwise_job)
                                return job_finish_times[it].second;
			}

		  // RV: Job_index instead of JobID?
		  //     This function would return the provided parameter, assuming it is found.
		  // not used
			const Job_index get_pathwisejob_job(const Job_index pathwise_job) const
			{
				return pathwise_job;
			}

			const JobFinishTimes& get_pathwise_jobs() const
			{
				return job_finish_times;
			}

			// Either check whether the job_finish_times overlap.
			bool check_overlap(const JobFinishTimes& from_pwj) const
			{
				bool allJobsIntersect = true;
				// The JobFinishTimes vectors are sorted.
				// Check intersect for matching jobs.
				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();
				while (from_it != from_pwj.end() &&
					state_it != job_finish_times.end())
				{
					if (from_it->first == state_it->first)
					{
						if (!from_it->second.intersects(state_it->second))
						{
							allJobsIntersect = false;
							break;
						}
						from_it++;
						state_it++;
					}
					else if (from_it->first < state_it->first)
						from_it++;
					else
						state_it++;
				}
				return allJobsIntersect;
			}

			void widen_overlap(const JobFinishTimes& from_pwj)
			{
				// The JobFinishTimes vectors are sorted.
				// Assume check_overlap() is true.
				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();
				while (from_it != from_pwj.end() &&
					state_it != job_finish_times.end())
				{
					if (from_it->first == state_it->first)
					{
						state_it->second.widen(from_it->second);
						from_it++;
						state_it++;
					}
					else if (from_it->first < state_it->first)
						from_it++;
					else
						state_it++;
				}
			}
		  

			// End of functions for pathwise self-suspending exploration

			private:

			const unsigned int num_jobs_scheduled;

			// set of jobs that have been dispatched (may still be running)
			const Index_set scheduled_jobs;

			// imprecise set of certainly running jobs
			std::vector<std::pair<Job_index, Interval<Time>>> certain_jobs;

			// system availability intervals
			std::vector<Interval<Time>> core_avail;

			const hash_value_t lookup_key;

			// no accidental copies
			Schedule_state(const Schedule_state& origin)  = delete;
		};

		template<class Time> class Schedule_node
		{
			private:

			Time earliest_pending_release;

			Job_set scheduled_jobs;
			hash_value_t lookup_key;
			Interval<Time> finish_time;

			// no accidental copies
			Schedule_node(const Schedule_node& origin)  = delete;

			typedef Schedule_state<Time> State;

			struct eft_compare 
			{
				bool operator() (State* x, State* y) const
				{
					return x->earliest_finish_time() < y->earliest_finish_time();
				}
			};

			typedef typename std::multiset<State*, eft_compare> State_ref_queue;
			State_ref_queue states;

			public:

			// initial node
			Schedule_node()
			: lookup_key{0}
			, finish_time{0,0}
			, earliest_pending_release{0}
			{
			}

			// transition: new node by scheduling a job in an existing state
			Schedule_node(
				const Schedule_node& from,
				const State& from_state,
				const Job<Time>& j,
				std::size_t idx,
				Interval<Time> ftimes,
				const Time next_earliest_release)
			: scheduled_jobs{from.scheduled_jobs, idx}
			, lookup_key{from.next_key(j)}
			, finish_time{ftimes}
			, earliest_pending_release{next_earliest_release}
			{
			}

			Time earliest_job_release() const
			{
				return earliest_pending_release;
			}

			hash_value_t get_key() const
			{
				return lookup_key;
			}

			const Job_set& get_scheduled_jobs() const
			{
				return scheduled_jobs;
			}

			bool matches(const Schedule_node& other) const
			{
				return lookup_key == other.lookup_key &&
					   scheduled_jobs == other.scheduled_jobs;
			}

			hash_value_t next_key(const Job<Time>& j) const
			{
				return get_key() ^ j.get_key();
			}

			const Interval<Time>& finish_range() const
			{
				return finish_time;
			}

			void add_state(State* s)
			{
				states.insert(s);
			}

			friend std::ostream& operator<< (std::ostream& stream,
											 const Schedule_node<Time>& n)
			{
				stream << "Node(" << n.states.size() << ")";
				return stream;
			}

			//return the number of states in the node
			int states_size() const
			{
				return states.size();
			}

			const State* get_first_state() const
			{
				auto first = states.begin();
				return *first;
			}

			const State* get_last_state() const
			{
				auto last = --(states.end());
				return *last;
			}

			const State_ref_queue* get_states() const
			{
				return &states;
			}

			bool merge_states(const Interval<Time> &new_st, const Schedule_state<Time> &from, const Job<Time>& sched_job)
			{
				// #NS# I am sorry for this merge function, but I am not sure how to make it better. This merge funciton is messy
				// because I am trying to merge if possible, before even creating a new state so that I don't have to create a new state
				// and then throw it away if it can be merged cause i do not know how to handle memory management in that case.
				// RV: instead of merging with only one state, try to merge with more states if possible.
				int merge_budget = states.size();
				int extra_budget = 0;  // Once merged, how many states should still be checked.

				static StatCollect stats = StatCollect("merge");
				stats.tick(merge_budget);

				bool result = false;
				for (auto& state : states)
				{
					if (merge_budget <= 0)
						break;
					if (new_st.intersects(state->finish_range()))
					{
						Interval<Time> ival{0,0};

						if (state->pathwisejob_exists(sched_job.get_job_index(),ival))
						{
							if (!new_st.intersects(ival))
								continue;
						}
						// First perform the check.
						if (state->check_or_widen(from.get_pathwise_jobs(), false)) {
							// When all matching jobs intersect, widen them.
							state->check_or_widen(from.get_pathwise_jobs(), true);
							// Find sched_job and widen its finish time interval
							state->widen_pathwise_job(sched_job.get_job_index(), new_st);

							//Widen the main finish range
							state->update_finish_range(new_st);

							result = true;

							// Try to merge with a few more states.
							// std::cerr << "Merged with " << merge_budget << " of " << states.size() << " states left.\n";
							merge_budget = extra_budget;
							extra_budget -= 2;
						}
					}
					merge_budget--;
				}

				stats.tick(result);
				stats.print();

				return result;
			}

		};

	}
}
