#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>

#include <set>

#include "index_set.hpp"
#include "jobs.hpp"
#include "cache.hpp"
#include "statistics.hpp"

namespace NP {

	// use pointers as a primitive form of unique ID

	namespace Uniproc {

		// Similar to global/state.hpp, an Job_index is the index into the Job array.
		// typedef std::size_t Job_index;
	  
		typedef Index_set Job_set;
		template<class Time> class Schedule_state;
		template<class Time> class Schedule_node;

		template<class Time> class Schedule_state
		{
			private:

			// holds tthe finish time interval of the state
			Interval<Time> finish_time;
			// holds information about the state being a dead end or not
			// It is first intialized to be a deadend, if a job is found that can be scheduled from this state then
			// will converted to false.
			bool deadend = true;

			// no accidental copies
			Schedule_state(const Schedule_state& origin)  = delete;


			// #NS# job_finish_times holds the finish times of all the jobs that still have unscheduled successor
			// It does not need to remember this all the time, there can be a flag that says i don't need anything to be remembered
			// and in those cases these are not assigned any value.
		  // RV: this data structure is similar to the certain_jobs from global/set.
		  //     instead of unordered_map, using vector<
		  // typedef typename std::unordered_map<NP::Job_index, Interval<Time>> JobFinishTimes;
			typedef typename std::vector<std::pair<Job_index,Interval<Time>>> JobFinishTimes;
			JobFinishTimes job_finish_times;

			public:

			// initial state
			Schedule_state()
			: finish_time{0, 0}
			{
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_state(Interval<Time> ftimes)
			:finish_time{ftimes}
			{
			}

			bool is_deadend()
			{
				return deadend;
			}

			void not_deadend()
			{
				deadend = false;
			}

			Time earliest_finish_time() const
			{
				return finish_time.from();
			}

			Time latest_finish_time() const
			{
				return finish_time.until();
			}

			const Interval<Time>& finish_range() const
			{
				return finish_time;
			}

			void update_finish_range(const Interval<Time> &update)
			{
				assert(update.intersects(finish_time));
				finish_time.widen(update);
			}

		  // Unused
			void copy_state(const Interval<Time> &newst)
			{
				finish_time.equate(newst);
			}

			// #NS# all the following functions are purely to handle the job_finish_times

			// Find the offset in the job_finish_times vector where the index should be located.
			int jft_find(const Job_index pred_job) const
			{
				int start=0;
				int end=job_finish_times.size();
				while (start<end) {
					int mid=(start+end)/2;
					if (job_finish_times[mid].first == pred_job)
						return mid;
					else if (job_finish_times[mid].first < pred_job)
						start = mid+1;  // mid is too small, mid+1 might fit.
					else
						end = mid;
	       			}
				return start;
			}

		        void add_pred(const Job_index pred_job, Interval<Time> ft)
			{
				// keep sorted according to Job_index.
				int it = jft_find(pred_job);
				job_finish_times.insert(job_finish_times.begin()+it, std::make_pair(pred_job, ft));
				// job_finish_times.emplace(pred_job, ft);
			}

			void add_pred_list(JobFinishTimes jft_list)
			{
				for(auto e: jft_list)
				{
					add_pred(e.first, e.second);
				}
			}

			void del_pred(const Job_index pred_job)
			{
				auto it = jft_find(pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job)
					job_finish_times.erase(job_finish_times.begin()+it);
			}

			void widen_pathwise_job(const Job_index pred_job, const Interval<Time> ft)
			{
				int it = jft_find(pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					(job_finish_times[it].second).widen(ft);
				}
			}

			bool pathwisejob_exists(const Job_index pred_job) const
			{
				int it = jft_find(pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					return true;
				}
				return false;
			}

			// RV: similar to get_finish_times() in global/state.hpp .
			bool pathwisejob_exists(const Job_index pred_job, Interval<Time> &ft) const
			{
				int it = jft_find(pred_job);
				if (it < job_finish_times.size() && job_finish_times[it].first==pred_job) {
					ft = job_finish_times[it].second;
					return true;
				}
				return false;
			}

			// RV: only called after a check that pathwise_job exists.
			const Interval<Time>& get_pathwisejob_ft(const Job_index pathwise_job) const
			{
				int it = jft_find(pathwise_job);
				//if (it < job_finish_times.size() && job_finish_times[it].first == pathwise_job)
				return job_finish_times[it].second;
				//	else
				// return NULL;
			}

			
		  
		  // Unused.
			const Job_index get_pathwisejob_job(const Job_index pathwise_job) const
			{
				return 0; // job_finish_times.find(pathwise_job)->first;
			}

			const JobFinishTimes& get_pathwise_jobs() const
			{
				return job_finish_times;
			}

			// Either check whether the job_finish_times can be merged or merge them without checking.
			bool check_or_widen(const JobFinishTimes& from_pwj, bool widen)
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
						// Same jobs.
						if (widen)
						{
							state_it->second.widen(from_it->second);
						}
						else
						{
							if (!from_it->second.intersects(state_it->second))
							{
								allJobsIntersect = false;
								break;
							}
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
		  		  
			friend std::ostream& operator<< (std::ostream& stream,
											 const Schedule_state<Time>& s)
			{
				stream << "State(" << s.finish_range() << ")";
				return stream;
			}
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
				//#define MERGE_STATISTICS
				static StatCollect stat = StatCollect("merge");
				stat.tick(merge_budget);
#ifdef MERGE_STATISTICS
				static int stats_count[50] = { 0 };
				static int stats_total = 0;
				static int stats_max = 0;
				static int stats_merge_ok = 0;
				static int stats_merge_fail = 0;
				bool stats_print=false;
				stats_total++;
				stats_print = (stats_total%500 == 0);
				if (merge_budget < 50)
				  stats_count[merge_budget]++;
				if (merge_budget > stats_max)
				{
				  stats_max = merge_budget;
				  stats_print = true;
				}
				if (stats_print) {
				  std::cerr << "merge stats (" << stats_total << ", "<< stats_max << ") : " << stats_count[0];
				  for (int i=1; i<50 && i<=stats_max; i++)
				    {
				      std::cerr << ", " << stats_count[i];
				    }
				  std::cerr << "\n"; 
				}
#endif
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
				stat.tick(result);
				stat.print();
#ifdef MERGE_STATISTICS
				if (result)
				  stats_merge_ok++;
				else
				  stats_merge_fail++;
				if (stats_print) {
				  std::cerr << "merge stats ("<< stats_total << ")  ok: "<< stats_merge_ok << ". fail: "<<stats_merge_fail << "\n";
				}
#endif
				return result;
			}

		};

	}
}
