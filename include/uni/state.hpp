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
	  
		typedef Index_set Job_set;
		template<class Time> class Schedule_state;
		template<class Time> class Schedule_node;

		template<class Time> class Schedule_state
		{
			typedef std::vector<std::pair<const Job<Time>*, Interval<Time>>> Successors_list;
			typedef std::vector<Successors_list> Successors;

			private:

			// holds the finish time interval of the state
			Interval<Time> finish_time;

			// holds information about the state being a dead end or not
			// It is first intialized to be a deadend, if a job is found that can be scheduled from this state then
			// will converted to false.
			bool deadend = true;

			// keeps track of the earliest time a job with at least one predecessor becomes certainly ready
			Time earliest_certain_successor_jobs_ready_time;

			// job_finish_times holds the finish times of all the jobs that still have unscheduled successor
			typedef std::vector<std::pair<Job_index, Interval<Time>>> JobFinishTimes;
			JobFinishTimes job_finish_times;  

			// Find the offset in the job_finish_times vector where the index should be located.
			int jft_find(const Job_index pred_job) const
			{
				int start = 0;
				int end = job_finish_times.size();
				while (start < end) {
					int mid = (start + end) / 2;
					if (job_finish_times[mid].first == pred_job)
						return mid;
					else if (job_finish_times[mid].first < pred_job)
						start = mid + 1;  // mid is too small, mid+1 might fit.
					else
						end = mid;
				}
				return start;
			}

			// no accidental copies
			Schedule_state(const Schedule_state& origin) = delete;

			public:

			// initial state
			Schedule_state()
			: finish_time{0, 0}
			, earliest_certain_successor_jobs_ready_time{ Time_model::constants<Time>::infinity() }
			{
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_state(
				const Schedule_state<Time>& from, 
				const Job_index dispatched_j, 
				Interval<Time> ftime_interval, 
				const Job_set& scheduled_jobs,
				const Successors& successors_of)
			:finish_time{ftime_interval}
			{
				job_finish_times.reserve(from.job_finish_times.size() + 1);

				// updates the list of finish times of jobs with successors w.r.t. the previous system state
				// and calculates the earliest time a job with precedence constraints will become ready
				bool added_j = false;
				earliest_certain_successor_jobs_ready_time = Time_model::constants<Time>::infinity();
				for (auto ft : from.job_finish_times)
				{
					auto job = ft.first;
					auto lft = ft.second.max();

					if (!added_j && job > dispatched_j)
					{
						bool successor_pending = false;
						for (auto succ : successors_of[dispatched_j]) 
						{
							successor_pending = true;
							Time max_susp = succ.second.max();
							earliest_certain_successor_jobs_ready_time = 
								std::min(earliest_certain_successor_jobs_ready_time, 
									std::max(succ.first->latest_arrival(), ftime_interval.max() + max_susp));
						}
						if (successor_pending)
							job_finish_times.push_back(std::make_pair(dispatched_j, ftime_interval));
						added_j = true;
					}

					bool successor_pending = false;
					for (auto succ : successors_of[job]) {
						auto to_job = succ.first->get_job_index();
						if (!scheduled_jobs.contains(to_job))
						{
							successor_pending = true;
							Time max_susp = succ.second.max();
							earliest_certain_successor_jobs_ready_time = 
								std::min(earliest_certain_successor_jobs_ready_time, 
									std::max(succ.first->latest_arrival(), lft + max_susp));
						}
					}
					if (successor_pending)
						job_finish_times.push_back(std::make_pair(job, ft.second));
				}

				if (!added_j)
				{
					bool successor_pending = false;
					for (auto succ : successors_of[dispatched_j])
					{
						successor_pending = true;
						Time max_susp = succ.second.max();
						earliest_certain_successor_jobs_ready_time = 
							std::min(earliest_certain_successor_jobs_ready_time, 
								std::max(succ.first->latest_arrival(), ftime_interval.max() + max_susp));
					}
					if (successor_pending)
						job_finish_times.push_back(std::make_pair(dispatched_j, ftime_interval));
				}
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

			Time get_earliest_certain_successor_jobs_ready_time() const
			{
				return earliest_certain_successor_jobs_ready_time;
			}

			void update_successor_job_ready_time(Time t)
			{
				earliest_certain_successor_jobs_ready_time = std::max(earliest_certain_successor_jobs_ready_time, t);
			}

			// #NS# all the following functions are purely to handle the job_finish_times
			const JobFinishTimes& get_pathwise_jobs() const
			{
				return job_finish_times;
			}

			
			void widen_pathwise_job(const Job_index job, const Interval<Time> ft)
			{
				auto it = jft_find(job);
				if (it < job_finish_times.size() && job_finish_times[it].first == job)
					job_finish_times[it].second.widen(ft);
			}

			// Checks if the state kept information on the finishing time interval of job in the current system state,
			// and returns the finishing time interval in the variable 'ft' if that is the case
			bool pathwisejob_exists(const Job_index job, Interval<Time> &ft) const
			{
				auto it = jft_find(job);
				if (it < job_finish_times.size() && job_finish_times[it].first == job)
				{
					ft = job_finish_times[it].second;
					return true;
				}
				return false;
			}

			// RV: only called after a check that job exists.
			const Interval<Time>& get_pathwisejob_ft(const Job_index job) const
			{
				assert(jft_find(job) < job_finish_times.size());
				return job_finish_times[jft_find(job)].second;
			}

			// Check whether the job_finish_times intersect
			bool can_merge_with(const JobFinishTimes& from_pwj)
			{
				bool allJobsIntersect = true;
				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();

				while (from_it != from_pwj.end() &&	state_it != job_finish_times.end())
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

			// Check whether the job_finish_times can be merged and merge them if yes.
			bool try_and_merge(const JobFinishTimes& from_pwj)
			{
				if (!can_merge_with(from_pwj))
					return false;

				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();

				while (from_it != from_pwj.end() &&	state_it != job_finish_times.end())
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
				return true;	  
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
			Time next_certain_successor_job_ready_time;
			Time next_certain_source_job_release;

			Job_set scheduled_jobs;
			hash_value_t lookup_key;

			// finish_time in a node keeps track of the EFT and LFT of the last scheduled job 
			// in any state in the node
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
			Schedule_node(
				const Time next_earliest_release = 0,
				const Time next_certain_source_job_release = Time_model::constants<Time>::infinity() // the next time a job without predecessor is certainly released
			)
			: lookup_key{0}
			, finish_time{0,0}
			, earliest_pending_release{ next_earliest_release }
			, next_certain_successor_job_ready_time{ Time_model::constants<Time>::infinity() }
			, next_certain_source_job_release{ next_certain_source_job_release }
			{
			}

			// transition: new node by scheduling a job in an existing state
			Schedule_node(
				const Schedule_node& from,
				const Job<Time>& j,
				std::size_t idx,
				const Time next_earliest_release,
				const Time next_certain_source_job_release // the next time a job without predecessor is certainly released
			)
			: scheduled_jobs{from.scheduled_jobs, idx}
			, lookup_key{from.next_key(j)}
			, finish_time{0, Time_model::constants<Time>::infinity() }
			, earliest_pending_release{next_earliest_release}
			, next_certain_source_job_release{next_certain_source_job_release}
			, next_certain_successor_job_ready_time{ Time_model::constants<Time>::infinity() }
			{
			}

			~Schedule_node()
			{
				for (State* s : states)
					delete s;
			}

			Time earliest_job_release() const
			{
				return earliest_pending_release;
			}

			Time get_next_certain_source_job_release() const
			{
				return next_certain_source_job_release;
			}

			Time next_certain_job_ready_time() const
			{
				return std::min(next_certain_successor_job_ready_time, next_certain_source_job_release);
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

			Interval<Time> finish_range() const
			{
				return finish_time;
			}

			Time get_earliest_core_availability() const
			{
				return finish_time.min();
			}

			Time get_latest_core_availability() const
			{
				return finish_time.max();
			}

			void add_state(State* s)
			{
				if (states.empty())
				{
					finish_time = s->finish_range();
					next_certain_successor_job_ready_time = s->get_earliest_certain_successor_jobs_ready_time();
				}
				else
				{
					finish_time.widen(s->finish_range());
					next_certain_successor_job_ready_time = std::max(next_certain_successor_job_ready_time, s->get_earliest_certain_successor_jobs_ready_time());
				}	
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

			bool merge_states(const Schedule_state<Time>& s)
			{
				// RV: instead of merging with only one state, try to merge with more states if possible.
				int merge_budget = 1;
				static StatCollect stats = StatCollect("merge");
				stats.tick(merge_budget);

				Interval<Time> ft = s.finish_range();
				State* last_state_merged;

				bool result = false;
				for (auto& state : states)
				{
					if (ft.intersects(state->finish_range()))
					{
						if (result == false) // if we did not merge with any state yet
						{
							// When the finish time intervals of all jobs with unfinished successors intersect in both states, 
							// widen their finish time intervals in the existing state.
							if (state->try_and_merge(s.get_pathwise_jobs()))
							{
								//Widen the main finish range in the merged state and the node
								state->update_finish_range(ft);
								finish_time.widen(ft);
								//update the certain next job ready time in the merged state and node
								state->update_successor_job_ready_time(s.get_earliest_certain_successor_jobs_ready_time());
								next_certain_successor_job_ready_time = std::max(next_certain_successor_job_ready_time, s.get_earliest_certain_successor_jobs_ready_time());

								result = true;

								// Try to merge with a few more states.
								// std::cerr << "Merged with " << merge_budget << " of " << states.size() << " states left.\n";
								merge_budget--;
								if (merge_budget == 0)
									break;

								ft = state->finish_range();
								last_state_merged = state;
							}
						}
						else // if we already merged with one state at least
						{
							// When the finish time intervals of all jobs with unfinished successors intersect in both states, 
							// widen their finish time intervals in the existing state.
							if (state->try_and_merge(last_state_merged->get_pathwise_jobs()))
							{
								//Widen the main finish range in the merged state
								state->update_finish_range(last_state_merged->finish_range());
								//update the certain next job ready time in the merged state
								state->update_successor_job_ready_time(last_state_merged->get_earliest_certain_successor_jobs_ready_time());

								// the state was merged => we can thus remove the old one from the list of states
								states.erase(last_state_merged);
								delete last_state_merged;

								// Try to merge with a few more states.
								// std::cerr << "Merged with " << merge_budget << " of " << states.size() << " states left.\n";
								merge_budget--;
								if (merge_budget == 0)
									break;

								ft = state->finish_range();
								last_state_merged = state;
							}
						}
					}
				}
				stats.tick(result);
				stats.print();

				return result;
			}

		};

	}
}