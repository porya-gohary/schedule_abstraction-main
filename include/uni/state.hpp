#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>

#include <set>

#include "index_set.hpp"
#include "jobs.hpp"
#include "cache.hpp"

namespace NP {

	// use pointers as a primitive form of unique ID

	namespace Uniproc {

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

			typedef typename std::unordered_map<Job<Time>, Interval<Time>> JobFinishTimes;
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

			void copy_state(const Interval<Time> &newst)
			{
				finish_time.equate(newst);
			}

			void add_pred(const Job<Time> pred_job, const Interval<Time> ft)
			{
				job_finish_times[pred_job] = ft;
			}

			void add_pred_list(JobFinishTimes jft_list)
			{
				for(auto e: jft_list)
				{
					add_pred(e.first, e.second);
				}
			}

			void del_pred(const Job<Time> pred_job)
			{
				job_finish_times.erase(pred_job);
			}

			void widen_pathwise_job(const Job<Time> pred_job, const Interval<Time> ft)
			{
				auto widen_job = job_finish_times.find(pred_job);
				(widen_job->second).widen(ft);
			}

			const Interval<Time>& get_pathwisejob_ft(const Job<Time>& pathwise_job) const
			{
				return (job_finish_times.find(pathwise_job))->second;
			}

			const Job<Time>& get_pathwisejob_job(const Job<Time>& pathwise_job) const
			{
				return (job_finish_times.find(pathwise_job))->first;
			}

			const JobFinishTimes& get_pathwise_jobs() const
			{
				return job_finish_times;
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

			typedef typename std::set<State*, eft_compare> State_ref_queue;
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
				// typedef typename std::set<State *, eft_compare>::iterator State_ref;
				// State_ref s_ref = states.begin();

				// while (s_ref != states.end())
				// {
				// 	if (new_st.intersects((*s_ref)->finish_range()))
				// 	{
				// 		bool state_mergable = true;

				// 		for (const auto &job_finish_time_pair : (*s_ref)->get_pathwise_jobs())
				// 		{
				// 			const Job<Time> &pred_job = job_finish_time_pair.first;
				// 			const Interval<Time> &job_finish_time = job_finish_time_pair.second;

				// 			if (((from.get_pathwisejob_job(pred_job)).get_id() == pred_job.get_id() &&
				// 				 job_finish_time.intersects(from.get_pathwisejob_ft(pred_job))) ||
				// 				(pred_job.get_id() == sched_job.get_id() && job_finish_time.intersects(new_st)))
				// 			{
				// 				continue;
				// 			}
				// 			else
				// 			{
				// 				state_mergable = false;
				// 				break;
				// 			}
				// 		}

				// 		if (state_mergable)
				// 		{
				// 			State *merged_state = *s_ref;
				// 			merged_state->update_finish_range(new_st);
							
				// 			for (const auto &job_finish_time_pair : (merged_state)->get_pathwise_jobs())
				// 			{
				// 				const Job<Time> &pred_job = job_finish_time_pair.first;
				// 				const Interval<Time> &job_finish_time = job_finish_time_pair.second;

				// 				if ((from.get_pathwisejob_job(pred_job)).get_id() == pred_job.get_id())
				// 				{
				// 					(merged_state)->widen_pathwise_job(pred_job, job_finish_time);
				// 				}
				// 				else if (pred_job.get_id() == sched_job.get_id())
				// 				{
				// 					(merged_state)->widen_pathwise_job(pred_job,new_st);
				// 				}
				// 			}

				// 			s_ref = states.erase(s_ref);
				// 			states.insert(merged_state);
				// 			return true;
				// 		}
				// 		else
				// 		{
				// 			s_ref++;
				// 		}
				// 	}
				// 	else
				// 	{
				// 		s_ref++;
				// 	}
				// }

				return false;
			}

		};

	}
}