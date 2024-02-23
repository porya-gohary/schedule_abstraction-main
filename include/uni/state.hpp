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

			typedef typename std::unordered_map<JobID, Interval<Time>> JobFinishTimes;
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

			void add_pred(const JobID pred_job, Interval<Time> ft)
			{
				job_finish_times.emplace(pred_job, ft);
			}

			void add_pred_list(JobFinishTimes jft_list)
			{
				for(auto e: jft_list)
				{
					add_pred(e.first, e.second);
				}
			}

			void del_pred(const JobID pred_job)
			{
				job_finish_times.erase(pred_job);
			}

			void widen_pathwise_job(const JobID pred_job, const Interval<Time> ft)
			{
				auto it = job_finish_times.find(pred_job);
				if (it != job_finish_times.end()) {
					(it->second).widen(ft);
				}
			}

			bool pathwisejob_exists(const JobID pred_job) const
			{
				auto it = job_finish_times.find(pred_job);
				if (it != job_finish_times.end()) {
					return true;
				}
				return false;
			}

			const Interval<Time>& get_pathwisejob_ft(const JobID& pathwise_job) const
			{
				return job_finish_times.find(pathwise_job)->second;
			}

			const JobID& get_pathwisejob_job(const JobID& pathwise_job) const
			{
				return job_finish_times.find(pathwise_job)->first;
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
				for (auto& state : states)
				{
					if (new_st.intersects(state->finish_range()))
					{
						bool allJobsIntersect = true;
						for (const auto& job : from.get_pathwise_jobs())
						{
							if (state->pathwisejob_exists(job.first))
							{
								if (!job.second.intersects(state->get_pathwisejob_ft(job.first)))
								{
									allJobsIntersect = false;
									break;
								}
							}
						}
						if (state->pathwisejob_exists(sched_job.get_id()))
						{
							if (!new_st.intersects(state->get_pathwisejob_ft(sched_job.get_id())))
							{
								allJobsIntersect = false;
							}
						}

						if (allJobsIntersect)
						{
							//Widen finish time intervals for all pathwise jobs
							for (auto& job : from.get_pathwise_jobs())
							{
								state->widen_pathwise_job(job.first, job.second);
							}

							//Find sched_job and widen its finish time interval
							state->widen_pathwise_job(sched_job.get_id(), new_st);

							//Widen the main finish range
							state->update_finish_range(new_st);

							return true;
						}
					}
				}
				return false;
			}

		};

	}
}