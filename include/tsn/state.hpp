#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>

#include <set>

#include "index_set.hpp"
#include "jobs.hpp"
#include "cache.hpp"
#include "space.hpp"

namespace NP {

	// use pointers as a primitive form of unique ID

	namespace TSN {

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

			public:

			// initial state
			Schedule_state()
			: finish_time{0, 0}
			{
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_state(Interval<Time> ftimes)
			: finish_time{ftimes}
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
			std::vector<Time> trmax;
			const Job<Time>& job_0;


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

			// find trmax, find the first incomplete job in the list of jobs sorted by latest arrival
			void calculate_trmax(const std::deque<std::multimap<Time, const Job<Time>*>>& jobs_by_latest_arrival_priority) {
				trmax.reserve(jobs_by_latest_arrival_priority.size());
				for (int i = 0; i < jobs_by_latest_arrival_priority.size(); i++) {
					auto it = jobs_by_latest_arrival_priority[i].begin();
					for (; it != jobs_by_latest_arrival_priority[i].end(); ++it) {
						const Job<Time>& j = *(it->second);

						// not relevant if already scheduled
						if (scheduled_jobs.contains(index_of(j)))
							continue;

						trmax.push_back(j.latest_arrival());
						break;
					}
					if(it == jobs_by_latest_arrival_priority[i].end())
						trmax.push_back(Time_model::constants<Time>::infinity());
				}
			}

			std::size_t index_of(const Job<Time>& j)
			{
				return (std::size_t)(&j - &job_0);
			}

			public:

			// initial node
			Schedule_node(const Job<Time>& j0, const std::deque<std::multimap<Time, const Job<Time>*>>& jobs_by_latest_arrival_priority)
			: job_0{ j0 }
			, lookup_key{0}
			, finish_time{0,0}
			, earliest_pending_release{0}
			{
				calculate_trmax(jobs_by_latest_arrival_priority);
			}

			// transition: new node by scheduling a job in an existing state
			Schedule_node(
				const Schedule_node& from,
				const Job<Time>& j,
				Interval<Time> ftimes,
				const Time next_earliest_release,
				const std::deque<std::multimap<Time, const Job<Time>*>>& jobs_by_latest_arrival_priority)
			: job_0{from.job_0}
			, scheduled_jobs{from.scheduled_jobs, (std::size_t)(& j - &from.job_0)}
			, lookup_key{from.next_key(j)}
			, finish_time{ftimes}
			, earliest_pending_release{next_earliest_release}
			{
				calculate_trmax(jobs_by_latest_arrival_priority);
			}

			Time get_trmax(int priority) const
			{
				assert(priority < trmax.size());
				return trmax[priority];
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

			bool merge_states(const Interval<Time> &new_st)
			{
				typedef typename std::set<State*, eft_compare>::iterator State_ref;
				State_ref s_ref = states.begin();
				while(s_ref!=states.end())
				{
					if(new_st.intersects((*s_ref)->finish_range()))
					{
						State* merged_state = *s_ref;
						merged_state->update_finish_range(new_st);
						s_ref = states.erase(s_ref);
						states.insert(merged_state);
						DM("\nState merged\n");
						return true;
					}
					else
						s_ref++;
				}
				return false;
			}

		};

	}
}
