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

			Interval<Time> finish_time;

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

			public:

			typedef Schedule_state<Time> State;
			State st;

			typedef std::deque<State> States;
			States states;

			// initial node
			Schedule_node()
			: lookup_key{0}
			, earliest_pending_release{0}
			{
			}

			// transition: new node by scheduling a job in an existing state
			Schedule_node(
				const Schedule_node& from,
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

			void add_state(const State& s)
			{
				states.push_back(s);
			}

			friend std::ostream& operator<< (std::ostream& stream,
			                                 const Schedule_node<Time>& s)
			{
				stream << "Node(" << s.get_scheduled_jobs() << ")";
				return stream;
			}

			//return the number of states in the node
			int states_size() const
			{
				return states.size();
			}

			//return the 'i'th state
			State get_state(int i) const
			{
				st = states[i];
				return st;
			}

		};

	}
}
