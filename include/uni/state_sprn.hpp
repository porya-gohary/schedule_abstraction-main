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

		template<class Time> class Schedule_state
		{
			private:

			Interval<Time> finish_time;
			Schedule_node *ptrnode;

			// no accidental copies
			Schedule_state(const Schedule_state& origin)  = delete;

			public:

			// initial state
			Schedule_state(
				Schedule_node *nptr)
			: finish_time{0, 0}
			, ptrnode{nptr}
			{
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_state(
				std::size_t idx,
				Interval<Time> ftimes,
				Schedule_node *nptr)
			: finish_time{ftimes}
			, ptrnode{nptr}
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

			Time earliest_job_release() const
			{
				return ptrnode->earliest_job_release();
			}

			hash_value_t get_key() const
			{
				return ptrnode->get_key();
			}

			hash_value_t next_key(const Job<Time>& j) const
			{
				return ptrnode->get_key() ^ j.get_key();
			}

			const Job_set& get_scheduled_jobs() const
			{
				return ptrnode->get_scheduled_jobs();
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

			typedef Schedule_state<Time> State;
			typedef std::deque<State> States;

			// no accidental copies
			Schedule_node(const Schedule_node& origin)  = delete;

			public:

			// initial state
			Schedule_node()
			: lookup_key{0}
			, earliest_pending_release{0}
			{
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_node(
				const Schedule_node& from,
				const Job<Time>& j,
				std::size_t idx,
				std::deque<State> sts,
				const Time next_earliest_release)
			: scheduled_jobs{from.scheduled_jobs, idx}
			, lookup_key{from.next_key(j)}
			, earliest_pending_release{next_earliest_release}
			, States{sts}
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

			friend std::ostream& operator<< (std::ostream& stream,
			                                 const Schedule_node<Time>& s)
			{
				stream << "Node(" << s.get_scheduled_jobs() << ")";
				return stream;
			}
		};

	}
}
