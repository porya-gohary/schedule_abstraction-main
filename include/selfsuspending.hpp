#ifndef SELFSUSPENDING_HPP
#define SELFSUSPENDING_HPP

#include "jobs.hpp"

namespace NP {

	template<class Time>
	class Suspending_Task{

	public:

		Suspending_Task(JobID from,
						JobID to,
						Interval<Time> sus_times)
		: from(from)
		, to(to)
		, sus_times(sus_times)
		{
		}

		JobID get_fromID() const
		{
			return from;
		}

		JobID get_toID() const
		{
			return to;
		}

		Time get_minsus() const
		{
			return sus_times.from();
		}

		Time get_maxsus() const
		{
			return sus_times.until();
		}

	private:
		JobID from;
		JobID to;
		Interval<Time> sus_times;

	};

	template<class Time>
	void validate_susp_refs(std::vector<Suspending_Task<Time>>& susps,
	                        const typename Job<Time>::Job_set jobs)
	{
		for (auto st : susps) {
			lookup<Time>(jobs, st.get_fromID());
			lookup<Time>(jobs, st.get_toID());
		}
	}

}

#endif
