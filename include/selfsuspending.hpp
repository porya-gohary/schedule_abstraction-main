#ifndef SELFSUSPENDING_HPP
#define SELFSUSPENDING_HPP

#include <vector>

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

		Interval<Time> get_suspension() const
		{
			return sus_times;
		}

		void set_toIndex(Job_index index)
	  	{
			toIndex = index;
		}

		Job_index get_toIndex() const
		{
			return toIndex;
		}

		void set_fromIndex(Job_index index)
	  	{
			fromIndex = index;
		}

		Job_index get_fromIndex() const
		{
			return fromIndex;
		}

	private:
		JobID from;
		JobID to;
		// RV: toIndex and fromIndex are set during validation
		Job_index toIndex;
		Job_index fromIndex;

		Interval<Time> sus_times;

	};

	class Invalid_Self_Suspending_Parameter : public std::exception
	{
		public:

		Invalid_Self_Suspending_Parameter(const JobID& bad_id)
		: ref(bad_id)
		{}

		const JobID ref;

		virtual const char* what() const noexcept override
		{
			return "invalid self-suspending reference";
		}

	};

	template<class Time>
	void validate_susp_refs(std::vector<Suspending_Task<Time>>& susps,
							const typename Job<Time>::Job_set jobs)
	{
	
	  for (int idx=0; idx<susps.size(); idx++) {
			const Job<Time>& fromJob = lookup<Time>(jobs, susps[idx].get_fromID());
			const Job<Time>& toJob = lookup<Time>(jobs, susps[idx].get_toID());
			// RV: after validation, susps should not be changed.
			//     set toIndex and fromIndex here.
			susps[idx].set_toIndex((Job_index)(&toJob - &(jobs[0])));
			susps[idx].set_fromIndex((Job_index)(&fromJob - &(jobs[0])));
			if (susps[idx].get_maxsus() < susps[idx].get_minsus()) {
				throw Invalid_Self_Suspending_Parameter(susps[idx].get_fromID());
			}
		}
	}

}

#endif
