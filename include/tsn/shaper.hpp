#ifndef SHAPER_HPP
#define SHAPER_HPP

#include <ostream>
#include <vector>
#include <deque>
#include <algorithm> // for find
#include <functional> // for hash
#include <exception>
#include <cmath>

#include "time.hpp"
#include "config.h"

namespace NP {

	template<class Time> class Time_Aware_Shaper {

	public:
		typedef std::vector<Time_Aware_Shaper<Time>> TAS_set;
		typedef std::deque<Interval<Time>> Intervals;
		typedef Time Priority; // Make it a time value to support EDF
		typedef Time Period;

	private:
		Priority priority;
		Period period;
		Intervals tas_close_queue;

	public:

		Time_Aware_Shaper(Priority prio,
			Period per,
			Intervals tascq)
		: priority(prio), period(per), tas_close_queue(tascq)
		{
		}

		Time get_period() const
		{
			return period;
		}

		Time get_priority() const
		{
			return priority;
		}

		Time next_open(Time check) const
		{
			DM("Inside next_open: check "<<check<<"\n");
			Intervals gate_close = get_gates_close(check,check);

			for(auto gci: gate_close)
				DM("Gate:"<<gci<<" \n");

			for(Interval<Time> gc: gate_close)
			{
				if(check < gc.from())
					return check;
				else if(check >= gc.from() && check < gc.upto())
					return gc.upto();
			}
			return check;
		}

		Intervals get_gates_open(Time start, Time end) const
		{
			Intervals GO;
			Intervals GC = get_gates_close(start, end);
			Time start_period = floor(start/period);
			Time end_period = ceil(end/period);
			Time next_gate_open;

			if(start_period*period < GC[0].from())
				next_gate_open = start_period*period;
			else
				next_gate_open = GC[0].upto();

			for(Interval<Time> gc : GC)
			{
				GO.emplace_back(Interval<Time>{next_gate_open, gc.from()});
				next_gate_open = gc.upto();
			}

			if(next_gate_open < end_period*period)
				GO.emplace_back(Interval<Time>{next_gate_open, end_period*period});
			
			return GO;
		}

		Intervals get_gates_close(Time start, Time end) const
		{
			Intervals GC;
			Time start_period = floor(start/period);
			Time end_period = ceil(end/period);
			for(Time current_period = start_period; current_period <= end_period; current_period += 1 )
			{
				for(Interval<Time> gc: tas_close_queue)
				{
					GC.emplace_back(Interval<Time>{current_period*period + gc.from(), 
													current_period*period + gc.upto()});
				}
			}
			return GC;
		}
	};

	class InvalidTASParameter : public std::exception
	{
		public:

		InvalidTASParameter(const JobID& bad_id)
		: ref(bad_id)
		{}

		const JobID ref;

		virtual const char* what() const noexcept override
		{
			return "invalid TAS parameter";
		}

	};

	template<class Time>
	void validate_tas_refs(const typename Time_Aware_Shaper<Time>::TAS_set tasQueues,
	                         const typename Job<Time>::Job_set jobs)
	{
		if (tasQueues.size()>0)
		{
			for (auto job : jobs) {
				bool priority_check = false;
				for (auto& tas : tasQueues) {
					if (job.get_priority() == tas.get_priority()) {
						priority_check = true;
						break;
					}
				}
				if (priority_check == false)
					throw InvalidTASParameter(job.get_id());
			}
		}
	}

}

#endif
