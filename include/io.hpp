#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include <utility>

#include "interval.hpp"
#include "time.hpp"
#include "jobs.hpp"
#include "precedence.hpp"
#include "aborts.hpp"
#include "config.h"
#include "tsn/shaper.hpp"

namespace NP {

	inline void skip_over(std::istream& in, char c)
	{
		while (in.good() && in.get() != (int) c)
			/* skip */;
	}

	inline bool skip_one(std::istream& in, char c)
	{
		if (in.good() && in.peek() == (int) c) {
			in.get(); /* skip */
			return true;
		} else
			return false;
	}

	inline void skip_all(std::istream& in, char c)
	{
		while (skip_one(in, c))
			/* skip */;
	}

	inline bool more_data(std::istream& in)
	{
		in.peek();
		return !in.eof();
	}

	inline bool more_fields_in_line(std::istream& in)
	{
		if (in.peek() == (int)'\n' || in.peek() == (int)'\r')
			return false;
		else
			return true;
	}

	inline void next_field(std::istream& in)
	{
	  // RV: The added while will consume any sequence of spaces and separators.
	  //     Empty fields are not possible.
		while(in.peek() == ',' || in.peek()==' ')
		{
			// eat up any trailing spaces
			skip_all(in, ' ');
			// eat up field separator
			skip_one(in, ',');
		}
	}

	inline void next_line(std::istream& in)
	{
		skip_over(in, '\n');
	}

	inline JobID parse_job_id(std::istream& in)
	{
		unsigned long jid, tid;

		in >> tid;
		next_field(in);
		in >> jid;
		return JobID(jid, tid);
	}

	//Functions that help parse selfsuspending tasks file
	template<class Time>
	Precedence_constraint<Time> parse_precedence_constraint(std::istream &in)
	{
		unsigned long from_tid, from_jid, to_tid, to_jid;
		Time sus_min=0, sus_max=0;

		std::ios_base::iostate state_before = in.exceptions();

		in.exceptions(std::istream::failbit | std::istream::badbit);

		in >> from_tid;
		next_field(in);
		in >> from_jid;
		next_field(in);
		in >> to_tid;
		next_field(in);
		in >> to_jid;
		next_field(in);
		if (more_fields_in_line(in))
		{
			in >> sus_min;
			next_field(in);
			in >> sus_max;
		}

		in.exceptions(state_before);

		return Precedence_constraint<Time>{JobID{from_jid, from_tid},
											JobID{to_jid, to_tid},
		                          			Interval<Time>{sus_min, sus_max}};
	}

	template<class Time>
	std::vector<Precedence_constraint<Time>> parse_precedence_file(std::istream& in)
	{
		// skip column headers
		next_line(in);
		std::vector<Precedence_constraint<Time>> cstr;

		// parse all rows
		while (more_data(in)) {
			// each row contains one self-suspending constraint
			cstr.push_back(parse_precedence_constraint<Time>(in));
			next_line(in);
		}
		return cstr;
	}

	//Functions that help parse the main job workload
	template<class Time> Job<Time> parse_job(std::istream& in, std::size_t idx)
	{
		unsigned long tid, jid;

		std::ios_base::iostate state_before = in.exceptions();

		Time arr_min, arr_max, cost_min, cost_max, dl, prio;

		in.exceptions(std::istream::failbit | std::istream::badbit);

		in >> tid;
		next_field(in);
		in >> jid;
		next_field(in);
		in >> arr_min;
		next_field(in);
		in >> arr_max;
		next_field(in);
		in >> cost_min;
		next_field(in);
		in >> cost_max;
		next_field(in);
		in >> dl;
		next_field(in);
		in >> prio;

		in.exceptions(state_before);

		return Job<Time>{jid, Interval<Time>{arr_min, arr_max},
						Interval<Time>{cost_min, cost_max}, dl, prio, idx, tid};
	}

	template<class Time>
	typename Job<Time>::Job_set parse_file(std::istream& in)
	{
		// first row contains a comment, just skip it
		next_line(in);

		typename Job<Time>::Job_set jobs;
		std::size_t idx=0;

		while (more_data(in)) {
			jobs.push_back(parse_job<Time>(in, idx));
			idx++;
			// munge any trailing whitespace or extra columns
			next_line(in);
		}

		return jobs;
	}

	template<class Time> Time_Aware_Shaper<Time> parse_tas(std::istream& in)
	{
		Time prio, period, gb;
		Time gate_close, gate_open;
		Time curr_time = 0;
		bool tas, cbs, isvar;

		typename Time_Aware_Shaper<Time>::Intervals tas_queue;

		std::ios_base::iostate state_before = in.exceptions();

		in.exceptions(std::istream::failbit | std::istream::badbit);

		in >> prio;
		next_field(in);
		in >> period;
		next_field(in);
		in >> cbs;
		next_field(in);
		in >> tas;
		next_field(in);
		in >> isvar;
		next_field(in);
		in >> gb;
		next_field(in);

		if(tas == true)
		{
			int interval_num = 0;
			while(more_fields_in_line(in)) {
				in >> gate_close;
				next_field(in);
				in >> gate_open;
				if( in.peek() != '\n' && in.peek() != '\r')
					next_field(in);

				gate_close = curr_time + gate_close;
				curr_time = gate_close;

				gate_open = curr_time + gate_open;
				curr_time = gate_open;

				tas_queue.push_back(Interval<Time>{gate_close,gate_open});
			}
		}

		in.exceptions(state_before);
		return Time_Aware_Shaper<Time>{prio, period, tas, cbs, isvar, gb, tas_queue};
	}

	template<class Time>
	typename Time_Aware_Shaper<Time>::TAS_set parse_tas_file(std::istream& in)
	{
		// first row contains the header of each column, so skip it
		next_line(in);

		typename Time_Aware_Shaper<Time>::TAS_set TAS_queues;

		int i = 0;
		while (more_data(in)) {
			TAS_queues.push_back(parse_tas<Time>(in));
			// munge any trailing whitespace or extra columns
			next_line(in);
		}

		return TAS_queues;
	}

	//Functions that help parse the abort actions file
	template<class Time>
	Abort_action<Time> parse_abort_action(std::istream& in)
	{
		unsigned long tid, jid;
		Time trig_min, trig_max, cleanup_min, cleanup_max;

		std::ios_base::iostate state_before = in.exceptions();

		in.exceptions(std::istream::failbit | std::istream::badbit);

		in >> tid;
		next_field(in);
		in >> jid;
		next_field(in);
		in >> trig_min;
		next_field(in);
		in >> trig_max;
		next_field(in);
		in >> cleanup_min;
		next_field(in);
		in >> cleanup_max;

		in.exceptions(state_before);

		return Abort_action<Time>{JobID{jid, tid},
		                          Interval<Time>{trig_min, trig_max},
		                          Interval<Time>{cleanup_min, cleanup_max}};
	}


	template<class Time>
	std::vector<Abort_action<Time>> parse_abort_file(std::istream& in)
	{
		// first row contains a comment, just skip it
		next_line(in);

		std::vector<Abort_action<Time>> abort_actions;

		while (more_data(in)) {
			abort_actions.push_back(parse_abort_action<Time>(in));
			// munge any trailing whitespace or extra columns
			next_line(in);
		}

		return abort_actions;
	}


}

#endif
