#include "doctest.h"

#include <iostream>
#include <sstream>

#include "problem.hpp"
#include "io.hpp"
#include "interval.hpp"
#include "tsn/space.hpp"

using namespace NP;

static const auto inf = Time_model::constants<dtime_t>::infinity();

const std::string TAS_constGB_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"   1,  1,    0,    0,    1,    2,  10,    0\n"
"   1,  2,   10,   10,    1,    2,  20,    0\n"
"   1,  3,   20,   20,    1,    2,  30,    0\n"
"   2,  1,   10,   12,    5,    7,  30,    1\n"
"   3,  1,   25,   25,   15,   19,  60,    2\n"
"   4,  1,   10,   12,    5,    7,  10,    0\n";

const std::string TAS_constGB_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     0,   1,    2,  2, 0, 0\n"
"   1,  30,   0,   1,     0,   1,   11,  4\n"
"   2,  60,   0,   1,     0,   1,   20, 17\n";

TEST_CASE("[TSN] Constant Guardband for TAS") {
	auto jobs_in = std::istringstream(TAS_constGB_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_constGB_shaper_file);

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(jobs_in),
		parse_dag_file(dag_in),
		parse_abort_file<dtime_t>(aborts_in),
		parse_tas_file<dtime_t>(shaper_in),
		1
	};

	Analysis_options opts;
	opts.early_exit = false;

	auto space = TSN::State_space<dtime_t>::explore(prob, opts);
	CHECK(!space.is_schedulable());

	CHECK(space.get_finish_times(prob.jobs[2]).min() == 21);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 31);

};

const std::string TAS_varGB_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"   1,  1,    0,    0,    1,    2,  10,    0\n"
"   1,  2,   10,   10,    1,    2,  20,    0\n"
"   1,  3,   20,   20,    1,    2,  30,    0\n"
"   2,  1,   10,   12,    5,    7,  30,    1\n"
"   3,  1,   25,   25,   15,   19,  60,    2\n";

const std::string TAS_varGB_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     1,   0,    2,  2, 0, 0\n"
"   1,  30,   0,   1,     1,   0,   11,  4\n"
"   2,  60,   0,   1,     1,   0,   20, 17\n";

TEST_CASE("[TSN] Variable Guardband for TAS") {
	auto jobs_in = std::istringstream(TAS_varGB_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_varGB_shaper_file);

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(jobs_in),
		parse_dag_file(dag_in),
		parse_abort_file<dtime_t>(aborts_in),
		parse_tas_file<dtime_t>(shaper_in),
		1
	};

	Analysis_options opts;
	opts.early_exit = false;

	auto space = TSN::State_space<dtime_t>::explore(prob, opts);
	CHECK(space.is_schedulable());

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 6);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 7);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 16);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 17);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 26);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 27);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 21);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 24);
	CHECK(space.get_finish_times(prob.jobs[4]).min() == 53);
	CHECK(space.get_finish_times(prob.jobs[4]).max() == 57);

};

const std::string od_check_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"  1,   1,    0,    0,    1,    2,  10,    0\n";

const std::string od_check_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     1,   0,    2,  2, 0, 0\n";

TEST_CASE("[TSN] Check for Overlap Delete") {
	auto jobs_in = std::istringstream(od_check_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(od_check_shaper_file);

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(jobs_in),
		parse_dag_file(dag_in),
		parse_abort_file<dtime_t>(aborts_in),
		parse_tas_file<dtime_t>(shaper_in),
		1
	};

	Analysis_options opts;
	opts.early_exit = false;

	auto space = TSN::State_space<dtime_t>::explore(prob, opts);
	
	// typedef std::deque<Interval<Time>> Intervals;

	std::deque<Interval<dtime_t>> m1,rem1,res1;
	m1.emplace_back(Interval<dtime_t>{10, 20});
	rem1.emplace_back(Interval<dtime_t>{16, 18});
	res1.emplace_back(Interval<dtime_t>{10, 15});
	res1.emplace_back(Interval<dtime_t>{19, 20});
	CHECK(space.od(m1,rem1) == res1);

	std::deque<Interval<dtime_t>> m2,rem2,res2;
	m2.emplace_back(Interval<dtime_t>{10, 20});
	rem2.emplace_back(Interval<dtime_t>{16, 20});
	res2.emplace_back(Interval<dtime_t>{10, 15});
	CHECK(space.od(m2,rem2) == res2);

	std::deque<Interval<dtime_t>> m3,rem3,res3;
	m3.emplace_back(Interval<dtime_t>{10, 20});
	rem3.emplace_back(Interval<dtime_t>{10, 18});
	res3.emplace_back(Interval<dtime_t>{19, 20});
	CHECK(space.od(m3,rem3) == res3);

	std::deque<Interval<dtime_t>> m4,rem4,res4;
	m4.emplace_back(Interval<dtime_t>{10, 20});
	rem4.emplace_back(Interval<dtime_t>{7, 18});
	res4.emplace_back(Interval<dtime_t>{19, 20});
	CHECK(space.od(m4,rem4) == res4);

	std::deque<Interval<dtime_t>> m5,rem5,res5;
	m5.emplace_back(Interval<dtime_t>{10, 20});
	rem5.emplace_back(Interval<dtime_t>{16, 25});
	res5.emplace_back(Interval<dtime_t>{10, 15});
	CHECK(space.od(m5,rem5) == res5);

	std::deque<Interval<dtime_t>> m6,rem6,res6;
	m6.emplace_back(Interval<dtime_t>{10, 20});
	rem6.emplace_back(Interval<dtime_t>{10, 20});
	res6 = space.od(m6,rem6);
	CHECK(res6.size() == 0);
};