#include "doctest.h"

#include <iostream>
#include <sstream>

#include "problem.hpp"
#include "io.hpp"
#include "interval.hpp"
#include "tsn/space.hpp"

using namespace NP;

static const auto inf = Time_model::constants<dtime_t>::infinity();

// SET OF PACKET AND SHAPER FILES THAT CHECKS IF BRANCHING AND MERGING PROPERTIES 
// PERFORM CORRECTLY. DEADLINE MISSES STATES AND EXAMPLES WITH NO DEADLINE MISS
// ALSO CHECKED
// =========================================================================

// PACKET FILES
// -----------------------------------------------------------------------

// A packet set that when run with the gates "TAS_constGB0_shaper_file", 
// "TAS_varGB_ndm_shaper_file" and "TAS_varGB_dm_shaper_file" produces a 
// graph with no branches in it
const std::string TAS_nobranch_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"   1,  1,    0,    0,    1,    2,  10,    0\n"
"   1,  2,   10,   10,    1,    2,  20,    0\n"
"   1,  3,   20,   20,    1,    2,  30,    0\n"
"   2,  1,   10,   12,    5,    7,  30,    1\n"
"   3,  1,   25,   25,   15,   19,  60,    2\n";

// A packet set that when run with the gates "TAS_constGB0_shaper_file", 
// "TAS_varGB_ndm_shaper_file" and "TAS_varGB_dm_shaper_file" produces a 
// graph that shows branching and merging
const std::string TAS_branch_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"   1,  1,    0,    0,    1,    2,  10,    0\n"
"   1,  2,   10,   10,    1,    2,  20,    0\n"
"   1,  3,   20,   20,    1,    2,  40,    0\n"
"   1,  4,   25,   25,    1,    2,  40,    0\n"
"   2,  1,   10,   12,    5,    7,  30,    1\n"
"   2,  2,   20,   22,    5,    7,  50,    1\n"
"   3,  1,   25,   25,   15,   19,  60,    2\n";


// SHAPER FILES
// ----------------------------------------------------------------------

// A shaper file that supports TAS while having a constant guardband of 0
const std::string TAS_constGB0_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     0,   0,    2,  2, 0, 0\n"
"   1,  30,   0,   1,     0,   0,   11,  4\n"
"   2,  60,   0,   1,     0,   0,   20, 17\n";

// A shaper file that assumes a variable guardband that 
// produces a deadline miss when run with the packet set "TAS_nobranch_packets_file"
const std::string TAS_varGB_ndm_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     1,   0,    2,  2, 0, 0\n"
"   1,  30,   0,   1,     1,   0,   11,  4\n"
"   2,  60,   0,   1,     1,   0,   20, 17\n";

// A shaper file that assumes a variable guardband that 
// does not produce deadline miss when run with the packet
// set "TAS_nobranch_packets_file" 
const std::string TAS_varGB_dm_shaper_file =
"Prio, Per, CBS, TAS, isVar, CGB, Ints\n"
"   0,  10,   0,   1,     1,   0,    2,  2, 0, 0\n"
"   1,  30,   0,   1,     1,   0,   11,  2\n"
"   2,  60,   0,   1,     1,   0,   20, 17\n";


TEST_CASE("[TSN] Constant Guardband for TAS with no branching and IPG=0") {
	auto jobs_in = std::istringstream(TAS_nobranch_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_constGB0_shaper_file);

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

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 1);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 2);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 11);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 12);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 21);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 26);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 20);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 22);
	CHECK(space.get_finish_times(prob.jobs[4]).min() == 52);
	CHECK(space.get_finish_times(prob.jobs[4]).max() == 56);

};

TEST_CASE("[TSN] Constant Guardband for TAS with no branching and IPG=1") {
	auto jobs_in = std::istringstream(TAS_nobranch_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_constGB0_shaper_file);

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(jobs_in),
		parse_dag_file(dag_in),
		parse_abort_file<dtime_t>(aborts_in),
		parse_tas_file<dtime_t>(shaper_in),
		1
	};

	Analysis_options opts;
	opts.c_ipg = 3;
	opts.early_exit = false;

	auto space = TSN::State_space<dtime_t>::explore(prob, opts);
	CHECK(space.is_schedulable());

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 1);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 2);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 11);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 12);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 25);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 27);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 20);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 22);
	CHECK(space.get_finish_times(prob.jobs[4]).min() == 52);
	CHECK(space.get_finish_times(prob.jobs[4]).max() == 56);

};

TEST_CASE("[TSN] Constant Guardband for TAS with branching and IPG=0") {
	auto jobs_in = std::istringstream(TAS_branch_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_constGB0_shaper_file);

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

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 1);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 2);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 11);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 12);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 21);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 31);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 27);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 33);
	CHECK(space.get_finish_times(prob.jobs[4]).min() == 20);
	CHECK(space.get_finish_times(prob.jobs[4]).max() == 22);
	CHECK(space.get_finish_times(prob.jobs[5]).min() == 26);
	CHECK(space.get_finish_times(prob.jobs[5]).max() == 30);
	CHECK(space.get_finish_times(prob.jobs[6]).min() == 52);
	CHECK(space.get_finish_times(prob.jobs[6]).max() == 56);

};

TEST_CASE("[TSN] Variable Guardband for TAS with no deadline miss and IPG=0") {
	auto jobs_in = std::istringstream(TAS_nobranch_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_varGB_ndm_shaper_file);

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

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 5);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 6);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 15);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 16);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 25);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 26);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 20);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 23);
	CHECK(space.get_finish_times(prob.jobs[4]).min() == 52);
	CHECK(space.get_finish_times(prob.jobs[4]).max() == 56);

};

TEST_CASE("[TSN] Variable Guardband for TAS with deadline miss and IPG=0") {
	auto jobs_in = std::istringstream(TAS_nobranch_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_varGB_dm_shaper_file);

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

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 5);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 6);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 18);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 20);

};

// ================================================================================

// PACKET AND GATE SETS THAT CHECK IF THE DELETION OF ELIGIBLE TRANSMISSION
// INTERVAL BY HIGHER PRIORITY PACKETS FUNCTIONS AS INTENDED
// ================================================================================

// PACKET FILE
const std::string TAS_hp_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,    DL, Prio\n"
"   0,  0,   63,   66,    0,   63,  7063,    0\n"
"   0,  1, 7063, 7066,    0,   63, 14063,    0\n"
"   1,  0,  728,  764,    0,  728,  8728,    1\n"
"   1,  1, 8728, 8764,    0,  728, 16728,    1\n";

// SHAPER FILE
const std::string TAS_hp_shaper_file =
"Prio,  Per, CBS, TAS, isVar, CGB, Ints\n"
"   0, 1000,   0,   1,     1, 728,    0,  0,  89, 911\n"
"   1, 1000,   0,   1,     1, 728,    0, 89, 911,   0\n";

TEST_CASE("[TSN] Higher Prioirty Open Gate deletion from Eligible Transmission Interval") {
	auto jobs_in = std::istringstream(TAS_hp_packets_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");
	auto shaper_in = std::istringstream(TAS_hp_shaper_file);

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

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 1000);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 1063);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 8000);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 8063);
	CHECK(space.get_finish_times(prob.jobs[2]).min() == 1089);
	CHECK(space.get_finish_times(prob.jobs[2]).max() == 1817);
	CHECK(space.get_finish_times(prob.jobs[3]).min() == 9089);
	CHECK(space.get_finish_times(prob.jobs[3]).max() == 9817);

};

// ================================================================================

// PACKET AND GATE SETS THAT ARE USED AS DUMMY FILES IN ORDER TO CHECK THE 
// OVERLAP DELETE FUNCTION
// ================================================================================

// PACKET FILE
const std::string od_check_packets_file =
"TID, JID, Rmin, Rmax, Cmin, Cmax,  DL, Prio\n"
"  1,   1,    0,    0,    1,    2,  10,    0\n";

// SHAPER FILE
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