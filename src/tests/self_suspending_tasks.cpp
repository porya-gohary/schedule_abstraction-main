#include "doctest.h"

#include <iostream>
#include <sstream>

#include "io.hpp"
#include "uni/space.hpp"

using namespace NP;

#define NOSUSP 0
#define GENERAL_SUSP 1
#define PATHWISE_SUSP 2


const std::string uniproc_sn_gen_jobs_file =
"   Task ID,     Job ID,          Arrival min,          Arrival max,             Cost min,             Cost max,             Deadline,             Priority\n"
"1, 1,  1,  2, 3,  4,  8,  1\n"
"2, 1,  0,  0, 1,  1, 15,  2\n";


const std::string uniproc_sn_gen_susp_dag_file =
"Predecessor TID,	Predecessor JID,	Successor TID, Successor JID, Min_Sus, Max_Sus\n"
" 1,  1,  2,  1,  2,  5\n";

TEST_CASE("[susp] Uniproc Supernode General") {
	auto susp_dag_in = std::istringstream(uniproc_sn_gen_susp_dag_file);
	auto in = std::istringstream(uniproc_sn_gen_jobs_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(in),
		parse_suspending_file<dtime_t>(susp_dag_in)};

	Analysis_options opts;

	opts.be_naive = true;
	auto nspace = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(nspace.is_schedulable());

	opts.be_naive = false;
	auto space = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(space.is_schedulable());

	for (const Job<dtime_t>& j : prob.jobs) {
		CHECK(nspace.get_finish_times(j) == space.get_finish_times(j));
		CHECK(nspace.get_finish_times(j).from() != 0);
	}

	CHECK(space.get_finish_times(prob.jobs[0]).min() == 4);
	CHECK(space.get_finish_times(prob.jobs[0]).max() == 6);
	CHECK(space.get_finish_times(prob.jobs[1]).min() == 7);
	CHECK(space.get_finish_times(prob.jobs[1]).max() == 12);
}
