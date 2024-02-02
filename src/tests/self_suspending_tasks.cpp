#include "doctest.h"

#include <iostream>
#include <sstream>

#include "io.hpp"
#include "uni/space.hpp"

using namespace NP;

#define NOSUSP 0
#define GENERAL_SUSP 1
#define PATHWISE_SUSP 2

const std::string uniproc_sn_susp_jobs_file =
"Task ID, Job ID, Arrival min, Arrival max, Cost min, Cost max, Deadline, Priority\n"
"1      , 1      , 0          , 0          , 1        , 2        , 10      , 1      \n"
"1      , 2      , 10         , 10         , 1        , 2        , 20      , 2      \n"
"1      , 3      , 20         , 20         , 1        , 2        , 50      , 3      \n"
"1      , 4      , 30         , 30         , 1        , 2        , 60      , 4      \n"
"1      , 5      , 40         , 40         , 1        , 2        , 70      , 5      \n"
"1      , 6      , 50         , 50         , 1        , 2        , 80      , 6      \n"
"2      , 7      , 0          , 0          , 7        , 7        , 30      , 8      \n"
"2      , 8      , 30         , 30         , 7        , 7        , 100     , 9      \n"
"3      , 9      , 0          , 0          , 3        , 13       , 60      , 7      ";

const std::string uniproc_sn_susp_dag_file =
"Predecessor TID, Predecessor JID, Successor TID, Successor JID, Min_Sus, Max_Sus\n"
"2               , 7              , 1            , 3            , 10     , 12     ";

TEST_CASE("[susp] Uniproc Supernode Self-Suspensions") {
	auto susp_dag_in = std::istringstream(uniproc_sn_susp_dag_file);
	auto in = std::istringstream(uniproc_sn_susp_jobs_file);
	auto dag_in  = std::istringstream("\n");
	auto aborts_in = std::istringstream("\n");

	Scheduling_problem<dtime_t> prob{
		parse_file<dtime_t>(in),
		parse_suspending_file<dtime_t>(susp_dag_in)};

	Analysis_options opts;

	opts.use_self_suspensions = GENERAL_SUSP;

	opts.be_naive = true;
	auto gen_nspace = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(gen_nspace.is_schedulable());

	opts.be_naive = false;
	auto gen_space = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(gen_space.is_schedulable());

	opts.use_self_suspensions = PATHWISE_SUSP;

	opts.be_naive = true;
	auto pw_nspace = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(pw_nspace.is_schedulable());

	opts.be_naive = false;
	auto pw_space = Uniproc::State_space<dtime_t>::explore(prob, opts);
	CHECK(pw_space.is_schedulable());

	for (const Job<dtime_t>& j : prob.jobs) {
		CHECK(gen_nspace.get_finish_times(j) == gen_space.get_finish_times(j));
		CHECK(pw_nspace.get_finish_times(j) == pw_space.get_finish_times(j));
		CHECK(pw_nspace.get_finish_times(j).upto() <= gen_space.get_finish_times(j).upto());
	}
}
