#ifndef NP_PROBLEM_HPP
#define NP_PROBLEM_HPP

#include "jobs.hpp"
#include "precedence.hpp"
#include "aborts.hpp"
#include "selfsuspending.hpp"

namespace NP {

	// Description of a non-preemptive scheduling problem
	template<class Time>
	struct Scheduling_problem {

		typedef typename Job<Time>::Job_set Workload;
		typedef typename std::vector<Abort_action<Time>> Abort_actions;
		typedef typename std::vector<Suspending_Task<Time>> Suspending_Tasks;

		// ** Description of the workload:
		// (1) a set of jobs
		Workload jobs;
		// (2) a set of precedence constraints among the jobs
		Precedence_constraints dag;
		// (3) abort actions for (some of) the jobs
		Abort_actions aborts;
		// (4) for self-suspending tasks
		Suspending_Tasks sts;

		// ** Platform model:
		// on how many (identical) processors are the jobs being
		// dispatched (globally, in priority order)
		unsigned int num_processors;

		// Classic default setup: no abort actions
		Scheduling_problem(Workload jobs, Precedence_constraints dag,
		                   unsigned int num_processors = 1)
		: num_processors(num_processors)
		, jobs(jobs)
		, dag(dag)
		{
			assert(num_processors > 0);
			validate_prec_refs<Time>(dag, jobs);
		}

		// For self-suspending tasks
		Scheduling_problem(Workload jobs, Suspending_Tasks sts,
		                   unsigned int num_processors = 1)
		: num_processors(num_processors)
		, jobs(jobs)
		, sts(sts)
		{
			assert(num_processors > 0);
			validate_susp_refs<Time>(sts, jobs);
		}

		// Constructor with abort actions and precedence constraints
		Scheduling_problem(Workload jobs, Precedence_constraints dag,
		                   Abort_actions aborts,
		                   unsigned int num_processors)
		: num_processors(num_processors)
		, jobs(jobs)
		, dag(dag)
		, aborts(aborts)
		{
			assert(num_processors > 0);
			validate_prec_refs<Time>(dag, jobs);
			validate_abort_refs<Time>(aborts, jobs);
		}

		// Full Constructor 
		Scheduling_problem(Workload jobs, Precedence_constraints dag,
						   Suspending_Tasks sts,
		                   Abort_actions aborts,
		                   unsigned int num_processors)
		: num_processors(num_processors)
		, jobs(jobs)
		, dag(dag)
		, sts(sts)
		, aborts(aborts)
		{
			assert(num_processors > 0);
			validate_prec_refs<Time>(dag, jobs);
			validate_abort_refs<Time>(aborts, jobs);
			validate_susp_refs<Time>(sts, jobs);
		}

		// Convenience constructor: no DAG, no abort actions
		Scheduling_problem(Workload jobs,
		                   unsigned int num_processors = 1)
		: jobs(jobs)
		, num_processors(num_processors)
		{
			assert(num_processors > 0);
		}
	};

	// Common options to pass to the analysis engines
	struct Analysis_options {
		// After how many seconds of CPU time should we give up?
		// Zero means unlimited.
		double timeout;

		// After how many scheduling decisions (i.e., depth of the
		// schedule graph) should we terminate the analysis?
		// Zero means unlimited.
		unsigned int max_depth;

		// Should we terminate the analysis upon encountering the first
		// deadline miss?
		bool early_exit;

		// Should we use state-merging techniques or naively explore the
		// whole state space in a brute-force manner (only useful as a
		// baseline).
		bool be_naive;

		// Implementation-specific: how large should the lookup table
		// of the main workload index be?
		std::size_t num_buckets;

		// If using self-suspending tasks, then which implementation to take
		unsigned int use_self_suspensions;

		// If using supernodes
		bool use_supernodes;

		Analysis_options()
		: timeout(0)
		, max_depth(0)
		, early_exit(true)
		, num_buckets(1000)
		, be_naive(false)
		, use_self_suspensions(0)
		, use_supernodes(true)
		{
		}
	};
}

#endif
