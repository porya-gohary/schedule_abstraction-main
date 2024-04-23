#ifndef SCHEDULE_SPACE_H
#define SCHEDULE_SPACE_H

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>

#include <iostream>
#include <ostream>
#include <cassert>

#include "config.h"
#include "problem.hpp"
#include "jobs.hpp"
#include "precedence.hpp"
#include "selfsuspending.hpp"
#include "clock.hpp"

#include "uni/state.hpp"

#define NOSUSP 0
#define GENERAL_SUSP 1
#define PATHWISE_SUSP 2

namespace NP {

	namespace Uniproc {

		template<class Time> class State_space
		{
			public:

			typedef Scheduling_problem<Time> Problem;
			typedef typename Scheduling_problem<Time>::Workload Workload;
			typedef typename Scheduling_problem<Time>::Abort_actions Abort_actions;
			typedef typename Scheduling_problem<Time>::Suspending_Tasks Suspending_Tasks;
			typedef Schedule_state<Time> State;
			typedef Schedule_node<Time> Node;

			// This is the explore() that is called from the nptest file.
			// This function processes the Problem and based on the analysis options determines 
			// a State_space 's' that is returned. the cpu_time computes the time taken to run the analysis
			static State_space explore(
					const Problem& prob,
					const Analysis_options& opts)
			{
				// this is a uniprocessor analysis
				assert(prob.num_processors == 1);

				auto s = State_space(prob.jobs, prob.dag, prob.sts, prob.aborts,
									 opts.timeout, opts.max_depth,
									 opts.num_buckets, opts.early_exit, opts.use_self_suspensions, opts.use_supernodes);
				s.cpu_time.start();
				if (opts.be_naive)
					s.explore_naively();
				else
					s.explore();

				s.cpu_time.stop();
				return s;
			}

			// convenience interface for tests
			static State_space explore_naively(const Workload& jobs)
			{
				Problem p{jobs};
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static State_space explore(const Workload& jobs)
			{
				Problem p{jobs};
				Analysis_options o;
				return explore(p, o);
			}

			// get_finish_times searches rta (stores the response time of all the jobs) for job j 
			// and returns the finish time interval if j is found else returns the interval[0,infinity]
			Interval<Time> get_finish_times(const Job<Time>& j) const
			{
				auto rt_item = rta[j.get_job_index()];
				if (!rt_item.valid) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rt_item.rt;
				}
			}
		  
			Interval<Time> get_finish_times(Job_index j) const
			{
				auto rt_item = rta[j];
				if (!rt_item.valid) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rt_item.rt;
				}
			}

			bool is_schedulable() const
			{
				return !aborted && !observed_deadline_miss;
			}

			bool was_timed_out() const
			{
				return timed_out;
			}

			unsigned long number_of_nodes() const
			{
				return num_nodes;
			}

			unsigned long number_of_states() const
			{
				return num_states;
			}

			unsigned long number_of_edges() const
			{
				return num_edges;
			}

			unsigned long max_exploration_front_width() const
			{
				return width;
			}

			double get_cpu_time() const
			{
				return cpu_time;
			}

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH

			// The structure of an edge and all edge related functionalities are only run if CONFIG_COLLECT_SCHEDULE_GRAPH
			// is set. This is because in the analysis edge has no meaning or importance. 
			// There is no information stored in an edge that is necessary for building the graph or furthering the analysis.
			// It is only required when a DOT file is to be generated so that the graphical version of the
			// abstraction graph can be represented with edges.
		  // RV:  replaced State with Node
			struct Edge {
				const Job<Time>* scheduled;
				const Node* source;  
				const Node* target;
				const Interval<Time> finish_range;

				Edge(const Job<Time>* s, const Node* src, const Node* tgt,
					 const Interval<Time>& fr)
				: scheduled(s)
				, source(src)
				, target(tgt)
				, finish_range(fr)
				{
				}

				bool deadline_miss_possible() const
				{
					return scheduled->exceeds_deadline(finish_range.upto());
				}

				Time earliest_finish_time() const
				{
					return finish_range.from();
				}

				Time latest_finish_time() const
				{
					return finish_range.upto();
				}

				Time earliest_start_time() const
				{
					return finish_range.from() - scheduled->least_cost();
				}

				Time latest_start_time() const
				{
					return finish_range.upto() - scheduled->maximal_cost();
				}

			};

			const std::deque<Edge>& get_edges() const
			{
				return edges;
			}

		  // RV: replaced State with Node
			const std::deque<Node>& get_nodes() const
			{
				return nodes;
			}

#endif
			private:

			typedef Job_set Scheduled;

			// All the nodes created through the making of the schedule abstrcation graph are stored in Nodes.
			// This is done in the function new_node()
			typedef std::deque<Node> Nodes;

			typedef Node* Node_ref;
			typedef State* State_ref;

			// The Nodes_map typedef allows for storing a key-pair of hash_value_t and Node_ref. The hash_value_t is the
			// typedef defined in jobs.hpp that allows for creating a unique key for each node. A variable called the
			// nodes_by_key is created using Nodes_map. When a new node is made it is added to the map.
			// The use of Nodes_map can be found in new_node(), done_with_current_state() ad schedule().
		  // RV: unordered_multimap might not be optimal. Some sorted list of tree structure might be useful.
			typedef std::unordered_multimap<hash_value_t, Node_ref> Nodes_map;

			typedef const Job<Time>* Job_ref;

			// By_time_map stores all the jobs along with a particular time. there are three variables of type By_time_map.
			// The three variables jobs_by_latest_arrival, jobs_by_earliest_arrival and jobs_by_deadline each store the job
			// along with the latest arrival time, earliest arrival time and the deadline respectively.
			typedef std::multimap<Time, Job_ref> By_time_map;

			// The Todo_queue stores the references of all the nodes that are edge nodes of the graph that are yet to
			// be processed to check if the graph can be built further. There are usually 2 Todo_queues defined. This is 
			// because the graph is built in a breadth first manner and so the one todo_queue is for the new breadth that 
			// is made while building the graph and the other todo_queue stores the nodes left to be processed in the 
			// current breadth.
			typedef std::deque<Node_ref> Todo_queue;

			// Response_times stores and updates the finish time intervals of the jobs as the graph is being built
			// RV: Response_times was there before supernodes. How does this work in global?
			//     global/space.hpp uses JobID here. Not sure whether Job_index would work better.
			//     The rta structure will contain the response times of all jobs.
			//     Instead of using a map with JobID and search, a vector can be used with Job_index as index.
			struct Response_time_item  {
				bool valid;  // RV: is rt valid for the given Job_index?
				Interval<Time> rt;

				Response_time_item()
					: valid(false)
					, rt(0,0)
				{
				}
			};
			typedef std::vector<Response_time_item> Response_times;

			typedef std::vector<std::size_t> Job_precedence_set;


#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			typedef std::deque<Edge> Edges;
			typedef typename std::deque<Edge>::iterator Edge_ref;

			Edges edges;
#endif

			Response_times rta;
			bool aborted;
			bool timed_out;

			const Workload& jobs;

			// RV: in global, these maps have a const version, to prevent accidental changes.
			//     Added here as well.
			std::vector<Job_precedence_set> _job_precedence_sets;

			By_time_map _jobs_by_latest_arrival_with_susp;
			By_time_map _jobs_by_latest_arrival_without_susp;		  
			By_time_map _jobs_by_earliest_arrival;
			By_time_map _jobs_by_deadline;
			// RV: global/space.hpp also has jobs_by_win .
			//     global_space.hpp doesn't have abort_actions. 
			std::vector<const Abort_action<Time>*> abort_actions;

			const By_time_map& jobs_by_latest_arrival_with_susp;
			const By_time_map& jobs_by_latest_arrival_without_susp;
			const By_time_map& jobs_by_earliest_arrival;
			const By_time_map& jobs_by_deadline;
			const std::vector<Job_precedence_set>& job_precedence_sets;


			// #NS# The following variables are used for self-suspending tasks, and this is a very very messy way to store things 
			// #NS# but it works for now.
			
			// In order to store the componenets of the self-suspending tasks, we create a vector of a vector of references
			// predecessors_of: for each job, constaints a list of predeccessors 
			// from/_suspending_tasks: for each job, contains a list of successors and the suspension time until each successor becomes ready
			// RV: the suspending_tasks_list might include Job_index, next to JobID.
			//     the predecessors_of and successors_of might have the const protection.
			typedef std::vector<const Suspending_Task<Time>*> Predecessors_list;
			typedef std::vector<std::pair<Job_ref, Interval<Time>>> Successors_list;
			
		public:
			typedef std::vector<Predecessors_list> Predecessors;
			typedef std::vector<Successors_list> Successors;

		private:
			Predecessors _predecessors_of;
			Successors _successors_of;

			const Predecessors& predecessors_of;
			const Successors& successors_of;


			Nodes nodes;
			unsigned long num_nodes, num_states, num_edges, width;
			Nodes_map nodes_by_key;

			static const std::size_t num_todo_queues = 3;

			Todo_queue todo[num_todo_queues];
			int todo_idx;
			unsigned long current_job_count;

			Processor_clock cpu_time;
			double timeout;

			unsigned int max_depth;

			bool early_exit;
			bool observed_deadline_miss;
			unsigned int want_self_suspensions;
			bool use_supernodes;

			// Constructor of the State_space class
			State_space(const Workload& jobs,
						const Precedence_constraints &dag_edges,
						const Suspending_Tasks &susps,
						const Abort_actions& aborts,
						double max_cpu_time = 0,
						unsigned int max_depth = 0,
						std::size_t num_buckets = 1000,
						bool early_exit = true,
						unsigned int use_self_suspensions = 0,
						bool use_supernodes = false)
			: jobs(jobs)
			, aborted(false)
			, timed_out(false)
			, timeout(max_cpu_time)
			, max_depth(max_depth)
			, num_nodes(0)
			, num_states(0)
			, num_edges(0)
			, width(0)
			, todo_idx(0)
			, current_job_count(0)
			, rta(jobs.size())
			, jobs_by_latest_arrival_with_susp(_jobs_by_latest_arrival_with_susp)
			, jobs_by_latest_arrival_without_susp(_jobs_by_latest_arrival_without_susp)
			, jobs_by_earliest_arrival(_jobs_by_earliest_arrival)
			, jobs_by_deadline(_jobs_by_deadline)
			, _job_precedence_sets(jobs.size())
			, _predecessors_of(jobs.size())
			, _successors_of(jobs.size())
			, job_precedence_sets(_job_precedence_sets)
			, predecessors_of(_predecessors_of)
			, successors_of(_successors_of)
			, early_exit(early_exit)
			, want_self_suspensions(use_self_suspensions)
			, observed_deadline_miss(false)
			, abort_actions(jobs.size(), NULL)
			, use_supernodes(use_supernodes)
			{
				for (auto e : dag_edges) {
					const Job<Time>& from = lookup<Time>(jobs, e.first);
					const Job<Time>& to   = lookup<Time>(jobs, e.second);
					_job_precedence_sets[to.get_job_index()].push_back(from.get_job_index());
				}
				for (const Suspending_Task<Time>& st : susps) {
					_predecessors_of[st.get_toIndex()].push_back(&st);
					_job_precedence_sets[st.get_toIndex()].push_back(st.get_fromIndex());
					_successors_of[st.get_fromIndex()].push_back({ &jobs[st.get_toIndex()], st.get_suspension() });
				}
				for (const Job<Time>& j : jobs) {
					if (_predecessors_of[j.get_job_index()].size()>0) {
						_jobs_by_latest_arrival_with_susp.insert({j.latest_arrival(), &j});
					} else {
						_jobs_by_latest_arrival_without_susp.insert({j.latest_arrival(), &j});
					}
					_jobs_by_earliest_arrival.insert({j.earliest_arrival(), &j});
					_jobs_by_deadline.insert({j.get_deadline(), &j});
				}
				for (const Abort_action<Time>& a : aborts) {
					const Job<Time>& j = lookup<Time>(jobs, a.get_id());
					abort_actions[j.get_job_index()] = &a;
				}
			}

			private:

		  // Not used
			// For all the jobs of Workload, the maximum deadline amongst all the jobs is returned
			static Time max_deadline(const Workload &jobs)
			{
				Time dl = 0;
				for (auto j : jobs)
					dl = std::max(dl, j.get_deadline());
				return dl;
			}

			// When a path in the graph has a job scheduled on it, the job and the finish time interval of the job on that 
			// path are passed as parameters to update_finish_times. If this job already exists in rts, then the finish time
			// stored in rta is widened with the new finish time 'range'. If no such job is present in rta, the job is added
			// to rta along with its finish time interval
			void update_finish_times(const Job<Time>& j, Interval<Time> range)
			{
				Job_index jidx = j.get_job_index();
				Response_time_item* rt_item = &(rta[jidx]);
				if (!rt_item->valid) {
					rt_item->valid=true;
					rt_item->rt = range;
					if (j.exceeds_deadline(range.upto()))
						observed_deadline_miss = true;
				} else {
					rt_item->rt.widen(range);
					if (j.exceeds_deadline(rt_item->rt.upto()))
						observed_deadline_miss = true;
				}
				DM("      New finish time range for " << j
				   << ": " << rt_item->rt << std::endl);

				if (early_exit && observed_deadline_miss)
					aborted = true;
			}

		  	// Old: the address of j is subtracted from the first element of jobs to find the distance between the two ie., index
			// New: job has an attribute for job_index, which is the position within the workload file.
		  // not used
			std::size_t index_of(const Job<Time>& j) const
			{
				// return (std::size_t) (&j - &(jobs[0]));
				return j.get_job_index();
			}

			// incomplete() checks to see that a job j jas not been scheduled yet.
			bool incomplete(const Scheduled &scheduled, const Job<Time>& j) const
			{
				return !scheduled.contains(j.get_job_index());
			}

			bool incomplete(const Scheduled &scheduled, Job_index j) const
			{
				return !scheduled.contains(j);
			}

			bool incomplete(const Node& n, const Job<Time>& j) const
			{
				return incomplete(n.get_scheduled_jobs(), j);
			}

			bool incomplete(const Node& n, Job_index j) const
			{
				return incomplete(n.get_scheduled_jobs(), jobs[j]);
			}

			// find next time by which a job is certainly ready in system state 's'
			Time next_certain_job_ready_time(const Node& n, const State &s)
			{
				Time ncjr_wos = n.get_next_certain_source_job_release();
				Time ncjr_ws = s.get_earliest_certain_successor_jobs_ready_time();
				return std::min(ncjr_wos, ncjr_ws);
			}

			// Find next time by which a source job (i.e., a job without predecessors) of higher priority than the reference_job
			// is certainly released in any state in the node 'n'. 
			Time next_certain_higher_priority_source_job_release(
				const Node& n,
				const Job<Time>& reference_job,
				Time until = Time_model::constants<Time>::infinity())
			{
				Time nejr = until;

				for (auto it = jobs_by_latest_arrival_without_susp.lower_bound(n.get_next_certain_source_job_release());
					 it != jobs_by_latest_arrival_without_susp.end(); it++) {
					const Job<Time>& j = *(it->second);

					if (nejr < j.latest_arrival())
						break;

					// irrelevant if not of higher priority
					if (!j.higher_priority_than(reference_job))
						continue;

					// not relevant if already scheduled
					if (!incomplete(n, j))
						continue;

					nejr = j.latest_arrival();
					// Jobs ordered by latest_arrival, so next jobs are later. We can thus stop searching.
					break;
				}
				return nejr;
			}

			// Find next time by which a successor job (i.e., a job with predecessors) of higher priority than the reference_job
			// is certainly released in system state 's' at or before a time 'max'.
			Time next_certain_higher_priority_successor_job_ready_time(
				const Node& n,
				const State& s,
				const Job<Time>& reference_job,
				Time until =Time_model::constants<Time>::infinity())
			{
				Time nejr = until;

				// a higer priority successor job cannot be ready before 
				// a job of any priority is released
				Time t_earliest = n.earliest_job_release();
				for (auto it = jobs_by_latest_arrival_with_susp.lower_bound(t_earliest);
					 it != jobs_by_latest_arrival_with_susp.end(); it++) 
				{
					const Job<Time>& j = *(it->second);

					if(nejr < j.latest_arrival())
						break;

					// irrelevant if not of higher priority
					if (!j.higher_priority_than(reference_job))
						continue;

					// not relevant if already scheduled or not ready (i.e., with non-completed predecessors)
					if (!incomplete(n, j) || !ready(n, j))
						continue;

					auto t = latest_ready_time_hp(n, s, j);
					
					if(nejr > t)
						nejr = t;
						// No break, as later jobs might have less suspension.
				}
				return nejr;
			}

			// The combined version, where optimization is not possible
			Time next_certain_higher_priority_job_ready_time(
				const Node& n,
				const State& s,
				const Job<Time>& reference_job,
				Time until = Time_model::constants<Time>::infinity())
			{
				Time nchpjr_wos = next_certain_higher_priority_source_job_release(n,reference_job, until);
				Time nchpjr_ws = next_certain_higher_priority_successor_job_ready_time(n,s,reference_job, nchpjr_wos);
				return nchpjr_ws;
			}

		  
// define a couple of iteration helpers

// For all the jobs that arrives after the earliest job release of a node n, check if these jobs are still incomplete 
#define foreach_possibly_pending_job(ppj_macro_local_n, ppj_macro_local_j) 	\
	for (auto ppj_macro_local_it = jobs_by_earliest_arrival			\
					 .lower_bound((ppj_macro_local_n).earliest_job_release()); \
		 ppj_macro_local_it != jobs_by_earliest_arrival.end() 				\
			&& (ppj_macro_local_j = ppj_macro_local_it->second); 	\
		 ppj_macro_local_it++) \
		if (incomplete(ppj_macro_local_n, *ppj_macro_local_j))

// Iterate over all incomplete jobs that are released no later than ppju_macro_local_until
// For all the jobs that arrive after the earliest job release of a node n, and before the ppju_macro_local_until time,
// check if they are incomplete
#define foreach_possbly_pending_job_until(ppju_macro_local_n, ppju_macro_local_j, ppju_macro_local_until) 	\
	for (auto ppju_macro_local_it = jobs_by_earliest_arrival			\
					 .lower_bound((ppju_macro_local_n).earliest_job_release()); \
		 ppju_macro_local_it != jobs_by_earliest_arrival.end() 				\
			&& (ppju_macro_local_j = ppju_macro_local_it->second, ppju_macro_local_j->earliest_arrival() <= (ppju_macro_local_until)); 	\
		 ppju_macro_local_it++) \
		if (incomplete(ppju_macro_local_n, *ppju_macro_local_j))

// Iterare over all incomplete jobs that are certainly released no later than
// cpju_macro_local_until remove state HERE 
// for_each_possibly_pending_job also check if their latest arrival occurs before the cpju_macro_local_until
#define  foreach_certainly_pending_job_until(cpju_macro_local_n, cpju_macro_local_j, cpju_macro_local_until) \
	foreach_possbly_pending_job_until(cpju_macro_local_n, cpju_macro_local_j, (cpju_macro_local_until)) \
		if (cpju_macro_local_j->latest_arrival() <= (cpju_macro_local_until))

			// returns true if there is certainly some pending job of higher
			// priority at the given time ready to be scheduled		  
			bool exists_certainly_released_higher_prio_job(
				const Node& n,
				const State& s,
				const Job<Time>& reference_job,
				Time at)
			{
				auto ts_min = s.earliest_finish_time();
				assert(at >= ts_min);

				DM("    ? " << __FUNCTION__ << " at " << at << ": "
				   << reference_job << std::endl);
				// consider all possibly pending jobs and check if they have a higher priority than the 
				// reference job while also having no predecessor jobs pending.
				const Job<Time>* jp;
				foreach_certainly_pending_job_until(n, jp, at) {
					const Job<Time>& j = *jp;
					DM("        - considering " << j << std::endl);
					// skip reference job
					if (&j == &reference_job)
						continue;
					// check priority
					if (!j.higher_priority_than(reference_job))
						continue;
					
					if (ready_at_hp(n, s, j, at)) {
						DM("          => found one: " << j << " <<HP<< "
						   << reference_job << std::endl);
						return true;
					}
				}
				return false;
			}

			// Find the earliest possible job release of all jobs in a node except for the ignored job
			Time earliest_possible_job_release(
				const Node& n,
				const Job<Time>& ignored_job)
			{
				DM("      - looking for earliest possible job release starting from: "
				   << n.earliest_job_release() << std::endl);
				const Job<Time>* jp;
				foreach_possibly_pending_job(n, jp) {
					const Job<Time>& j = *jp;

					DM("         * looking at " << j << std::endl);

					// skip if it is the one we're ignoring
					if (&j == &ignored_job)  //RV: is this the corrrect test?
						continue;

					DM("         * found it: " << j.earliest_arrival() << std::endl);

					// it's incomplete and not ignored => found the earliest
					return j.earliest_arrival();
				}

				DM("         * No more future releases" << std::endl);

				return Time_model::constants<Time>::infinity();
			}

			// Find the earliest possible certain job release of all source jobs (i.e., without predecessors) 
			// in a node except for the ignored job
			Time earliest_certain_source_job_release(
				const Node& n,
				const Job<Time>& ignored_job)
			{
				DM("      - looking for earliest certain source job release starting from: "
					<< n.get_next_certain_source_job_release() << std::endl);
				
				for (auto it = jobs_by_latest_arrival_without_susp.lower_bound(n.get_next_certain_source_job_release()); it != jobs_by_latest_arrival_without_susp.end(); it++) 
				{
					const Job<Time>* jp = it->second;
					DM("         * looking at " << *jp << std::endl);

					// skip if it is the one we're ignoring or the job was dispatched already
					if (jp == &ignored_job || !incomplete(n, *jp))  //RV: is this the corrrect test?
						continue;

					DM("         * found it: " << jp->latest_arrival() << std::endl);

					// it's incomplete and not ignored => found the earliest
					return jp->latest_arrival();
				}

				DM("         * No more future releases" << std::endl);
				return Time_model::constants<Time>::infinity();
			}

			// Ensure that there is no higher prioirty job that has certainly released
			bool priority_eligible(const Node& n, const State &s, const Job<Time> &j, Time t)
			{
				return !exists_certainly_released_higher_prio_job(n, s, j, t);
			}

			//checks if the job j can be the next potential job for a state 
			bool potentially_next(const Node &n, const State &s, const Job<Time> &j)
			{
				auto t_latest = s.latest_finish_time();

				// if t_latest >=  j.earliest_arrival(), then the
				// job is trivially potentially next, so check the other case.

				if (t_latest < j.earliest_arrival()) {
					Time r = next_certain_job_ready_time(n, s);
					// if something else is certainly released before j, then j can't
					// possibly be next
					if (r < j.earliest_arrival())
						return false;
				}
				return true;
			}

			// returns true if all predecessors of j have completed in state s
			bool ready(const Node &n, const Job<Time> &j)
			{
				const Job_precedence_set &preds =
					job_precedence_sets[j.get_job_index()];
				// check that all predecessors have completed
				return n.get_scheduled_jobs().includes(preds);
			}

			// Returns true if all predecessors of job 'j' have completed in state 's' no later than 'at'
			// assumes the job to check is a higher priority job than the job we try to disptach
			bool ready_at_hp(const Node &n, const State &s, const Job<Time> &j, Time at)
			{
				if (!ready(n, j))
					return false;
				
				// if job 'j' is ready in state 's', we check whether it will be ready no later than 'at'
				if (latest_ready_time_hp(n,s,j) > at)
					return false;
				else
					return true;
			}

			// Checks if job 'j' may be scheduled next in system state 's'
			bool is_eligible(const Node &n, const State &s, const Job<Time> &j)
			{
				// has been scheduled already, then false
				if (!incomplete(n, j)) {
					DM("  --> already complete"  << std::endl);
					return false;
				}

				// has predecessors yet to be scheduled, then false
				if (!ready(n, j)) {
					DM("  --> not ready"  << std::endl);
					return false;
				}

				// Job j is not the highest priority job to be scheduled, then false
				auto t_s = next_earliest_start_time(n, s, j);
				if (!priority_eligible(n, s, j, t_s)) {
					DM("  --> not prio eligible"  << std::endl);
					return false;
				}

				// if not potentially next, then false
				if (!potentially_next(n, s, j)) {
					DM("  --> not potentially next" <<  std::endl);
					return false;
				}

				// passes all the checks above then true
				return true;
			}

			void make_initial_node()
			{
				// construct initial node
				Node& n = new_node(jobs_by_earliest_arrival.begin()->first, jobs_by_latest_arrival_without_susp.begin()->first);
				State& s = new_state();
				n.add_state(&s);
			}

			// create a new node by adding an elemen to nodes, creating a new state and adding this new state to the node,
			// add this new node to the todo queue as it now a leaf node in the graph, yet to be processed. Add the node
			// along with its key to nodes_by_key. 
			template <typename... Args>
			Node& new_node(Args&&... args)
			{
				nodes.emplace_back(std::forward<Args>(args)...);
				Node_ref n_ref = &(*(--nodes.end()));
				auto njobs = n_ref->get_scheduled_jobs().size();
				assert (
					(!njobs && num_nodes == 0) // initial state
				    || (njobs == current_job_count + 1) // normal State
				    || (njobs == current_job_count + 2 && aborted) // deadline miss
				);
				auto idx = njobs % num_todo_queues;
				todo[idx].push_back(n_ref);
				nodes_by_key.insert(std::make_pair(n_ref->get_key(), n_ref));
				num_nodes++;
				width = std::max(width, (unsigned long) todo[idx].size() - 1);
				return *n_ref;
			}
			
			template <typename... Args>
			Node& new_node(Interval<Time> ftimes, const Node& from_node, const State& from, const Job<Time>& sched_job, Args&&... args)
			{
				nodes.emplace_back(from_node, sched_job, std::forward<Args>(args)...);
				Node_ref n_ref = &(*(--nodes.end()));

				auto njobs = n_ref->get_scheduled_jobs().size();
				assert (
					(!njobs && num_nodes == 0) // initial state
				    || (njobs == current_job_count + 1) // normal State
				    || (njobs == current_job_count + 2 && aborted) // deadline miss
				);
				auto idx = njobs % num_todo_queues;
				todo[idx].push_back(n_ref);
				nodes_by_key.insert(std::make_pair(n_ref->get_key(), n_ref));
				num_nodes++;
				width = std::max(width, (unsigned long) todo[idx].size() - 1);
				return *n_ref;
			}

			// Create a new state, either the first intial state with no parameters or a new state in the graph with 
			// finish times ftimes. States created by adding an element to states
			State& new_state()
			{
				State* s = new State();
				num_states++;
				return *s;
			}

			State& new_state(Interval<Time>& ftimes, const State& from, const Job<Time>& sched_job, const Job_set& scheduled_jobs)
			{
				State* s = new State(from, sched_job.get_job_index(), ftimes, scheduled_jobs,successors_of);
				num_states++;
				return *s;
			}

			// The todo queue contains the nodes that have to be processed in order to build the graph further,
			// there are currently two todo queues being used, one is for the nodes that have 'current_job_count' 
			// scheduled jobs and the other queue is for nodes that have 'current_job_count+1' scheduled jobs. 
			// In the function below, the first queue with 'current_job_count' scheduled jobs is checked to see if 
			// it is empty or not. If it is empty then the second queue with 'current_job_count++' scheduled jobs 
			// is checked to see if it empty or not.If both queues are empty then all nodes have been processed and 
			// the function returns false

			bool not_done()
			{
				// if the curent queue is empty, move on to the next
				if (todo[todo_idx].empty()) {
					current_job_count++;
					todo_idx = current_job_count % num_todo_queues;
					return !todo[todo_idx].empty();
				} else
					return true;
			}

			const Node& next_node()
			{
				auto n = todo[todo_idx].front();
				return *n;
			}

			bool in_todo(Node_ref n)
			{
				for (auto it : todo)
					if (it == next_earliest_job_abortion)
						return true;
				return false;
			}

			void check_cpu_timeout()
			{
				if (timeout && get_cpu_time() > timeout) {
					aborted = true;
					timed_out = true;
				}
			}

			void check_depth_abort()
			{
				if (max_depth && current_job_count == max_depth
					&& todo[todo_idx].empty()) {
					aborted = true;
				}
			}

			// If the a node has been processed, it can be removed from the todo queue as the todo queue only contains
			// nodes that are yet to be processed.
			void done_with_current_node()
			{
				Node_ref n = todo[todo_idx].front();
				// remove from TODO list
				todo[todo_idx].pop_front();

#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
				// If we don't need to collect all states, we can remove
				// all those that we are done with, which saves a lot of
				// memory.

				// remove from lookup map
				auto matches = nodes_by_key.equal_range(n->get_key());
				bool deleted = false;
				for (auto it = matches.first; it != matches.second; it++)
					if (it->second == n) {
						nodes_by_key.erase(it);
						deleted = true;
						break;
					}
				assert(deleted);

				// delete from master sequence to free up memory
				assert(&(*(nodes.begin())) == n);
				nodes.pop_front();
#endif
			}

			//////////////////////////////////////
			// Rules for finding the next state //
			//////////////////////////////////////

			// Returns the earliest time job 'j' may be ready in state 's' 
			// accounting for precedence and suspension time
			Time earliest_ready_time(const Node &n, const State &s, const Job<Time>& j)
			{
				Time seft = j.earliest_arrival();
				assert(ready(n, j));

				for(auto e : predecessors_of[j.get_job_index()])
				{
					Interval<Time> rbounds = (want_self_suspensions == PATHWISE_SUSP) ?
						s.get_pathwisejob_ft(e->get_fromIndex()) :
						get_finish_times(jobs[e->get_fromIndex()]);

					if(seft <= (rbounds.from() + e->get_minsus()))
						seft = rbounds.from() + e->get_minsus();
				}
				return seft;
			}

			// Returns the latest time job 'j' may be ready in state 's'
			// accounting for precedence and suspension time
			Time latest_ready_time(const Node &n, const State &s, const Job<Time>& j)
			{
				Time slft = j.latest_arrival();
				Time avail_max = s.latest_finish_time();

				for(auto e : predecessors_of[j.get_job_index()])
				{
					Time susp_max = e->get_maxsus();
					// if there is no suspension, and the predecessor has been dispatched already,
					// then the precedence constraint is certainly resolved when the processor become available
					if (susp_max == 0)
						slft = std::max(slft, avail_max);
					else
					{
						Interval<Time> rbounds = (want_self_suspensions == PATHWISE_SUSP) ?
							s.get_pathwisejob_ft(e->get_fromIndex()) :
							get_finish_times(e->get_fromIndex());

						slft = std::max(slft, rbounds.until() + susp_max);
					}
				}
				return slft;
			}

			// assumes the job to check is a higher priority job than the job we try to disptach
			Time latest_ready_time_hp(const Node& n, const State& s, const Job<Time>& j)
			{
				Time slft = j.latest_arrival();
				Time avail_min = s.earliest_finish_time();

				for (auto e : predecessors_of[j.get_job_index()])
				{
					Time susp_max = e->get_maxsus();
					// if there is no suspension, and the predecessor has been dispatched already,
					// then the precedence constraint is certainly resolved when the processor become available
					if (susp_max == 0)
						slft = std::max(slft, avail_min);
					else
					{
						Interval<Time> rbounds = (want_self_suspensions == PATHWISE_SUSP) ?
							s.get_pathwisejob_ft(e->get_fromIndex()) :
							get_finish_times(e->get_fromIndex());

						slft = std::max(slft, rbounds.until() + susp_max);
					}
				}
				return slft;
			}

			Time next_earliest_start_time(const Node &n, const State &s, const Job<Time>& j)
			{
				// t_S in paper, see definition 6.
				return std::max(s.earliest_finish_time(), earliest_ready_time(n, s, j));
			}

			Time next_earliest_finish_time(const Job<Time>& j, Time est)
			{
				// e_k, equation 5
				return est + j.least_cost();
			}

			Time next_latest_finish_time(const Job<Time>& j, Time lst)
			{
				return lst + j.maximal_cost();
			}

			Time next_earliest_job_abortion(const Abort_action<Time> &a)
			{
				return a.earliest_trigger_time() + a.least_cleanup_cost();
			}

			Time next_latest_job_abortion(const Abort_action<Time> &a)
			{
				return a.latest_trigger_time() + a.maximum_cleanup_cost();
			}

			Interval<Time> next_finish_times(const Node &n, const State &s, const Job<Time> &j, Time est, Time lst)
			{
				// #NS# This function I messed with again, it takes in lst as a parameter which I don't think it should
				// maybe a better way to write this funciton?
				auto i = j.get_job_index();

				if (abort_actions[i]) {
					// complicated case -- need to take aborts into account

					auto et = abort_actions[i]->earliest_trigger_time();

					// Rule: if we're certainly past the trigger, the job is
					//       completely skipped.

					if (s.earliest_finish_time() >= et)
						// job doesn't even start, is skipped immediately
						return s.finish_range();

					// Otherwise, it might start execution. Let's compute the
					// regular and aborted completion times.

					auto eft = next_earliest_finish_time(j, est);
					auto lft = next_latest_finish_time(j, lst);

					auto eat = next_earliest_job_abortion(*abort_actions[i]);
					auto lat = next_latest_job_abortion(*abort_actions[i]);

					return Interval<Time>{
						std::min(eft, eat),
						std::min(lft, lat)
					};
				}
				else {
					// standard case -- this job is never aborted or skipped
					return Interval<Time>{
						next_earliest_finish_time(j, est),
						next_latest_finish_time(j, lst)
					};
				}
			}

			// creates a new edge from an existing to another node that can either be an existing node
			// in case of merge or a new node when merge is not possible
			void process_new_edge(
				const Node& from,
				const Node& to,
				const Job<Time>& j,
				const Interval<Time>& finish_range)
			{
				// update response times
				update_finish_times(j, finish_range);
				// update statistics
				num_edges++;
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
				edges.emplace_back(&j, &from, &to, finish_range);
#endif
			}

			// naive: no state merging
			void schedule_job(const Node &n, const State& s, const Job<Time> &j)
			{
				Time t_high = next_certain_higher_priority_job_ready_time(n, s, j);
				Time t_wc = std::max(s.latest_finish_time(), next_certain_job_ready_time(n, s));
				Time lst = std::min(t_wc, t_high-Time_model::constants<Time>::epsilon());
				Time est = next_earliest_start_time(n, s, j);
								
				Interval<Time> finish_range = next_finish_times(n, s, j, est, lst);
				Node& next =
					new_node(finish_range, n, s, j, j.get_job_index(),
							  earliest_possible_job_release(n, j),
							  earliest_certain_source_job_release(n, j));
				/*DM("      -----> N" << (nodes.end() - nodes.begin())
				   << std::endl);*/
				State& next_state = new_state(finish_range, s, j, next.get_scheduled_jobs());
				next.add_state(&next_state);
				process_new_edge(n, next, j, next_state.finish_range());
			}

			// Start by making the initial node. While the graph has not finished building and is not aborted, we check to
			// see how the states within each node is manipulated
			void explore_naively()
			{
				make_initial_node();

				while (not_done() && !aborted) {
					const Node& n = next_node();

					DM("\n==================================================="
					   << std::endl);
					/*DM("Looking at: N"
					   << (todo[todo_idx].front() - nodes.begin() + 1)
					   << " " << n << std::endl);*/
					const auto *n_states = n.get_states();

					// for each state in the node n
					for(State *s : *n_states)
					{
						// Identify relevant interval for next job
						// relevant job buckets
						auto ts_min = s->earliest_finish_time();
						auto rel_min = n.earliest_job_release();
						auto t_l = std::max(next_certain_job_ready_time(n, *s), s->latest_finish_time());

						Interval<Time> next_range{std::min(ts_min, rel_min), t_l};

						DM("=> next range = " << next_range << std::endl);

						bool found_at_least_one = false;

						DM("\n---\nChecking for pending and later-released jobs:"
						   << std::endl);
						const Job<Time>* jp;
						foreach_possbly_pending_job_until(n, jp, next_range.upto()) {
							const Job<Time>& j = *jp;
							DM("+ " << j << std::endl);
							// if it can be scheduled next...
							if (is_eligible(n, *s, j)) {
								DM("  --> can be next "  << std::endl);
								// create the relevant state and continue
								schedule_job(n, *s, j);
								found_at_least_one = true;
							}
						}

						DM("---\nDone iterating over all jobs." << std::endl);

						// check for a dead end
						if (!found_at_least_one &&
							n.get_scheduled_jobs().size() != jobs.size()) {
							// out of options and we didn't schedule all jobs
							observed_deadline_miss = true;
							if (early_exit)
								aborted = true;
						}

						check_cpu_timeout();
						check_depth_abort();
					}
					done_with_current_node();
				}
			}

			// Adds a new state into an existing matching node
			void schedule(const Node& n, const Node_ref match, const State& s, const Job<Time> &j, const Time est, const Time lst)
			{
				Interval<Time> finish_range = next_finish_times(n, s, j, est, lst);
				State& st = new_state(finish_range, s, j, match->get_scheduled_jobs());

				if (match->merge_states(st))
				{
					delete& st;
					num_states--;
				}
				else
				{
					if(use_supernodes)
					{
						DM("State not merged but added to the node");
						match->add_state(&st);
					}
					else
						schedule(n, s, j, est, lst);
				}
				process_new_edge(n, *match, j, finish_range);
			}

			// Creates a new node for a new state
			Node_ref schedule(const Node& n, const State& s, const Job<Time> &j, const Time est, const Time lst)
			{
				Interval<Time> finish_range = next_finish_times(n, s, j, est, lst);//requires est and lst
				DM("Creating a new node"<<std::endl);
				Node& next =
					new_node(finish_range, n, s, j, j.get_job_index(),
							  earliest_possible_job_release(n, j),
							  earliest_certain_source_job_release(n,j));
				//DM("      -----> N" << (nodes.end() - nodes.begin()) << " " <<(todo[todo_idx].front() - nodes.begin() + 1)
				//   << std::endl);
				State& next_state = new_state(finish_range, s, j, next.get_scheduled_jobs());
				next.add_state(&next_state);
				process_new_edge(n, next, j, finish_range);
				Node_ref n_ref = &(*(--nodes.end()));
				return n_ref;
			}

			void explore()
			{
				// #NS# The main workflow of building the graph starts from here
				make_initial_node();

				while (not_done() && !aborted) 
				{
					const Node& n = next_node();

					DM("\n==================================================="
					   << std::endl);
					/*DM("Looking at: N"
					   << (todo[todo_idx].front() - nodes.begin() + 1)
					   << " " << n << std::endl);*/

					// Obtain all the states in Node n
					const auto *n_states = n.get_states();

					// #NS# All statemenst below are to be calculated once per node and they are done that way 
					// The latest of all the lfts of of each state in node n
					auto lft_max = n.get_latest_core_availability();
					// among the incomplete jobs of node n, the time at which the earliest certain release occurs.
					auto nxt_ready_job = n.next_certain_job_ready_time();

					auto upbnd_twc = std::max(lft_max, nxt_ready_job);
					DM("=> upper-bound on twc at node level = "<< upbnd_twc << std::endl);

					// Keep track of the number of states that have led to new states. This is needed to detect
					// a deadend wherein a state cannot be expanded. The intial value is set to 0 and as we find states that 
					// can be expanded, we incremet this variable
					int num_states_expanded = 0;

					// Go through all the jobs that are not dispatched yet and are released at or before the latest time 
					// at which any job will certainly be scheduled in any state in the node n. 
					for (auto it = jobs_by_earliest_arrival.lower_bound(n.earliest_job_release()); 
							it != jobs_by_earliest_arrival.end() && it->first <= upbnd_twc;
							it++)
					{
						const Job<Time>& j = *(it->second);
						// j has been scheduled already, it will not be scheduled again
						if (!incomplete(n, j))
						{
							continue;
						}

						// j has predecessors yet to be scheduled, then it cannot be the next job scheduled 
						// (incomplete predecessors must be scheduled first)
						if (!ready(n, j)) {
							continue;
						}

						// nxt_node_created keeps track of whether a node in which we may add any state resulting from dispatching j already exists.
						// All states in node 'n' for which the job 'j' is eligible will be added to that same node. 
						bool nxt_node_created = false;
						// If a node was already created, we keep a reference to it
						Node_ref next_node;

						// this part of t_high is computed once per job, and is common to all the states explored
						Time t_high_wos = next_certain_higher_priority_source_job_release(n, j, upbnd_twc + 1);

						// if there is a higher priority job that is certainly ready before job j is released at the earliest, 
						// then j will never be the next job dispached by the scheduler
						if (t_high_wos <= j.earliest_arrival())
							continue;

						DM("\n---\nChecking for states to expand to for job " << j.get_id() << std::endl);
						for (State* s : *n_states)
						{
							// calculate the earliest time job j may start executing in state s
							Time est = next_earliest_start_time(n, *s, j);

							// Calculate lst, i.e., the latest time at which job j must start in state s to be 
							// the next job dispatched by the scheduler, using t_wc and t_high
							Time t_wc = std::max(s->latest_finish_time(), next_certain_job_ready_time(n, *s));
							Time t_high_ws = next_certain_higher_priority_successor_job_ready_time(n, *s, j, t_wc + 1);
							Time t_high = std::min(t_high_ws, t_high_wos);
							Time lst = std::min(t_wc, t_high - Time_model::constants<Time>::epsilon());

							// check if j may be scheduled next in system state s (i.e., the earliest time it may start at
							// is no later than the latest time it must start at)
							if (est <= lst)
							{
								if (nxt_node_created == false)
								{
									// Find the key of the next state from state s connected by an edge with job j
									// Find all states that have the same key
									// #NS# I wrote three different schedule functions to handle the three different cases
									// but this can be made better I guess cause they all do very similar things
									// The three similar functions are schedule_merge_node, schedule_merge_edge and schedule_new
									auto k = n.next_key(j);
									auto r = nodes_by_key.equal_range(k);

									if (r.first != r.second) 
									{
										Job_set sched_jobs{ n.get_scheduled_jobs(), j.get_job_index() };
										for (auto it = r.first; it != r.second; it++)
										{
											Node& found = *it->second;
 
											if (found.get_scheduled_jobs() != sched_jobs)
												continue;

											// If we have reached here, it means that we have found an existing node with the same 
											// set of scheduled jobs than the new state resuting from scheduling job j in system state s.
											// Thus, our new state can be added to that existing node.
											next_node = it->second;
											schedule(n, next_node, *s, j, est, lst);
											nxt_node_created = true;
											break;
										}									
									}

									// if there is no already existing nodes in which the new state can be added, then a new node is created
									if (nxt_node_created == false)
									{
										next_node = schedule(n, *s, j, est, lst);
										nxt_node_created = true;
									}
								}
								else
								{
									schedule(n, next_node, *s, j, est, lst);
								}

								// If the state is not a deadend, then we already know that the state has been expanded
								// so there is no need to increase the num_states_expanded variable. But if the state
								// is found to be a deadend, then this state can be set to not a deadend because the 
								// state has just been expanded.
								if (s->is_deadend())
								{
									s->not_deadend();
									num_states_expanded++;
								}
							}
						}
					}

					
					DM("---\nDone iterating over all jobs." << std::endl);

					// check for a dead end, check that all the states in the node has been expanded
					if (num_states_expanded != n.states_size() && n.get_scheduled_jobs().size() != jobs.size()) {
						// out of options and we didn't schedule all jobs
						observed_deadline_miss = true;
						if (early_exit)
							aborted = true;
						DM(":: Didn't find any possible successors." << std::endl);
					}

					check_cpu_timeout();
					check_depth_abort();
					done_with_current_node();
				}
			}

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			friend std::ostream& operator<< (std::ostream& out,
											 const State_space<Time>& space)
			{
					std::map<const Schedule_node<Time>*, unsigned int> node_id;
					unsigned int i = 1;
					out << "digraph {" << std::endl;
					for (const Schedule_node<Time>& n : space.get_nodes()) {
						node_id[&n] = i++;
						out << "\tN" << node_id[&n]
							<< "[label=\"N" << node_id[&n] << ": {";
						const auto *n_states = n.get_states();

						for(State *s: *n_states)
						{
							out << "["
								<< s->earliest_finish_time()
								<< ", "
								<< s->latest_finish_time()
								<< "]\\n";
						}
						out << "}"
							<< "\\nER=";
						if (n.earliest_job_release() ==
							Time_model::constants<Time>::infinity()) {
							out << "N/A";
						} else {
							out << n.earliest_job_release();
						}
						out << "\"];"
							<< std::endl;
					}
					for (auto e : space.get_edges ()) {
						out << "\tN" << node_id[e.source]
							<< " -> "
							<< "N" << node_id[e.target]
							<< "[label=\""
							<< "T" << e.scheduled->get_task_id()
							<< " J" << e.scheduled->get_job_id()
							<< "\\nDL=" << e.scheduled->get_deadline()
							<< "\\nES=" << e.earliest_start_time()
							<< "\\nLS=" << e.latest_start_time()
							<< "\\nEF=" << e.earliest_finish_time()
							<< "\\nLF=" << e.latest_finish_time()
							<< "\"";
						if (e.deadline_miss_possible()) {
							out << ",color=Red,fontcolor=Red";
						}
						out << ",fontsize=8" << "]"
							<< ";"
							<< std::endl;
						if (e.deadline_miss_possible()) {
							out << "N" << node_id[e.target]
								<< "[color=Red];"
								<< std::endl;
						}
					}
					out << "}" << std::endl;
				return out;
			}
#endif
		};

	}
}

// #include "uni_iip/iip.hpp"

namespace std
{
	template<class Time> struct hash<NP::Uniproc::Schedule_node<Time>>
	{
		std::size_t operator()(NP::Uniproc::Schedule_node<Time> const& n) const
		{
			return n.get_key();
		}
	};
}


#endif
