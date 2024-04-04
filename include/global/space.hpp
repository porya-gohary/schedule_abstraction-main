#ifndef GLOBAL_SPACE_H
#define GLOBAL_SPACE_H

#include <unordered_map>
#include <map>
#include <vector>
#include <deque>
#include <forward_list>
#include <algorithm>

#include <iostream>
#include <ostream>
#include <cassert>

#include "config.h"

#ifdef CONFIG_PARALLEL
#include "tbb/concurrent_hash_map.h"
#include "tbb/enumerable_thread_specific.h"
#include "tbb/parallel_for.h"
#endif

#include "problem.hpp"
#include "clock.hpp"

#include "global/state.hpp"

#define NOSUSP 0
#define GENERAL_SUSP 1
#define PATHWISE_SUSP 2

namespace NP {

	namespace Global {

		template<class Time> class State_space
		{
			public:

			typedef Scheduling_problem<Time> Problem;
			typedef typename Scheduling_problem<Time>::Workload Workload;
			typedef typename Scheduling_problem<Time>::Suspending_Tasks Suspending_Tasks;
			typedef Schedule_state<Time> State;
		  typedef typename std::vector<Interval<Time>> CoreAvailability;

		  typedef Schedule_node<Time> Node;

			static State_space explore(
					const Problem& prob,
					const Analysis_options& opts)
			{
				// doesn't yet support exploration after deadline miss
				assert(opts.early_exit);

				auto s = State_space(prob.jobs, prob.dag, prob.sts, prob.num_processors, opts.timeout,
									 opts.max_depth, opts.num_buckets, opts.use_self_suspensions);
				s.be_naive = opts.be_naive;
				s.cpu_time.start();
				s.explore(); // ISSUE
				s.cpu_time.stop();
				return s;

			}

			// convenience interface for tests
			static State_space explore_naively(
				const Workload& jobs,
				unsigned int num_cpus)
			{
				Problem p{jobs, num_cpus};
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static State_space explore(
				const Workload& jobs,
				unsigned int num_cpus)
			{
				Problem p{jobs, num_cpus};
				Analysis_options o;
				return explore(p, o);
			}

			Interval<Time> get_finish_times(const Job<Time>& j) const
			{
				auto rbounds = rta.find(j.get_job_index());
				if (rbounds == rta.end()) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rbounds->second;
				}
			}

			Interval<Time> get_finish_times(Job_index j) const
			{
				auto rbounds = rta.find(j);
				if (rbounds == rta.end()) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rbounds->second;
				}
			}

			bool is_schedulable() const
			{
				return !aborted;
			}

			bool was_timed_out() const
			{
				return timed_out;
			}

			//currently unused, only required to compile nptest.cpp correctly
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

		  typedef std::deque<Node> Nodes;
			typedef std::deque<State> States;

#ifdef CONFIG_PARALLEL
			typedef tbb::enumerable_thread_specific< States > Split_states;
			typedef std::deque<Split_states> States_storage;
			typedef tbb::enumerable_thread_specific< Nodes > Split_nodes;
			typedef std::deque<Split_nodes> Nodes_storage;
#else
			typedef std::deque< States > States_storage;
			typedef std::deque< Nodes > Nodes_storage;
#endif

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH

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

			const States_storage& get_states() const
			{
				return states_storage;
			}

		  const Nodes_storage& get_nodes() const
		  {
		    return nodes_storage;
		  }
		  

#endif
			private:

		  typedef Node* Node_ref;
		  //typedef typename std::deque<Node>::iterator Node_ref;
		  typedef typename std::forward_list<Node_ref> Node_refs;
		        typedef State* State_ref;  // RV:  upstream:master
		  //typedef typename std::deque<State>::iterator State_ref;
			typedef typename std::forward_list<State_ref> State_refs;

#ifdef CONFIG_PARALLEL
			typedef tbb::concurrent_hash_map<hash_value_t, State_refs> States_map;
			typedef typename States_map::accessor States_map_accessor;
#else
			typedef std::unordered_map<hash_value_t, State_refs> States_map;
#endif
		  typedef std::unordered_map<hash_value_t, Node_refs> Nodes_map;
		  
		  
			typedef const Job<Time>* Job_ref;
			typedef std::multimap<Time, Job_ref> By_time_map;

		  // typedef std::deque<State_ref> Todo_queue;
		  typedef std::deque<Node_ref> Todo_queue;

			typedef Interval_lookup_table<Time, Job<Time>, Job<Time>::scheduling_window> Jobs_lut;

			typedef std::unordered_map<Job_index, Interval<Time> > Response_times;

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			std::deque<Edge> edges;
#endif

			Response_times rta;

#ifdef CONFIG_PARALLEL
			tbb::enumerable_thread_specific<Response_times> partial_rta;
#endif

			bool aborted;
			bool timed_out;

			const unsigned int max_depth;

			bool be_naive;

			const Workload& jobs;

			// not touched after initialization
			Jobs_lut _jobs_by_win;
			By_time_map _jobs_by_latest_arrival_with_susp;
			By_time_map _jobs_by_latest_arrival_without_susp;
			By_time_map _jobs_by_earliest_arrival;
			By_time_map _jobs_by_deadline;
			std::vector<Job_precedence_set> _predecessors;

			// use these const references to ensure read-only access
			const Jobs_lut& jobs_by_win;
			const By_time_map& jobs_by_latest_arrival_with_susp;
			const By_time_map& jobs_by_latest_arrival_without_susp;
			const By_time_map& jobs_by_earliest_arrival;
			const By_time_map& jobs_by_deadline;
			const std::vector<Job_precedence_set>& predecessors;

			// In order to store the components of self-suspending tasks, we create a vector of a vector of references
			// to_suspending_tasks: for each successor, a list of predecessors are present
			// from_suspending_tasks: for each predecessor, a list of successors are present
		        // RV: how is this used? Similar to definitions above, is a const reference useful?
		        //     It might be useful to store the Job_index in the suspending_task_list.
			typedef std::vector<const Suspending_Task<Time>*> suspending_task_to_list;
		  typedef std::vector<Job_index> suspending_task_from_list;
		  
			std::vector<suspending_task_to_list> _to_suspending_tasks;
			std::vector<suspending_task_from_list> _from_suspending_tasks;

			const std::vector<suspending_task_to_list>& to_suspending_tasks;
			const std::vector<suspending_task_from_list>&  from_suspending_tasks;

			States_storage states_storage;
		  Nodes_storage nodes_storage;

		  //Nodes nodes;
		  //States states;
		  Nodes_map nodes_by_key;
			States_map states_by_key;
			// updated only by main thread
			unsigned long num_nodes, num_states, width;
			unsigned long current_job_count;
			unsigned long num_edges;

		  static const std::size_t num_todo_queues = 3;
		  Todo_queue todo[num_todo_queues];
		  int todo_idx;

#ifdef CONFIG_PARALLEL
			tbb::enumerable_thread_specific<unsigned long> edge_counter;
#endif
			Processor_clock cpu_time;
			const double timeout;
			unsigned int wants_self_suspensions;

			const unsigned int num_cpus;

			State_space(const Workload& jobs,
						const Precedence_constraints &dag_edges,
						const Suspending_Tasks &susps,
						unsigned int num_cpus,
						double max_cpu_time = 0,
						unsigned int max_depth = 0,
						std::size_t num_buckets = 1000,
						unsigned int use_self_suspensions = NOSUSP)
			: _jobs_by_win(Interval<Time>{0, max_deadline(jobs)},
						   max_deadline(jobs) / num_buckets)
			, jobs(jobs)
			, aborted(false)
			, timed_out(false)
			, be_naive(false)
			, timeout(max_cpu_time)
			, max_depth(max_depth)
			, num_nodes(0)
			, num_states(0)
			, num_edges(0)
			, width(0)
			, current_job_count(0)
			, num_cpus(num_cpus)
			, jobs_by_latest_arrival_with_susp(_jobs_by_latest_arrival_with_susp)
			, jobs_by_latest_arrival_without_susp(_jobs_by_latest_arrival_without_susp)
			, jobs_by_earliest_arrival(_jobs_by_earliest_arrival)
			, jobs_by_deadline(_jobs_by_deadline)
			, jobs_by_win(_jobs_by_win)
			, _predecessors(jobs.size())
			, predecessors(_predecessors)
			, _to_suspending_tasks(jobs.size())
			, _from_suspending_tasks(jobs.size())
			, to_suspending_tasks(_to_suspending_tasks)
			, from_suspending_tasks(_from_suspending_tasks)
			, wants_self_suspensions(use_self_suspensions)
			{
				for (auto e : dag_edges) {
					const Job<Time>& from = lookup<Time>(jobs, e.first);
					const Job<Time>& to   = lookup<Time>(jobs, e.second);
					_predecessors[to.get_job_index()].push_back(from.get_job_index());
				}

				// RV: initialization of to_suspending_tasks and from_suspending_tasks.
				//     Is st.get_toID() equal to index_of(to) ?
				//     Is _predecessors[] related to to_suspending_tasks[] ? 
				for (const Suspending_Task<Time>& st : susps) {
					_to_suspending_tasks[st.get_toIndex()].push_back(&st);
					_from_suspending_tasks[st.get_fromIndex()].push_back(st.get_toIndex());
				}

				for (const Job<Time>& j : jobs) {
				  if (_to_suspending_tasks[j.get_job_index()].size()>0) {
					_jobs_by_latest_arrival_with_susp.insert({j.latest_arrival(), &j});
				  } else {
					_jobs_by_latest_arrival_without_susp.insert({j.latest_arrival(), &j});
				  }
					_jobs_by_earliest_arrival.insert({j.earliest_arrival(), &j});
					_jobs_by_deadline.insert({j.get_deadline(), &j});
					_jobs_by_win.insert(j);
				}

			}

			private:

			void count_edge()
			{
#ifdef CONFIG_PARALLEL
				edge_counter.local()++;
#else
				num_edges++;
#endif
			}

			static Time max_deadline(const Workload &jobs)
			{
				Time dl = 0;
				for (const auto& j : jobs)
					dl = std::max(dl, j.get_deadline());
				return dl;
			}

			void update_finish_times(Response_times& r, const Job_index id,
									 Interval<Time> range)
			{
				auto rbounds = r.find(id);
				if (rbounds == r.end()) {
					r.emplace(id, range);
				} else {
					rbounds->second |= range;
				}
				DM("RTA " << id << ": " << r.find(id)->second << std::endl);
			}

			void update_finish_times(
				Response_times& r, const Job<Time>& j, Interval<Time> range)
			{
				update_finish_times(r, j.get_job_index(), range);
				if (j.exceeds_deadline(range.upto()))
					aborted = true;
			}

			void update_finish_times(const Job<Time>& j, Interval<Time> range)
			{
				Response_times& r =
#ifdef CONFIG_PARALLEL
					partial_rta.local();
#else
					rta;
#endif
				update_finish_times(r, j, range);
			}


			std::size_t index_of(const Job<Time>& j) const
			{
				return (std::size_t) (&j - &(jobs[0]));
			}
			const Job<Time>& reverse_index_of(std::size_t index) const
			{
				return jobs[index];
			}

			const Job_precedence_set& predecessors_of(const Job<Time>& j) const
			{
			  return predecessors[j.get_job_index()];
			}

		  // RV: it seems that jobs can fall through the selection algorithm.
		  //     This function seems state specific due to core_availability.
		  void check_for_deadline_misses(const Node& old_n, const Node& new_n)
			{
			  auto check_from = old_n.get_first_state()->core_availability().min();
			  auto earliest   = new_n.get_last_state()->core_availability().min();

			  // check if we skipped any jobs that are now guaranteed
			  // to miss their deadline
			  for (auto it = jobs_by_deadline.lower_bound(check_from);
			       it != jobs_by_deadline.end(); it++) {
			    const Job<Time>& j = *(it->second);
			    if (j.get_deadline() < earliest) {
			      if (unfinished(new_n, j)) {
				DM("deadline miss: " << new_n << " -> " << j << std::endl);
							// This job is still incomplete but has no chance
							// of being scheduled before its deadline anymore.
							// Abort.
							aborted = true;
							// create a dummy node for explanation purposes
							auto frange = new_n.get_last_state()->core_availability() + j.get_cost();
							Node& next =
							  new_node(new_n, j, j.get_job_index(), 0);
							const CoreAvailability empty_cav={};
							State& next_s = new_state(*new_n.get_last_state(), j.get_job_index(), predecessors_of(j), frange, frange, empty_cav);  // ISSUE
							next.add_state(&next_s);

							// update response times
							update_finish_times(j, frange);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
							edges.emplace_back(&j, &new_n, &next, frange);
#endif
							count_edge();
							break;
						}
					} else
						// deadlines now after the next earliest finish time
						break;
				}
			}

		        // RV: nested for loops. Unclear whether job_info.first == index_of(lookup(jobs, job_check_id)) .
		        //     The code seems to use references to a Job and index_of a Job.
		        //     Is "b=reverse_index_of(a); c=b.get_id(); d=lookup(jobs,c); e=index_of(d);" equal to "e=a;"?
		        //     Can this function be called per node?
		  void remove_jobs_with_no_successors(State& s, const Node& n)
			{
				std::vector<Job<Time>> jobsToRemove;
				// remove jobs with no successors
				for (auto job_info : s.get_pathwise_jobs()) {
					Job_index job_check = job_info.first;
					bool successor_pending = false;
					//unsigned long jif=job_info.first;
					//unsigned long iol=index_of(lookup<Time>(jobs, job_check_id));
					//std::cerr << "JIF: " << jif << ", iol: " << iol << "\n";
					if (from_suspending_tasks[job_check].size() > 0) {
						for (auto to_jobs : from_suspending_tasks[job_check]) {
							if(n.job_incomplete(to_jobs)) {
								successor_pending = true;
								break;
							}
						}
					}

					if (!successor_pending) {
						jobsToRemove.push_back(job_check);
					}
				}

				for (auto job : jobsToRemove) {
				  s.del_pred(job.get_job_index());
				}
			}

			void make_initial_node(unsigned num_cores)
			{
				// construct initial state
				nodes_storage.emplace_back();
				states_storage.emplace_back();
				Node& n=new_node(num_cores);
				State& s=new_state(num_cores);
				n.add_state(&s);
			}

			States& states()
			{
#ifdef CONFIG_PARALLEL
				return states_storage.back().local();
#else
				return states_storage.back();
#endif
			}

			Nodes& nodes()
			{
#ifdef CONFIG_PARALLEL
				return nodes_storage.back().local();
#else
				return nodes_storage.back();
#endif
			}

			template <typename... Args>
			State_ref alloc_state(Args&&... args)
			{
			  states().emplace_back(std::forward<Args>(args)...); // ISSUE
				State_ref s = &(*(--states().end()));   // RV:  gn:master
				//State_ref s = --states().end();

				// RV: requires Node& parameter.
				//remove_jobs_with_no_successors(*s);

				// make sure we didn't screw up...
				//auto njobs = n->number_of_scheduled_jobs();
				//assert (
				//	(!njobs && num_states == 0) // initial state
				//	|| (njobs == currecomiplent_job_count + 1) // normal State
				//	|| (njobs == current_job_count + 2 && aborted) // deadline miss
				//);

				return s;
			}

			void dealloc_state(State_ref s)
			{
			        assert(&(*(--states().end())) == s);   // RV:  gn:master
				//assert(--states().end() == s);
				states().pop_back();
			}

			template <typename... Args>
			Node_ref alloc_node(Args&&... args)
			{
				nodes().emplace_back(std::forward<Args>(args)...);
				Node_ref n = &(*(--nodes().end()));   // RV:  gn:master
				//Node_ref n = --nodes().end();

				//remove_jobs_with_no_successors(*s);

				// make sure we didn't screw up...
				auto njobs = n->number_of_scheduled_jobs();
				assert (
					(!njobs && num_states == 0) // initial state
					|| (njobs == current_job_count + 1) // normal State
					|| (njobs == current_job_count + 2 && aborted) // deadline miss
				);

				return n;
			}

			void dealloc_node(Node_ref n)
			{
			        assert(&(*(--nodes().end())) == n);   // RV:  gn:master
				//assert(--nodes().end() == n);
				nodes().pop_back();
			}

			template <typename... Args>
			State& new_state(Args&&... args)
			{
			  return *alloc_state(std::forward<Args>(args)...); // ISSUE
			}

			template <typename... Args>
			State& new_or_merged_state(Args&&... args)
			{
			  State_ref s_ref = alloc_state(std::forward<Args>(args)...);

				// try to merge the new state into an existing state
				State_ref s = merge_or_cache_state(s_ref);
				if (s != s_ref) {
					// great, we merged!
					// clean up the just-created state that we no longer need
					dealloc_state(s_ref);
				}
				return *s;
			}

			template <typename... Args>
			Node& new_node(Args&&... args)
			{
				Node_ref n = alloc_node(std::forward<Args>(args)...);
				//State& st = new_state(num_cpus);
				//n->add_state(&st);

				auto njobs = n->get_scheduled_jobs().size();
				auto idx = njobs % num_todo_queues;
				todo[idx].push_back(n);
				//RV:  add node to nodes_by_key map.
				cache_node(n);
				num_nodes++;
				width = std::max(width, (unsigned long) todo[idx].size()-1);
				return *n;
			}

			template <typename... Args>
			Node& new_or_merged_node(Args&&... args)
			{
				Node_ref n_ref = alloc_node(std::forward<Args>(args)...);

				// try to merge the new state into an existing state
				Node_ref n = merge_or_cache_node(n_ref);
				if (n != n_ref) {
					// great, we merged!
					// clean up the just-created state that we no longer need
					dealloc_node(n_ref);
				}
				return *n;
			}

#ifdef CONFIG_PARALLEL
#warning  "Parallel code is not updated for supernodes."
			// make state available for fast lookup
			void insert_cache_state(States_map_accessor &acc, State_ref s)
			{
				assert(!acc.empty());

				State_refs& list = acc->second;
				list.push_front(s);
			}

			// returns true if state was merged
			State_ref merge_or_cache_state(State_ref s)
			{
				States_map_accessor acc;

				while (true) {
					// check if key exists
					if (states_by_key.find(acc, s->get_key())) {
						for (State_ref other : acc->second)
							if (other->try_to_merge(*s))
								return other;
						// If we reach here, we failed to merge, so go ahead
						// and insert it.
						insert_cache_state(acc, s);
						return s;
					// otherwise, key doesn't exist yet, let's try to create it
					} else if (states_by_key.insert(acc, s->get_key())) {
						// We created the list, so go ahead and insert our state.
						insert_cache_state(acc, s);
						return s;
					}
					// if we raced with concurrent creation, try again
				}
			}

			// make node available for fast lookup
			void insert_cache_node(Nodes_map_accessor &acc, Node_ref n)
			{
				assert(!acc.empty());

				Node_refs& list = acc->second;
				list.push_front(n);
			}

			// returns true if node was merged
			Node_ref merge_or_cache_node(Node_ref n)
			{
				Nodes_map_accessor acc;

				while (true) {
					// check if key exists
					if (nodes_by_key.find(acc, n->get_key())) {
						for (Node_ref other : acc->second)
							if (other->try_to_merge(*n))
								return other;
						// If we reach here, we failed to merge, so go ahead
						// and insert it.
						insert_cache_node(acc, n);
						return n;
					// otherwise, key doesn't exist yet, let's try to create it
					} else if (nodes_by_key.insert(acc, n->get_key())) {
						// We created the list, so go ahead and insert our node.
						insert_cache_node(acc, n);
						return n;
					}
					// if we raced with concurrent creation, try again
				}
			}

#else

			void cache_state(State_ref s)
			{
				// create a new list if needed, or lookup if already existing
				auto res = states_by_key.emplace(
					std::make_pair(0x45454, State_refs()));

				auto pair_it = res.first;
				State_refs& list = pair_it->second;

				list.push_front(s);
			}


			State_ref merge_or_cache_state(State_ref s_ref)
			{
				State& s = *s_ref;

				const auto pair_it = states_by_key.find(0x45454);

				// cannot merge if key doesn't exist
				if (pair_it != states_by_key.end())
					for (State_ref other : pair_it->second)
						if (other->try_to_merge(*s_ref))
							return other;
				// if we reach here, we failed to merge
				cache_state(s_ref);
				return s_ref;
			}
		  
			void cache_node(Node_ref n)
			{
				// create a new list if needed, or lookup if already existing
				auto res = nodes_by_key.emplace(
					std::make_pair(n->get_key(), Node_refs()));

				auto pair_it = res.first;
				Node_refs& list = pair_it->second;

				list.push_front(n);
			}


			Node_ref merge_or_cache_node(Node_ref n_ref)
			{
				Node& n = *n_ref;

				const auto pair_it = nodes_by_key.find(n.get_key());

				// cannot merge if key doesn't exist
				if (pair_it != nodes_by_key.end())
					for (Node_ref other : pair_it->second)
						if (other->try_to_merge(*n_ref))
							return other;
				// if we reach here, we failed to merge
				cache_node(n_ref);
				return n_ref;
			}
#endif

			void check_cpu_timeout()
			{
				if (timeout && get_cpu_time() > timeout) {
					aborted = true;
					timed_out = true;
				}
			}

			void check_depth_abort()
			{
				if (max_depth && current_job_count > max_depth)
					aborted = true;
			}

			bool unfinished(const Node& n, const Job<Time>& j) const
			{
			  return n.job_incomplete(j.get_job_index());
			}

		  // when a job has no self-suspension, state is not required.
		  bool ready(const Node& n, const Job<Time>& j) const
			{
			  return unfinished(n, j) && n.job_ready(predecessors_of(j));
			}

		  // when a job has self-suspension, susp_ready() should be called.
		  // however, it doesn't check the job_finish_times.
		  // susp_ready_at() is not used.
		  bool ready(const Node& n, const State& s, const Job<Time>& j) const
			{
			  return unfinished(n, j) && n.job_ready(predecessors_of(j)) && susp_ready(n,j);
			}

		        // RV: if lookup() is time consuming and jobs is fixed, store the result in e.
		        //     susp_ready() is used in assert() calls below, so performance might be relevant.
		  bool susp_ready(const Node&n, const Job<Time>& j) const
			{
			  const suspending_task_to_list fsusps = to_suspending_tasks[j.get_job_index()];

				for (auto e : fsusps)
				{
					if (n.job_incomplete(e->get_fromIndex()))
					{
						return false;
					}
				}
				return true;
			}

		  // RV:  Not passing node:  susp_ready(n,s,j) is assumed to be true.
		  Time get_seft(const State& s, const Job<Time>& j) const
			{
			  const suspending_task_to_list fsusps = to_suspending_tasks[j.get_job_index()];
			  Time seft = 0;  // RV:  what if the list is empty?
				// assert(susp_ready(n, s, j)); //RV:  susp_ready requires Node.

				for (auto e : fsusps)
				{
					Interval<Time> rbounds = get_finish_times(e->get_fromIndex());
					if (wants_self_suspensions == PATHWISE_SUSP)
						rbounds = s.get_pathwisejob_ft(e->get_fromIndex());

					if (seft <= (rbounds.from() + e->get_minsus()))
					{
						seft = rbounds.from() + e->get_minsus();
					}
				}

				return seft;
			}

		  Time get_slft(const State& s, const Job<Time>& j) const
			{
			  const suspending_task_to_list fsusps = to_suspending_tasks[j.get_job_index()];
				Time slft = 0;
				// assert(susp_ready(n,s,j)); // RV: susp_ready() requires Node 

				for(auto e: fsusps)
				{
					Interval<Time> rbounds = get_finish_times(e->get_fromIndex());
					if(wants_self_suspensions == PATHWISE_SUSP)
					  rbounds = s.get_pathwisejob_ft(e->get_fromIndex());

					if(slft <= (rbounds.until() + e->get_maxsus()))
					{
						slft = rbounds.until() + e->get_maxsus();
					}
				}

				return slft;
			}

		  // RV: function is never used. Assume susp_ready() as precondition.
		  bool susp_ready_at(const Node& n, const State& s, const Job<Time>& j, const Time at) const
			{
				const suspending_task_to_list fsusps = to_suspending_tasks[j.get_job_index()];

				for(auto e: fsusps)
				{
				  if(get_slft(s,j) > at)
				    return false;
				}
				return true;
			}

			bool all_jobs_scheduled(const Node& n) const
			{
				return n.number_of_scheduled_jobs() == jobs.size();
			}

			// assumes j is ready, including susp_ready()
		  Interval<Time> ready_times(const State& s, const Job<Time>& j) const
			{
				Interval<Time> r = j.arrival_window();
				for (auto pred : predecessors_of(j)) {
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				// RV: oldFrom, newFrom, oldUntil and newUntil are unused.
				//     are get_seft() and get_slft() results reused?
				Time seft = get_seft(s, j);
				Time slft = get_slft(s, j);
				if (r.from() < seft)
				{
				  //Time oldFrom = r.from();
					r.lower_bound(seft);
					//Time newFrom = r.from();
				}
				if (r.until() < slft)
				{
				  //Time oldUntil = r.until();
					r.extend_to(slft);
					//Time newUntil = r.until();
				}
				return r;
			}

			// assumes j is ready
			Interval<Time> ready_times(
						   const State& s, const Job<Time>& j,
				const Job_precedence_set& disregard) const
			{
				Interval<Time> r = j.arrival_window();
				for (auto pred : predecessors_of(j)) {
					// skip if part of disregard
					if (contains(disregard, pred))
						continue;
					Interval<Time> ft{0, 0};
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				// RV: similar to previous ready_times()
				Time seft = get_seft(s, j);
				Time slft = get_slft(s, j);
				if(r.from() < seft)
					r.lower_bound(seft);
				if(r.until() < slft)
					r.extend_to(slft);
				return r;
			}

		  Time latest_ready_time(const State& s, const Job<Time>& j) const
			{
			  return ready_times(s, j).max();
			}

		  // When a job doesn't depend on suspension, is node information sufficient?
		  Time latest_ready_time(const Node& n, const Job<Time>& j) const
			{
			  return j.arrival_window().max();
			  // return ready_times(s, j).max();
			}

		  Time earliest_ready_time(const State& s, const Job<Time>& j) const
			{
			  return ready_times(s, j).min();
			}

		  Time earliest_ready_time(const Node& s, const Job<Time>& j) const
			{
			  return j.arrival_window().min();
			  //return ready_times(s, j).min();
			}

			Time latest_ready_time(
				const State& s, Time earliest_ref_ready,
				const Job<Time>& j_hp, const Job<Time>& j_ref) const
			{
			  auto rt = ready_times(s, j_hp, predecessors_of(j_ref));
			  return std::max(rt.max(), earliest_ref_ready);
			}

		  //RV: for the no_suspension functions, which don't have state.
		  Time latest_ready_time(
				const Node& n, Time earliest_ref_ready,
				const Job<Time>& j_hp, const Job<Time>& j_ref) const
			{
			  Time t = j_hp.arrival_window().max();
			  return std::max(t, earliest_ref_ready);
			}

		  //RV: called next_certain_higher_priority_job_release*() in uni/space.hpp
		  //    what is the difference between ready and release?
			// Find next time by which any job is certainly released.
			// Note that this time may be in the past.  (RV: how?)
			Time next_higher_prio_job_ready_no_suspension(
							const Node& n,
				const Job<Time> &reference_job,
				const Time t_earliest) const
			{
				auto ready_min = earliest_ready_time(n, reference_job);
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
					if (ready(n, j)
						&& j.higher_priority_than(reference_job)) {
					  //RV: latest_ready_time might depend on state due to predecessors.
						when = std::min(when,
								latest_ready_time(n, ready_min, j, reference_job));
					}

				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival_without_susp
							   .lower_bound(t_earliest);
					 it != jobs_by_latest_arrival_without_susp.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(n, j)
						&& j.higher_priority_than(reference_job)) {
						// does it beat what we've already seen?
					  //RV: does latest_ready_time depend on state ?
						when = std::min(when,
							latest_ready_time(n, ready_min, j, reference_job));
					}
				}

				return when;
			}

			// Find next time by which any job is certainly released.
			// Note that this time may be in the past.
			Time next_higher_prio_job_ready_with_suspension(
							const Node& n,
				const State& s,
				const Job<Time> &reference_job,
							const Time t_earliest,
							Time max=Time_model::constants<Time>::infinity()) const
			{
			  auto ready_min = earliest_ready_time(s, reference_job);
				Time when = max;
				// check everything that overlaps with t_earliest
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
				  if (ready(n, s, j)
						&& j.higher_priority_than(reference_job)) {
						when = std::min(when,
								latest_ready_time(s, ready_min, j, reference_job));
					}

				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival_with_susp
							   .lower_bound(t_earliest);
					 it != jobs_by_latest_arrival_with_susp.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(n, j)
						&& j.higher_priority_than(reference_job)) {
						// does it beat what we've already seen?
						when = std::min(when,
								latest_ready_time(s, ready_min, j, reference_job));
					}
				}

				return when;
			}

		  Time next_higher_prio_job_ready(
						  const Node& n,
						  const State& s,
						  const Job<Time> &reference_job,
						  const Time t_earliest) const
		  {
		    Time t_ws = next_higher_prio_job_ready_no_suspension(n, reference_job, t_earliest);
		    Time t_wos = next_higher_prio_job_ready_with_suspension(n,s,reference_job, t_earliest, t_ws);
		    return t_wos;
		  }

			// Find next time by which any job is certainly released.
			// Note that this time may be in the past.
		  Time next_job_ready_no_suspension(const Node& n, const Time t_earliest) const
			{
				Time when = Time_model::constants<Time>::infinity();

				// check everything that overlaps with t_earliest
				// RV: should state be given to latest_ready_time?
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
					if (ready(n, j)){
						when = std::min(when, latest_ready_time(n, j));
					}


				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival_without_susp
							   .lower_bound(t_earliest);
					 it != jobs_by_latest_arrival_without_susp.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(n, j))
						// does it beat what we've already seen?
						when = std::min(when, latest_ready_time(n, j));
				}

				return when;
			}

		  Time next_job_ready_with_suspension(const Node& n, const State& s, const Time t_earliest, Time max=Time_model::constants<Time>::infinity()) const
			{
				Time when = max;

				// check everything that overlaps with t_earliest
				for (const Job<Time>& j : jobs_by_win.lookup(t_earliest))
				  if (ready(n, j)){
				    when = std::min(when, latest_ready_time(s, j));
					}


				// No point looking in the future when we've already
				// found one in the present.
				if (when <= t_earliest)
					return when;

				// Ok, let's look also in the future.
				for (auto it = jobs_by_latest_arrival_with_susp
							   .lower_bound(t_earliest);
					 it != jobs_by_latest_arrival_with_susp.end(); it++) {
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or blocked
					if (ready(n, j))
						// does it beat what we've already seen?
					  when = std::min(when, latest_ready_time(s, j));
				}

				return when;
			}

		  Time next_job_ready(const Node& n, const State& s, const Time t_earliest) const
		  {
		    Time t_ws = next_job_ready_no_suspension(n, t_earliest);
		    Time t_wos = next_job_ready_with_suspension(n,s,t_earliest, t_ws);
		    return t_wos;
		  }

			// assumes j is ready
			// NOTE: we don't use Interval<Time> here because the
			//       Interval c'tor sorts its arguments.
		  std::pair<Time, Time> start_times(
						    const State& s, const Job<Time>& j, Time t_wc, Time t_high) const
			{
			  auto rt = ready_times(s, j);  //RV: due to suspension, Node is needed. 
				auto at = s.core_availability();
				Time est = std::max(rt.min(), at.min());

				DM("rt: " << rt << std::endl
				<< "at: " << at << std::endl);

				// auto t_high = next_higher_prio_job_ready(n, s, j, at.min());
				Time lst    = std::min(t_wc,
					t_high - Time_model::constants<Time>::epsilon());

				DM("est: " << est << std::endl);
				DM("lst: " << lst << std::endl);

				return {est, lst};
			}

		  bool dispatch(const Node& n, const Job<Time>& j, Time t_wc_wos, Time t_high_wos)
			{
			  // check if this job has a feasible start-time interval
			  // Given a node and job, the new state would be added to the same
			  // next node. Some optimisation might be possible w.r.t. merging.
			  Node_ref next = nullptr;
			  const auto pair_it = nodes_by_key.find(n.next_key(j));
			  if (pair_it != nodes_by_key.end()) {
			    for (Node_ref other : pair_it->second) {
			      // Assume there are no key collisions for nodes.
			      // RV: If there are collisions, a node could adjust its key, such
			      //     that future nodes don't have collisions. 
			      next = other;
			      break;
			    }
			  }
			  const auto *n_states = n.get_states();
			  for (State *s : *n_states) {
			    Time t_wc_ws = next_job_ready_with_suspension(n, *s, t_wc_wos);
			    Time t_high_ws = next_higher_prio_job_ready_with_suspension(n,*s,j,t_high_wos);
			    
			    Time t_high = std::min(t_high_ws, t_high_wos);
			    Time t_wc = std::min(t_wc_ws, t_wc_wos);
			    
			    auto _st = start_times(*s, j, t_wc, t_high);  // RV: no node argument?
			    if (_st.first > _st.second)
			      return false; // nope, not schedulable, return.
			    
			    Interval<Time> st{_st};
			    
			    // yep, job j is a feasible successor in state s
			    
			    // compute range of possible finish times
			    Interval<Time> ftimes = st + j.get_cost();
			    
			    // update finish-time estimates
			    update_finish_times(j, ftimes);
			    
			    // expand the graph, merging if possible
			    // RV: compute loopup key
			    //     check whether a node exists
			    //     + merge with state in node or add extra state.
			    //     - create a new node.
			    // Compute core_avail for the state where j is started.
			    CoreAvailability cav;
			    s->next_core_avail(j.get_job_index(), predecessors_of(j), st, ftimes, cav);
			    // If there is no node yet, create one.
			    // If be_naive, a new node should be created time.
			    if (next == nullptr) {
			      // RV: earliest_pending_release?
			      next = &(new_node(n, j, index_of(j), 0));
			      // cache_node(next);  // done in new_node.
			      // num_nodes++;       // done in new_node.
			      // RV: Add an edge?
			    }
			    // next should always exist at this point, possibly without states
			    // Merge with existing states if they exist.
			    if (!next->merge_states(st, ftimes, *s, j, cav))
			      {
				// If merging didn't work, add a new state.
				State& new_s = new_state(*s, index_of(j), predecessors_of(j),
							 st, ftimes, cav);
				next->add_state(s);
			      }
			    
			    // make sure we didn't skip any jobs
			    // RV: is an easier detection possible?
			    //     can detection be done per node (and earliest core_available?)
			    if (be_naive) {
			      check_for_deadline_misses(n,*next); // ISSUE
			    }
			    
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			    // RV: should this be done for every state.
			    edges.emplace_back(&j, &s, next, ftimes);
#endif
			    count_edge();
			  }

			  if (next != nullptr)
			    check_for_deadline_misses(n, *next);

			  return true;
			}

			void explore(const Node& n)
			{
			  bool found_one = false;
			  
			  DM("----" << std::endl);
			  
			  // (0) define the window of interest
			  const State& s_first = *(n.get_first_state());
			  // earliest time a core is possibly available
			  auto t_min  = s_first.core_availability().min();
			  // latest time some unfinished job is certainly ready
			  auto t_job  = next_job_ready(n, s_first, t_min);
			  // latest time some core is certainly available
			  const State& s_last = *(n.get_last_state());
			  auto t_core = s_last.core_availability().max();
			  // latest time by which a work-conserving scheduler
			  // certainly schedules some job
			  auto t_wc   = std::max(t_core, t_job);

			  DM(s << std::endl);
			  DM("t_min: " << t_min << std::endl
			     << "t_job: " << t_job << std::endl
			     << "t_core: " << t_core << std::endl
			     << "t_wc: " << t_wc << std::endl);

			  Time t_wc_wos = next_job_ready_no_suspension(n,t_min);
			  DM("==== [1] ====" << std::endl);
			  // (1) first check jobs that may be already pending
			  for (const Job<Time>& j : jobs_by_win.lookup(t_min))
			    if (j.earliest_arrival() <= t_min && ready(n, j)) {
			      Time t_high_wos = next_higher_prio_job_ready_no_suspension(n,j, t_min);
			      found_one |= dispatch(n, j, t_wc_wos, t_high_wos);  // ISSUE
			    }

			  DM("==== [2] ====" << std::endl);
			  // (2) check jobs that are released only later in the interval
			  for (auto it = jobs_by_earliest_arrival.upper_bound(t_min);
			       it != jobs_by_earliest_arrival.end();
			       it++) {
			    const Job<Time>& j = *it->second;
			    DM(j << " (" << index_of(j) << ")" << std::endl);
			    // stop looking once we've left the window of interest
			    if (j.earliest_arrival() > t_wc_wos)
			      break;
			    
			    // Job could be not ready due to precedence constraints
			    if (!ready(n, j))
			      continue;

			    // Since this job is released in the future, it better
			    // be incomplete...
			    assert(unfinished(n, j));
			    Time t_high_wos = next_higher_prio_job_ready_no_suspension(n,j, t_min);

			    found_one |= dispatch(n, j, t_wc_wos, t_high_wos);
			  }

			  // check for a dead end
			  if (!found_one && !all_jobs_scheduled(n))
			    // out of options and we didn't schedule all jobs
			    aborted = true;
			}

			// naive: no state merging
			void explore_naively()
			{
				be_naive = true;
				explore();
			}

			void explore()
			{
				make_initial_node(num_cpus);

				while (current_job_count < jobs.size()) {
					unsigned long n;
#ifdef CONFIG_PARALLEL
					const auto& new_nodes_part = nodes_storage.back();
					n = 0;
					for (const Nodes& new_nodes : new_nodes_part) {
						n += new_nodes.size();
					}
#else
					Nodes& exploration_front = nodes();
					n = exploration_front.size();
#endif

					// allocate node space for next depth
					nodes_storage.emplace_back();
					states_storage.emplace_back();

					// keep track of exploration front width
					width = std::max(width, n);

					num_nodes += n;

					check_depth_abort();
					check_cpu_timeout();
					if (aborted)
						break;

#ifdef CONFIG_PARALLEL

					parallel_for(new_states_part.range(),
						[&] (typename Split_nodes::const_range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++) {
								const Nodes& new_nodes = *it;
								auto s = new_nodes.size();
								tbb::parallel_for(tbb::blocked_range<size_t>(0, s),
									[&] (const tbb::blocked_range<size_t>& r) {
										for (size_t i = r.begin(); i != r.end(); i++)
											explore(new_nodes[i]);
								});
							}
						});

#else
					for (const Node& n : exploration_front) {
					  explore(n); //  ISSUE
						check_cpu_timeout();
						if (aborted)
							break;
					}
#endif

					// clean up the state cache if necessary
					if (!be_naive)
						nodes_by_key.clear();

					current_job_count++;

#ifdef CONFIG_PARALLEL
					// propagate any updates to the response-time estimates
					for (auto& r : partial_rta)
						for (const auto& elem : r)
							update_finish_times(rta, elem.first, elem.second);
#endif

#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
					// If we don't need to collect all nodes, we can remove
					// all those that we are done with, which saves a lot of
					// memory.
#ifdef CONFIG_PARALLEL
					parallel_for(nodes_storage.front().range(),
						[] (typename Split_nodes::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					nodes_storage.pop_front();
					states_storage.pop_front();
#endif

				}


#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
				// clean out any remaining nodes
				while (!nodes_storage.empty()) {
#ifdef CONFIG_PARALLEL
					parallel_for(nodes_storage.front().range(),
						[] (typename Split_nodes::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					nodes_storage.pop_front();
				}
#endif


#ifdef CONFIG_PARALLEL
				for (auto &c : edge_counter)
					num_edges += c;
#endif
			}


#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			friend std::ostream& operator<< (std::ostream& out,
											 const State_space<Time>& space)
			{
					std::map<const Schedule_state<Time>*, unsigned int> state_id;
					unsigned int i = 0;
					out << "digraph {" << std::endl;
#ifdef CONFIG_PARALLEL
					for (const Split_states& states : space.get_states()) {
						for (const Schedule_state<Time>& s : tbb::flattened2d<Split_states>(states)) {
#else
					for (const auto& front : space.get_states()) {
						for (const Schedule_state<Time>& s : front) {
#endif
							state_id[&s] = i++;
							out << "\tS" << state_id[&s]
								<< "[label=\"S" << state_id[&s] << ": ";
							s.print_vertex_label(out, space.jobs);
							out << "\"];" << std::endl;
						}
					}
					for (const auto& e : space.get_edges()) {
						out << "\tS" << state_id[e.source]
							<< " -> "
							<< "S" << state_id[e.target]
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
							out << "S" << state_id[e.target]
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

namespace std
{
	template<class Time> struct hash<NP::Global::Schedule_state<Time>>
	{
		std::size_t operator()(NP::Global::Schedule_state<Time> const& s) const
		{
			return s.get_key();
		}
	};
}


#endif
