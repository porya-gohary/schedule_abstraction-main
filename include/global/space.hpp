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
				Problem p{ jobs, num_cpus };
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static State_space explore(
				const Workload& jobs,
				unsigned int num_cpus)
			{
				Problem p{ jobs, num_cpus };
				Analysis_options o;
				return explore(p, o);
			}

			Interval<Time> get_finish_times(const Job<Time>& j) const
			{
				return get_finish_times(j.get_job_index());
			}

			Interval<Time> get_finish_times(Job_index j) const
			{
				if (!rta[j].valid) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				}
				else {
					return rta[j].rt;
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
			typedef tbb::enumerable_thread_specific< Nodes > Split_nodes;
			typedef std::deque<Split_nodes> Nodes_storage;
#else
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

			// Similar to uni/space.hpp, make Response_times a vector of intervals.

			// typedef std::unordered_map<Job_index, Interval<Time> > Response_times;
			struct Response_time_item {
				bool valid;  // RV: is rt valid for the given Job_index?
				Interval<Time> rt;

				Response_time_item()
					: valid(false)
					, rt(0, 0)
				{
				}
			};
			typedef std::vector<Response_time_item> Response_times;



#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
			std::deque<Edge> edges;
#endif
			// Similar to uni/space.hpp, make rta a 

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
			// successors: for each predecessor, a list of successors are present
				// RV: how is this used? Similar to definitions above, is a const reference useful?
				//     It might be useful to store the Job_index in the suspending_task_list.
			typedef std::vector<const Suspending_Task<Time>*> suspending_task_to_list;
			typedef std::vector<std::pair<Job_index, Interval<Time>>> Successors_list;

			std::vector<suspending_task_to_list> _to_suspending_tasks;
			std::vector<Successors_list> _successors;

			const std::vector<suspending_task_to_list>& to_suspending_tasks;
			const std::vector<Successors_list>& successors;

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
				const Precedence_constraints& dag_edges,
				const Suspending_Tasks& susps,
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
				, rta(jobs.size())
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
				, _successors(jobs.size())
				, to_suspending_tasks(_to_suspending_tasks)
				, successors(_successors)
				, wants_self_suspensions(use_self_suspensions)
			{
				for (auto e : dag_edges) {
					const Job<Time>& from = lookup<Time>(jobs, e.first);
					const Job<Time>& to = lookup<Time>(jobs, e.second);
					_predecessors[to.get_job_index()].push_back(from.get_job_index());
				}

				// RV: initialization of to_suspending_tasks and from_suspending_tasks.
				//     Is st.get_toID() equal to index_of(to) ?
				//     Is _predecessors[] related to to_suspending_tasks[] ? 
				for (const Suspending_Task<Time>& st : susps) {
					_to_suspending_tasks[st.get_toIndex()].push_back(&st);
					_predecessors[st.get_toIndex()].push_back(st.get_fromIndex());
					_successors[st.get_fromIndex()].push_back({ st.get_toIndex(), st.get_suspension() });
				}

				for (const Job<Time>& j : jobs) {
					if (_to_suspending_tasks[j.get_job_index()].size() > 0) {
						_jobs_by_latest_arrival_with_susp.insert({ j.latest_arrival(), &j });
					}
					else {
						_jobs_by_latest_arrival_without_susp.insert({ j.latest_arrival(), &j });
					}
					_jobs_by_earliest_arrival.insert({ j.earliest_arrival(), &j });
					_jobs_by_deadline.insert({ j.get_deadline(), &j });
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

			static Time max_deadline(const Workload& jobs)
			{
				Time dl = 0;
				for (const auto& j : jobs)
					dl = std::max(dl, j.get_deadline());
				return dl;
			}

			void update_finish_times(Response_times& r, const Job_index id,
				Interval<Time> range)
			{
				if (!r[id].valid) {
					r[id].valid = true;
					r[id].rt = range;
				}
				else {
					r[id].rt |= range;
				}
				DM("RTA " << id << ": " << r[id].rt << std::endl);
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
				return (std::size_t)(&j - &(jobs[0]));
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
				auto earliest = new_n.get_last_state()->core_availability().min();

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
								new_node(new_n, j, j.get_job_index(), 0, 0);
							const CoreAvailability empty_cav = {};
							State& next_s = new_state(*new_n.get_last_state(), j.get_job_index(), predecessors_of(j), frange, frange, empty_cav, new_n.get_scheduled_jobs(), successors);
							next.add_state(&next_s);
							num_states++;

							// update response times
							update_finish_times(j, frange);
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
							edges.emplace_back(&j, &new_n, &next, frange);
#endif
							count_edge();
							break;
						}
					}
					else
						// deadlines now after the next earliest finish time
						break;
				}
			}

			void make_initial_node(unsigned num_cores)
			{
				// construct initial state
				nodes_storage.emplace_back();
				Node& n = new_node(num_cores, jobs_by_earliest_arrival.begin()->first, jobs_by_latest_arrival_without_susp.begin()->first);
				State& s = new_state(num_cores);
				n.add_state(&s);
				num_states++;
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
			Node_ref alloc_node(Args&&... args)
			{
				nodes().emplace_back(std::forward<Args>(args)...);
				Node_ref n = &(*(--nodes().end()));  
				
				// make sure we didn't screw up...
				auto njobs = n->number_of_scheduled_jobs();
				assert(
					(!njobs && num_states == 0) // initial state
					|| (njobs == current_job_count + 1) // normal State
					|| (njobs == current_job_count + 2 && aborted) // deadline miss
				);

				return n;
			}

			template <typename... Args>
			State& new_state(Args&&... args)
			{
				return *(new State(std::forward<Args>(args)...));
			}

			template <typename... Args>
			Node& new_node(Args&&... args)
			{
				Node_ref n = alloc_node(std::forward<Args>(args)...);

				DM("new node - global " << n << std::endl);
				auto njobs = n->number_of_scheduled_jobs();
				auto idx = njobs % num_todo_queues;
				todo[idx].push_back(n);
				//RV:  add node to nodes_by_key map.
				cache_node(n);
				num_nodes++;
				width = std::max(width, (unsigned long)todo[idx].size() - 1);
				return *n;
			}

#ifdef CONFIG_PARALLEL
			#warning  "Parallel code is not updated for supernodes."
				// make state available for fast lookup
				void insert_cache_state(States_map_accessor& acc, State_ref s)
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
					}
					else if (states_by_key.insert(acc, s->get_key())) {
						// We created the list, so go ahead and insert our state.
						insert_cache_state(acc, s);
						return s;
					}
					// if we raced with concurrent creation, try again
				}
			}

			// make node available for fast lookup
			void insert_cache_node(Nodes_map_accessor& acc, Node_ref n)
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
					}
					else if (nodes_by_key.insert(acc, n->get_key())) {
						// We created the list, so go ahead and insert our node.
						insert_cache_node(acc, n);
						return n;
					}
					// if we raced with concurrent creation, try again
				}
			}

#else
			void cache_node(Node_ref n)
			{
				// create a new list if needed, or lookup if already existing
				auto res = nodes_by_key.emplace(
					std::make_pair(n->get_key(), Node_refs()));

				auto pair_it = res.first;
				Node_refs& list = pair_it->second;

				list.push_front(n);
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
				return unfinished(n, j) && n.job_ready(predecessors_of(j)) && susp_ready(n, j);
			}

			// when a job has self-suspension, susp_ready() should be called.
			// however, it doesn't check the job_finish_times.
			// susp_ready_at() is not used.
			bool ready(const Node& n, const State& s, const Job<Time>& j) const
			{
				return unfinished(n, j) && n.job_ready(predecessors_of(j)) && susp_ready(n, j);
			}

			// RV: if lookup() is time consuming and jobs is fixed, store the result in e.
			//     susp_ready() is used in assert() calls below, so performance might be relevant.
			bool susp_ready(const Node& n, const Job<Time>& j) const
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
						s.get_finish_times(e->get_fromIndex(), rbounds); // ISSUE: segmentation fault: 

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

				for (auto e : fsusps)
				{
					Interval<Time> rbounds = get_finish_times(e->get_fromIndex());
					if (wants_self_suspensions == PATHWISE_SUSP)
						s.get_finish_times(e->get_fromIndex(), rbounds);

					if (slft <= (rbounds.until() + e->get_maxsus()))
					{
						slft = rbounds.until() + e->get_maxsus();
					}
				}

				return slft;
			}

			bool all_jobs_scheduled(const Node& n) const
			{
				return n.number_of_scheduled_jobs() == jobs.size();
			}

			// assumes j is ready, including susp_ready(), but
		  // time information of predecessors and suspension still needs to be checked.
			Interval<Time> ready_times(const State& s, const Job<Time>& j) const
			{
				assert(ready(j));
				Interval<Time> r = j.arrival_window();
				for (auto pred : predecessors_of(j)) {
					Interval<Time> ft{ 0, 0 };
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
					Interval<Time> ft{ 0, 0 };
					if (!s.get_finish_times(pred, ft))
						ft = get_finish_times(jobs[pred]);
					r.lower_bound(ft.min());
					r.extend_to(ft.max());
				}
				// RV: similar to previous ready_times()
				Time seft = get_seft(s, j);
				Time slft = get_slft(s, j);
				if (r.from() < seft)
					r.lower_bound(seft);
				if (r.until() < slft)
					r.extend_to(slft);
				return r;
			}

			Time latest_ready_time(const State& s, const Job<Time>& j) const
			{
				return ready_times(s, j).max();
			}

			Time latest_ready_time(
				const State& s, Time earliest_ref_ready,
				const Job<Time>& j_hp, const Job<Time>& j_ref) const
			{
				auto rt = ready_times(s, j_hp, predecessors_of(j_ref));
				return std::max(rt.max(), earliest_ref_ready);
			}

			Time earliest_ready_time(const State& s, const Job<Time>& j) const
			{
				return std::max(s.core_availability().min(), j.arrival_window().min());
			}

			// Find next time by which a source job (i.e., a job without predecessors) 
			// of higher priority than the reference_job
			// is certainly released in any state in the node 'n'. 
			Time next_certain_higher_priority_source_job_release(
				const Node& n,
				const Job<Time>& reference_job,
				Time until = Time_model::constants<Time>::infinity()) const
			{
				Time when = until;

				// a higher priority source job cannot be released before 
				// a source job of any priority is released
				Time t_earliest = n.get_next_certain_source_job_release();

				for (auto it = jobs_by_latest_arrival_without_susp.lower_bound(t_earliest);
					it != jobs_by_latest_arrival_without_susp.end(); it++) 
				{
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or not of higher priority
					if (ready(n, j) && j.higher_priority_than(reference_job)) 
					{
						when = j.latest_arrival();
						// Jobs are ordered by latest_arrival, so next jobs are later. 
						// We can thus stop searching.
						break;
					}
				}
				return when;
			}

			// Find next time by which a successor job (i.e., a job with predecessors) 
			// of higher priority than the reference_job
			// is certainly released in system state 's' at or before a time 'max'.
			Time next_certain_higher_priority_successor_job_ready_time(
				const Node& n,
				const State& s,
				const Job<Time>& reference_job,
				Time until = Time_model::constants<Time>::infinity()) const
			{
				auto ready_min = earliest_ready_time(s, reference_job);
				Time when = until;

				// a higer priority successor job cannot be ready before 
				// a successor job of any priority is ready
				Time t_earliest = s.get_earliest_certain_successor_jobs_ready_time();
				for (auto it = jobs_by_latest_arrival_with_susp.lower_bound(t_earliest);
					it != jobs_by_latest_arrival_with_susp.end(); it++) 
				{
					const Job<Time>& j = *(it->second);

					// check if we can stop looking
					if (when < j.latest_arrival())
						break; // yep, nothing can lower 'when' at this point

					// j is not relevant if it is already scheduled or not of higher priority
					if (ready(n, j)	&& j.higher_priority_than(reference_job)) {
						// does it beat what we've already seen?
						when = std::min(when,
							latest_ready_time(s, ready_min, j, reference_job));
						// No break, as later jobs might have less suspension.
					}
				}
				return when;
			}

			// find next time by which a job is certainly ready in system state 's'
			Time next_certain_job_ready_time(const Node& n, const State& s) const
			{
				Time t_ws = s.get_earliest_certain_successor_jobs_ready_time();
				Time t_wos = n.get_next_certain_source_job_release();
				return std::min(t_wos, t_ws);
			}

			// assumes j is ready
			// NOTE: we don't use Interval<Time> here because the
			//       Interval c'tor sorts its arguments.
			std::pair<Time, Time> start_times(
				const State& s, const Job<Time>& j, Time t_wc, Time t_high) const
			{
				auto rt = earliest_ready_time(s, j);  //RV: due to suspension, Node is needed. 
				auto at = s.core_availability().min();
				Time est = std::max(rt, at);

				DM("rt: " << rt << std::endl
					<< "at: " << at << std::endl);

				// auto t_high = next_higher_prio_job_ready(n, s, j, at.min());
				Time lst = std::min(t_wc,
					t_high - Time_model::constants<Time>::epsilon());

				DM("est: " << est << std::endl);
				DM("lst: " << lst << std::endl);

				return { est, lst };
			}

			// Find the earliest possible job release of all jobs in a node except for the ignored job
			Time earliest_possible_job_release(
				const Node& n,
				const Job<Time>& ignored_job)
			{
				DM("      - looking for earliest possible job release starting from: "
					<< n.earliest_job_release() << std::endl);

				for(auto it = jobs_by_earliest_arrival.lower_bound(n.earliest_job_release());
					it != jobs_by_earliest_arrival.end(); 	it++)
				{
					const Job<Time>& j = *(it->second);

					DM("         * looking at " << j << std::endl);

					// skip if it is the one we're ignoring or if it was dispatched already
					if (&j == &ignored_job || !unfinished(n, j))
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

				for (auto it = jobs_by_latest_arrival_without_susp.lower_bound(n.get_next_certain_source_job_release()); 
					it != jobs_by_latest_arrival_without_susp.end(); it++)
				{
					const Job<Time>* jp = it->second;
					DM("         * looking at " << *jp << std::endl);

					// skip if it is the one we're ignoring or the job was dispatched already
					if (jp == &ignored_job || !unfinished(n, *jp))
						continue;

					DM("         * found it: " << jp->latest_arrival() << std::endl);
					// it's incomplete and not ignored => found the earliest
					return jp->latest_arrival();
				}

				DM("         * No more future releases" << std::endl);
				return Time_model::constants<Time>::infinity();
			}

			bool dispatch(const Node& n, const Job<Time>& j, Time t_wc_wos, Time t_high_wos)
			{				
				// All states in node 'n' for which the job 'j' is eligible will 
				// be added to that same node. 
				// If such a node already exists, we keep a reference to it
				Node_ref next = nullptr;
				DM("--- global:dispatch() " << n << ", " << j << ", " << t_wc_wos << ", " << t_high_wos << std::endl);
				
				const auto* n_states = n.get_states();
				for (State* s : *n_states) 
				{
					// Calculate t_wc and t_high
					Time t_wc = std::max(s->core_availability().max(), next_certain_job_ready_time(n, *s));
					Time t_high_ws = next_certain_higher_priority_successor_job_ready_time(n, *s, j, t_wc+1);
					Time t_high = std::min(t_high_ws, t_high_wos);

					DM("=== t_high = " << t_high << ", t_wc = " << t_wc << std::endl);
					auto _st = start_times(*s, j, t_wc, t_high);  // RV: no node argument?
					if (_st.first > _st.second)
						return false; // nope, not next job that can be dispatched in state s, return.

					Interval<Time> st{ _st };

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
					
					// If be_naive, a new node and a new state should be created for each new job dispatch.
					if (be_naive)
						next = &(new_node(n, j, index_of(j), earliest_possible_job_release(n, j), earliest_certain_source_job_release(n, j)));
					
					if (next == nullptr)
					{
						const auto pair_it = nodes_by_key.find(n.next_key(j));
						if (pair_it != nodes_by_key.end()) {
							Job_set new_sched_jobs{ n.get_scheduled_jobs(), j.get_job_index() };
							for (Node_ref other : pair_it->second) {
								if (other->get_scheduled_jobs() == new_sched_jobs)
								{
									next = other;
									DM("=== dispatch: next exists." << std::endl);
									break;
								}
							}
						}
					}
					// If there is no node yet, create one.
					if (next == nullptr)
						next = &(new_node(n, j, index_of(j), earliest_possible_job_release(n, j), earliest_certain_source_job_release(n, j)));

					// next should always exist at this point, possibly without states

					// Compute core_avail for the state where j is started.
					CoreAvailability cav;
					s->next_core_avail(j.get_job_index(), predecessors_of(j), st, ftimes, cav);
					// create a new state resulting from scheduling j in state s.
					State& new_s = new_state(*s, index_of(j), predecessors_of(j),
						st, ftimes, cav, next->get_scheduled_jobs(), successors);

					// try to merge the new state with existing states.
					if (next->merge_states(new_s))
						delete &new_s; // if we could merge no need to keep track of the new state anymore
					else {
						next->add_state(&new_s); // else add the new state to the node
						num_states++;
					}

					// make sure we didn't skip any jobs
					// RV: is an easier detection possible?
					//     can detection be done per node (and earliest core_available?)
					if (be_naive) {
						check_for_deadline_misses(n, *next); // ISSUE
					}

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
					// RV: should this be done for every state.
					edges.emplace_back(&j, &s, next, ftimes);
#endif
					count_edge();
				}

				// if we are not using the naive exploration, we check for deadline misses only per job dispatched
				if (!be_naive && next != nullptr)
					check_for_deadline_misses(n, *next);

				return true;
			}

			void explore(const Node& n)
			{
				bool found_one = false;

				DM("---- global:explore(node)" << n.finish_range() << std::endl);
				
				// (0) define the window of interest
				auto t_min = n.earliest_job_release();
				// latest time some unfinished job is certainly ready
				auto nxt_ready_job = n.next_certain_job_ready_time();
				// latest time some core is certainly available
				auto avail_max = n.finish_range().max();
				// latest time by which a work-conserving scheduler
				// certainly schedules some job
				auto upbnd_t_wc = std::max(avail_max, nxt_ready_job);

				DM(n << std::endl);
				DM("t_min: " << t_min << std::endl
					<< "nxt_ready_job: " << nxt_ready_job << std::endl
					<< "avail_max: " << avail_max << std::endl
					<< "upbnd_t_wc: " << upbnd_t_wc << std::endl);

				DM("==== [1] ====" << std::endl);
				// (1) first check jobs that may be already pending
				for (const Job<Time>& j : jobs_by_win.lookup(t_min))
				{
					if (j.earliest_arrival() <= t_min && ready(n, j)) 
					{
						Time t_high_wos = next_certain_higher_priority_source_job_release(n, j, upbnd_t_wc+1);
						// if there is a higher priority job that is certainly ready before job j is released at the earliest, 
						// then j will never be the next job dispached by the scheduler
						if (t_high_wos <= j.earliest_arrival())
							continue;
						found_one |= dispatch(n, j, upbnd_t_wc, t_high_wos);
					}
				}

				DM("==== [2] ====" << std::endl);
				// (2) check jobs that are released only later in the interval
				for (auto it = jobs_by_earliest_arrival.upper_bound(t_min);
					it != jobs_by_earliest_arrival.end();
					it++) 
				{
					const Job<Time>& j = *it->second;
					DM(j << " (" << index_of(j) << ")" << std::endl);
					// stop looking once we've left the window of interest
					if (j.earliest_arrival() > upbnd_t_wc)
						break;

					// Job could be not ready due to precedence constraints
					if (!ready(n, j))
						continue;

					// Since this job is released in the future, it better
					// be incomplete...
					assert(unfinished(n, j));
					
					Time t_high_wos = next_certain_higher_priority_source_job_release(n, j, upbnd_t_wc+1);
					// if there is a higher priority job that is certainly ready before job j is released at the earliest, 
					// then j will never be the next job dispached by the scheduler
					if (t_high_wos <= j.earliest_arrival())
						continue;
					found_one |= dispatch(n, j, upbnd_t_wc, t_high_wos);
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

					// keep track of exploration front width
					width = std::max(width, n);

					check_depth_abort();
					check_cpu_timeout();
					if (aborted)
						break;

#ifdef CONFIG_PARALLEL

					parallel_for(new_states_part.range(),
						[&](typename Split_nodes::const_range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++) {
								const Nodes& new_nodes = *it;
								auto s = new_nodes.size();
								tbb::parallel_for(tbb::blocked_range<size_t>(0, s),
									[&](const tbb::blocked_range<size_t>& r) {
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
						[](typename Split_nodes::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					nodes_storage.pop_front();
#endif

				}


#ifndef CONFIG_COLLECT_SCHEDULE_GRAPH
				// clean out any remaining nodes
				while (!nodes_storage.empty()) {
#ifdef CONFIG_PARALLEL
					parallel_for(nodes_storage.front().range(),
						[](typename Split_nodes::range_type& r) {
							for (auto it = r.begin(); it != r.end(); it++)
								it->clear();
						});
#endif
					nodes_storage.pop_front();
				}
#endif


#ifdef CONFIG_PARALLEL
				for (auto& c : edge_counter)
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
