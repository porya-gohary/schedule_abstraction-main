#ifndef TSN_SPACE_H
#define TSN_SPACE_H

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
#include "clock.hpp"
#include "shaper.hpp"

#include "tsn/state.hpp"

namespace NP {

	namespace TSN {

		template<class Time> class State_space
		{
			public:

			typedef Scheduling_problem<Time> Problem;
			typedef typename Scheduling_problem<Time>::Workload Workload;
			typedef typename Scheduling_problem<Time>::Abort_actions Abort_actions;
			typedef typename Scheduling_problem<Time>::TAS_queues TAS_queues;
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

				auto s = State_space(prob.jobs, prob.dag, prob.aborts, prob.tasQueues,
				                     opts.timeout, opts.max_depth,
				                     opts.num_buckets, opts.early_exit);
				s.cpu_time.start();
				if (opts.be_naive)
					DM("Naive or no merge option currently not supported for TSN");
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
				auto rbounds = rta.find(j.get_id());
				if (rbounds == rta.end()) {
					return Interval<Time>{0, Time_model::constants<Time>::infinity()};
				} else {
					return rbounds->second;
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

			const std::deque<Node>& get_nodes() const
			{
				return nodes;
			}

#endif
			private:

			typedef Job_set Scheduled;

			typedef std::deque<Node> Nodes;

			// All the states created through the making of the schedule abstrcation graph are stored in States.
			// This is done in the function new_state()
			typedef std::deque<State> States;

			// Defining a deque of Interval<Time> so that it might contain the list of intervals when the gate
			// is open or closed. Or might contain a list of intervals that contain the finish time inetrvals 
			// for the new state.
			typedef std::deque<Interval<Time>> Intervals;

			// Iterators are defined for the to iterate over the elemnts in the Nodes and States deque 
			typedef typename std::deque<Node>::iterator Node_ref;
			typedef typename std::deque<State>::iterator State_ref;

			// The Nodes_map typedef allows for storing a key-pair of hash_value_t and Node_ref. The hash_value_t is the
			// typedef defined in jobs.hpp that allows for creating a unique key for each node. A variable called the
			// nodes_by_key is created using Nodes_map. When a new node is made it is added to the map.
			// The use of Nodes_map can be found in new_node(), done_with_current_state() ad schedule().
			typedef std::unordered_multimap<hash_value_t, Node_ref> Nodes_map;

			typedef const Job<Time>* Job_ref;

			// By_time_map stores all the jobs along with a partiular time. there are three variables of type By_time_map.
			// The three variables jobs_by_latest_arrival, jobs_by_earliest_arrival and jobs_by_deadline each store the job
			// along with the latest arrival time, earliest arrival time and the deadline respectively.
			typedef std::multimap<Time, Job_ref> By_time_map;

			// Defining deques of By_time_map so that the by_time_maps can be sorted by priority
			typedef std::deque<By_time_map> By_time_and_priority_map;

			//Add a sorted list of intervals sorted by the start of their interval
			typedef std::multimap<Time,Interval<Time>> Sorted_intervals;

			// The Todo_queue stores the references of all the nodes that are edge nodes of the graph that are yet to
			// be processed to check if the graph can be built further. There are usually 2 Todo_queues defined. This is 
			// because the graph is built in a breadth first manner and so the one todo_queue is for the new breadth that 
			// is made while building the graph and the other todo_queue stores the nodes left to be processed in the 
			// current breadth.
			typedef std::deque<Node_ref> Todo_queue;

			// Response_times stores and updates the finish time intervals of the jobs as the graph is being built
			typedef std::unordered_map<JobID, Interval<Time> > Response_times;

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

			std::vector<Job_precedence_set> job_precedence_sets;

			By_time_map jobs_by_latest_arrival;
			By_time_map jobs_by_earliest_arrival;
			By_time_map jobs_by_deadline;

			By_time_and_priority_map jobs_by_earliest_arrival_priority;
			By_time_and_priority_map jobs_by_latest_arrival_priority;

			std::vector<const Abort_action<Time>*> abort_actions;

			const typename Time_Aware_Shaper<Time>::TAS_set& tasQueues;
			const std::size_t num_fifo_queues;

			Nodes nodes;
			States states;
			unsigned long num_nodes, num_states, num_edges, width;
			Nodes_map nodes_by_key;

			static const std::size_t num_todo_queues = 3;

			//value of the inter_packet_gap is set here
			const Time c_ipg = 0;

			Todo_queue todo[num_todo_queues];
			int todo_idx;
			unsigned long current_job_count;

			Processor_clock cpu_time;
			double timeout;

			unsigned int max_depth;

			bool early_exit;
			bool observed_deadline_miss;

			// Constructor of the State_space class
			State_space(const Workload& jobs,
			            const Precedence_constraints &dag_edges,
			            const Abort_actions& aborts,
			            const TAS_queues& tasQueues,
			            double max_cpu_time = 0,
			            unsigned int max_depth = 0,
			            std::size_t num_buckets = 1000,
			            bool early_exit = true)
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
			, job_precedence_sets(jobs.size())
			, early_exit(early_exit)
			, observed_deadline_miss(false)
			, abort_actions(jobs.size(), NULL)
			, tasQueues(tasQueues)
			, num_fifo_queues(tasQueues.size())
			{
				for (int i = 0; i<num_fifo_queues; i++) {
					By_time_map to_add;
					jobs_by_latest_arrival_priority.emplace_back(to_add);
					jobs_by_earliest_arrival_priority.emplace_back(to_add);
				}
				for (const Job<Time>& j : jobs) {
					jobs_by_earliest_arrival_priority[j.get_priority()].insert({j.earliest_arrival(), &j});
					jobs_by_latest_arrival_priority[j.get_priority()].insert({j.latest_arrival(), &j});
				}
				for (const Job<Time>& j : jobs) {
					jobs_by_latest_arrival.insert({j.latest_arrival(), &j});
					jobs_by_earliest_arrival.insert({j.earliest_arrival(), &j});
					jobs_by_deadline.insert({j.get_deadline(), &j});
				}
				for (auto e : dag_edges) {
					const Job<Time>& from = lookup<Time>(jobs, e.first);
					const Job<Time>& to   = lookup<Time>(jobs, e.second);
					job_precedence_sets[index_of(to)].push_back(index_of(from));
				}
				for (const Abort_action<Time>& a : aborts) {
					const Job<Time>& j = lookup<Time>(jobs, a.get_id());
					abort_actions[index_of(j)] = &a;
				}
			}

			private:

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
				auto rbounds = rta.find(j.get_id());
				if (rbounds == rta.end()) {
					rta.emplace(j.get_id(), range);
					if (j.exceeds_deadline(range.upto()))
						observed_deadline_miss = true;
				} else {
					rbounds->second.widen(range);
					if (j.exceeds_deadline(rbounds->second.upto()))
						observed_deadline_miss = true;
				}
				DM("      New finish time range for " << j
				   << ": " << rta.find(j.get_id())->second << std::endl);

				if (early_exit && observed_deadline_miss)
					aborted = true;
			}

			// The address of j is subtracted from the first element of jobs to find the distance between the two ie., index
			std::size_t index_of(const Job<Time>& j) const
			{
				return (std::size_t) (&j - &(jobs[0]));
			}

			// incomplete() checks to see that a job j jas not been scheduled yet.
			bool incomplete(const Scheduled &scheduled, const Job<Time>& j) const
			{
				return !scheduled.contains(index_of(j));
			}

			bool incomplete(const Node& n, const Job<Time>& j) const
			{
				return incomplete(n.get_scheduled_jobs(), j);
			}

			// find trmax, find the first incomplete job in the list of jobs sorted by latest arrival
			Time get_trmax(const Node& n, Time priority) {
				const Scheduled& already_scheduled = n.get_scheduled_jobs();

			    for (auto it = jobs_by_latest_arrival_priority[priority].begin(); it != jobs_by_latest_arrival_priority[priority].end(); it++) {
			        const Job<Time>& j = *(it->second);

			        // not relevant if already scheduled
			        if (!incomplete(already_scheduled, j))
			            continue;

			        return j.latest_arrival();
			    }

			    return Time_model::constants<Time>::infinity();
			}

			// find if a job of a prioirty exists, in the list of jobs sorted by latest arrival
			bool unscheduled(const Node& n, Time priority) {
				const Scheduled& already_scheduled = n.get_scheduled_jobs();

			    for (auto it = jobs_by_latest_arrival_priority[priority].begin(); it != jobs_by_latest_arrival_priority[priority].end(); it++) {
			        const Job<Time>& j = *(it->second);

			        // not relevant if already scheduled
			        if (!incomplete(already_scheduled, j))
			            continue;

			        return true;
			    }

			    return false;
			}

			// find the t_wc value for TSN
			// for each queue, the maximum between the queue's trmax and the state's earliest certainly available
			// processor time is considered
			// the next instant the gate is open in comparison to this maximum is computed and stored. The minimum of all 
			// the stored values gives the t_wc
			Time 
			calc_t_wc(const Node& n, Time A_max)
			{
				Time t_wc = Time_model::constants<Time>::infinity();
				for(Time i = 0; i<num_fifo_queues; i++)
				{
					if(unscheduled(n,i))
					{
						t_wc = std::min(t_wc, tasQueues[i].next_open(std::max(A_max,get_trmax(n,i))));
						DM("t_wc mid:"<<t_wc<<"\n");
					}
				}
				return t_wc;
			}

// define a couple of iteration helpers

//custom macro/iteration helper for analysing TSN

// Loop through all the jobs that can possibly be on the top of the queue
#define foreach_possibly_fifo_top_job(pftj_macro_local_n, pftj_macro_local_j, pftj_macro_local_priority) \
	for (auto pftj_macro_local_it = jobs_by_earliest_arrival_priority[pftj_macro_local_priority]	\
									.lower_bound((pftj_macro_local_n).earliest_job_release());	\
			pftj_macro_local_it->first <= get_trmax(pftj_macro_local_n, pftj_macro_local_priority)	\
				&& pftj_macro_local_it != jobs_by_earliest_arrival_priority[pftj_macro_local_priority].end() \
				&& (pftj_macro_local_j = pftj_macro_local_it->second);	\
			pftj_macro_local_it++) \
		if(incomplete(pftj_macro_local_n, *pftj_macro_local_j))


// For all the jobs that arrives after the earliest job release of a node n, check if these jobs are still incomplete 
#define foreach_possibly_pending_job(ppj_macro_local_n, ppj_macro_local_j) 	\
	for (auto ppj_macro_local_it = jobs_by_earliest_arrival			\
                     .lower_bound((ppj_macro_local_n).earliest_job_release()); \
	     ppj_macro_local_it != jobs_by_earliest_arrival.end() 				\
	        && (ppj_macro_local_j = ppj_macro_local_it->second); 	\
	     ppj_macro_local_it++) \
		if (incomplete(ppj_macro_local_n, *ppj_macro_local_j))

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
					if (&j == &ignored_job)
						continue;

					DM("         * found it: " << j.earliest_arrival() << std::endl);

					// it's incomplete and not ignored => found the earliest
					return j.earliest_arrival();
				}

				DM("         * No more future releases" << std::endl);

				return Time_model::constants<Time>::infinity();
			}

			void make_initial_node()
			{
				// construct initial node
				new_node();
			}

			// create a new node by adding an elemen to nodes, creating a new state and adding this new state to the node,
			// add this new node to the todo queue as it now a leaf node in the graph, yet to be processed. Add the node
			// along with its key to nodes_by_key. 
			template <typename... Args>
			Node& new_node(Args&&... args)
			{
				DM("\nCreated Node\n");
				nodes.emplace_back(std::forward<Args>(args)...);
				Node_ref n_ref = --nodes.end();

				Interval<Time> fr_ipg = Interval<Time>{n_ref->finish_range().from() + c_ipg, n_ref->finish_range().upto() +c_ipg};

				State &st = (n_ref->get_key()==0) ?
					new_state() :
					new_state(fr_ipg);

				n_ref->add_state(&st);

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
			template <typename... Args>
			State& new_state()
			{
				states.emplace_back();
				State_ref s_ref = --states.end();
				num_states++;
				return *s_ref;
			}

			template <typename... Args>
			State& new_state(Interval<Time> ftimes)
			{
				states.emplace_back(ftimes);
				State_ref s_ref = --states.end();
				num_states++;
				return *s_ref;
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
				assert(nodes.begin() == n);
				nodes.pop_front();
#endif
			}

			//////////////////////////////////////
			// Rules for finding the next state //
			//////////////////////////////////////

			Intervals overlap_delete(Intervals main, Intervals remove)
			{
				// result hold the final result of overlap_delete
				Intervals result;

				Sorted_intervals intermediate;
				intermediate.insert({main[0].from(), main[0]});

				// i is responsible for tracking the elements in main
				// j is responsible for tracking the elements in remove
				Time i = 0;
				Time j = 0;

				while(j < remove.size())
				{
					if(intermediate.size() == 0)
					{
						i += 1;
						if(i >= main.size())
							break;

						intermediate.insert({main[i].from(),main[i]});
					}
					Interval<Time> a = intermediate.begin()->second;
					Interval<Time> b = remove[j];

					if(a.until() < b.from())
					{
						result.emplace_back(intermediate.begin()->second);
						auto it = intermediate.begin();
						intermediate.erase(it);
					}
					else if(b.until() < a.from())
						j+=1;
					else
					{
						auto it = intermediate.begin();
						intermediate.erase(it);

						if(a.from() < b.from())
						{
							Interval<Time> first_split = Interval<Time>{a.from(),b.from()-1};
							intermediate.insert({first_split.from(),first_split});	
						}

						if(a.until() > b.until())
						{
							Interval<Time> second_split = Interval<Time>{b.until()+1,a.until()};
							intermediate.insert({second_split.from(),second_split});	
						}
					}
				}
				for(auto inter: intermediate)
					result.emplace_back(inter.second);
				i+=1;
				while(i < main.size())
				{
					result.emplace_back(main[i]);
					i+=1;
				}
				return result;
			}

			Intervals next_finish_times(const Node &n, const Job<Time> &j, Time est, Time t_wc)
			{
				// CP= contains the times when the gates are closed for 'j' prioirty queue
				// HP= contains the times when the gates are open for all the prioirty queues larger than 'j'
				// ST= contains the interval of possible start times 
				// FT= contains the interval of possible finish times
				DM("Inside next_finish_time#############");
				DM("\n EST = "<<est<<"\t t_wc = "<<t_wc<<"\n");

				Intervals CP, HP, ST, FT;
				ST.emplace_back(Interval<Time>{est, t_wc});
				for(auto st:ST)
					DM("ST:"<<st<<"\n");

				CP = tasQueues[j.get_priority()].get_gates_close(est,t_wc);
				for(auto cp:CP)
					DM("CP:"<<cp<<"\n");
				ST = overlap_delete(ST, CP);
				for(auto st:ST)
					DM("ST:"<<st<<"\n");

				if(ST.size() == 0)
					return FT;

				for(Time i=0; i<j.get_priority(); i+=1)
				{
					if(get_trmax(n, i) > t_wc)
					{
						continue;
					}
					DM("HP entries: "<<get_trmax(n, i)<<" "<<t_wc<<"\n");
					HP =tasQueues[i].get_gates_open(get_trmax(n, i), t_wc);
					DM("HP Loop "<<i<<": ");
					for(auto hp: HP)
						DM(hp<<" ");
					ST = overlap_delete(ST, HP);
				}

				if(ST.size() == 0)
					return FT;

				DM("\nST: ");
				for(auto st: ST)
				{
					DM(st<<" ");
					FT.emplace_back(Interval<Time>{st.from() + j.least_cost(), st.upto() + j.maximal_cost()});
				}

				return FT;
			}

			// creates a new edge from an existing to another node that can either be an existing node
			// in case of merge or a new node when merge is not possible
			void process_new_edge(
				const Node& from,
				const Node& to,
				const Job<Time>& j,
				const Interval<Time>& fr)
			{
				// update statistics
				num_edges++;
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
				edges.emplace_back(&j, &from, &to, fr);
#endif
			}

			// The three schedule functions below are responsible for finding the finish time interval, adding the
			// state to a node and managing the edges for the nodes based on the newly created state

			// In schedule_merge_edge, the function handles the addition of a new state where,
			// not only is the state being merged into another node, but it also has a common edge with another job
			void schedule_merge_edge(const Node& n, const Node_ref match, const Job<Time> &j, Intervals finish_ranges)
			{
				for(auto fr: finish_ranges)
				{
					Interval<Time> fr_ipg = Interval<Time>{fr.from()+c_ipg,fr.upto()+c_ipg};
					if(!match->merge_states(fr_ipg))
					{
						DM("State not merged but added to the node");

						State &st = new_state(fr_ipg);
						match->add_state(&st);
					}
					update_finish_times(j, fr);
				}
			}

			// In schedule_merge_node, the function handles the addition of a new state that is to be merged into another 
			// node but has its own unique edge.
			void schedule_merge_node(const Node& n, const Node_ref match, const Job<Time> &j, Intervals finish_ranges)
			{
				for(auto fr: finish_ranges)
				{
					Interval<Time> fr_ipg = Interval<Time>{fr.from()+c_ipg,fr.upto()+c_ipg};
					if(!match->merge_states(fr_ipg))
					{
						DM("State not merged but added to the node");
						State &st = new_state(fr_ipg);
						match->add_state(&st);
					}
					update_finish_times(j, fr);
				}
				Interval<Time> f_r = Interval<Time>{finish_ranges[0].from(),finish_ranges[finish_ranges.size()-1].until()};
				process_new_edge(n, *match, j, f_r);
			}

			// In schedule_new, the funtion does not attempt to merge but instead creates a new node for the new state
			Node_ref schedule_new(const Node& n, const Job<Time> &j, Intervals finish_ranges)
			{
				DM("Creating a new node"<<std::endl);
				const Node& next =
					new_node(n, j, index_of(j),
					          finish_ranges[0],
					          earliest_possible_job_release(n, j));
				DM("      -----> N" << (nodes.end() - nodes.begin()) << " " <<(todo[todo_idx].front() - nodes.begin() + 1)
				   << std::endl);

				Node_ref n_ref = --nodes.end();


				for(auto fr: finish_ranges)
				{
					if(!n_ref->merge_states(fr))
					{
						DM("State not merged but added to the node");
						State &st = new_state(fr);
						n_ref->add_state(&st);
					}
					update_finish_times(j, fr);
				}

				Interval<Time> f_r = Interval<Time>{finish_ranges[0].from(),finish_ranges[finish_ranges.size()-1].until()};
				process_new_edge(n, next, j, f_r);

				const auto *n_states = next.get_states();
							
				return n_ref;
			}

			void explore()
			{
				make_initial_node();

				while (not_done() && !aborted) 
				{
					const Node& n = next_node();

					DM("\n==================================================="
					   << std::endl);
					DM("Looking at: N"
					   << (todo[todo_idx].front() - nodes.begin() + 1) 
					   << " " << n << std::endl);

					// Obtain all the states in Node n
					const auto *n_states = n.get_states();

					// Keep track of the number of states that have led to new states. This is needed to detect
					// a deadend wherein a state cannot be expanded. The intial value is set to 0 and as we find states that 
					// can be expanded, we increment this variable
					int num_states_expanded = 0;

					const Job<Time>* jp;
					for(int i = 0; i<num_fifo_queues; i++)
					{
						foreach_possibly_fifo_top_job(n, jp, i)
						{
							const Job<Time>& j = *jp;

							DM(j.least_cost());

							// atleast_one_node keeps track of whether atleast one node has been created with the current job
							// All other states that ensure that the job 'j' is eligible will merge into the same node. 
							bool atleast_one_node = false;

							// If atleast_one_node has been found then it has to be stored in the variable match
							Node_ref match;

							DM("\n---\nChecking for states to expand to for job "<< j.get_id()<<std::endl);
							for(State *s: *n_states)
							{
								Time t_wc = calc_t_wc(n,s->latest_finish_time());
								Time est = std::max(j.earliest_arrival(),s->earliest_finish_time());

								if(t_wc < est)
									continue;

								DM("t_wc and est calculated \n");
								Intervals finish_ranges = next_finish_times(n, j, est, t_wc);
								DM("finish ranges calculated \n");

								// check if there are finish ranges present in finish_ranges, if empty then there is no
								// eligible time when j can execute
								if(finish_ranges.size() > 0)
								{
									if(atleast_one_node == false)
									{
										// Find the key of the next state from state s connected by an edge with job j
										// Find all states that have the same key
										auto k = n.next_key(j);
										auto r = nodes_by_key.equal_range(k);

										if(r.first != r.second) {
											Job_set sched_jobs{n.get_scheduled_jobs(), index_of(j)};
											for(auto it = r.first; it != r.second; it++)
											{
												Node &found = *it->second;

												// If the found node and the node that is to be expanded have the same set of scheduled_jobs
												// only then can we find a match to merge
												if(found.get_scheduled_jobs() != sched_jobs)
													continue;

												// If we have reached here, it means that we have found a node where merge is possible.
												match = it->second;
												schedule_merge_node(n, match, j, finish_ranges);
												atleast_one_node = true;
												break;
											}
										}

										// if there exists no existing nodes where the new state can be merged into, then a new node is created
										if(atleast_one_node == false)
										{
											match = schedule_new(n, j, finish_ranges);
											atleast_one_node = true;
										}
									}								
									else
									{
										schedule_merge_edge(n, match, j, finish_ranges); //match is given here
									}

									// If the state is not a deadend, then we already know that the state has been expanded
									// so there is no need to increase the num_states_expanded variable. But if the state
									// is found to be a deadend, then this state can be set to not a deadend because the 
									// state has just been expanded.
									if(s->is_deadend())
									{
										s->not_deadend();
										num_states_expanded++;
									}
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
	template<class Time> struct hash<NP::TSN::Schedule_node<Time>>
    {
		std::size_t operator()(NP::TSN::Schedule_node<Time> const& n) const
        {
            return n.get_key();
        }
    };
}


#endif
