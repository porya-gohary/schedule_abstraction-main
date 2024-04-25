#ifndef GLOBAL_STATE_HPP
#define GLOBAL_STATE_HPP
#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>

#include <set>

#include "util.hpp"
#include "index_set.hpp"
#include "jobs.hpp"
#include "cache.hpp"
#include "statistics.hpp"

namespace NP {

	namespace Global {

		typedef Index_set Job_set;
		typedef std::vector<Job_index> Job_precedence_set;

		template<class Time> class Schedule_node;

		template<class Time> class Schedule_state
		{
		private:

			typedef std::vector<std::pair<Job_index, Interval<Time>>> JobFinishTimes;
			typedef std::vector<Interval<Time>> CoreAvailability;
			typedef std::vector<std::pair<const Job<Time>*, Interval<Time>>> Successors_list;
			typedef std::vector<Successors_list> Successors;

			// system availability intervals
			CoreAvailability core_avail;

			// keeps track of the earliest time a job with at least one predecessor becomes certainly ready
			Time earliest_certain_successor_jobs_ready_time;

			// imprecise set of certainly running jobs
			std::vector<std::pair<Job_index, Interval<Time>>> certain_jobs;

			// job_finish_times holds the finish times of all the jobs that still have unscheduled successor
			JobFinishTimes job_finish_times;

			// Find the offset in the JobFinishTimes vector where the index j should be located.
			int jft_find(const Job_index j) const
			{
				int start = 0;
				int end = job_finish_times.size();
				while (start < end) {
					int mid = (start + end) / 2;
					if (job_finish_times[mid].first == j)
						return mid;
					else if (job_finish_times[mid].first < j)
						start = mid + 1;  // mid is too small, mid+1 might fit.
					else
						end = mid;
				}
				return start;
			}

			// no accidental copies
			Schedule_state(const Schedule_state& origin) = delete;

		public:

			// initial state -- nothing yet has finished, nothing is running
			Schedule_state(unsigned int num_processors)
				: core_avail{ num_processors, Interval<Time>(Time(0), Time(0)) }
				, certain_jobs{}
				, earliest_certain_successor_jobs_ready_time{ Time_model::constants<Time>::infinity() }
			{
				assert(core_avail.size() > 0);
			}

			// transition: new state by scheduling a job in an existing state
			Schedule_state(
				const Schedule_state& from,
				Job_index j,
				const Job_precedence_set& predecessors,
				Interval<Time> start_times,
				Interval<Time> finish_times,
				const CoreAvailability& next_core_avail,
				const Job_set& scheduled_jobs,
				const Successors& successors_of)
			{
				auto est = start_times.min();
				auto lst = start_times.max();
				auto eft = finish_times.min();
				auto lft = finish_times.max();

				DM("est: " << est << std::endl
					<< "lst: " << lst << std::endl
					<< "eft: " << eft << std::endl
					<< "lft: " << lft << std::endl);

				int n_prec = 0;

				// update the set of certainly running jobs
				// keep them sorted to simplify merging
				bool added_j = false;
				for (const auto& rj : from.certain_jobs) 
				{
					auto running_job = rj.first;
					auto running_job_eft = rj.second.min();
					if (contains(predecessors, running_job)) 
					{
						n_prec++; // keep track of the number of predecessors of j that are certainly running
					}
					else if (lst <= running_job_eft) 
					{
						if (!added_j && running_job > j) 
						{
							// right place to add j
							certain_jobs.emplace_back(j, finish_times);
							added_j = true;
						}
						certain_jobs.emplace_back(rj);
					}
				}
				// if we didn't add it yet, add it at the back
				if (!added_j)
					certain_jobs.emplace_back(j, finish_times);

				// calculate the new core availability intervals
				if (next_core_avail.size() > 1)
					core_avail = next_core_avail;
				else
					from.next_core_avail(j, predecessors, start_times, finish_times, core_avail);

				assert(core_avail.size() > 0);

				// update the list of finish times of jobs with successors w.r.t. the previous system state
				// and calculate the earliest time a job with precedence constraints will become ready
				job_finish_times.reserve(from.job_finish_times.size() + 1);
				added_j = false;
				earliest_certain_successor_jobs_ready_time = Time_model::constants<Time>::infinity();
				for (auto ft : from.job_finish_times)
				{
					auto job = ft.first;
					auto job_eft = ft.second.min();
					auto job_lft = ft.second.max();
					// if there is a single core, then we know that 
					// jobs that were disptached in the past cannot have 
					// finished later than when our new job starts executing
					if (core_avail.size() == 1)
					{
						if (job_lft > lst)
							job_lft = lst;
					}

					if (!added_j && job > j)
					{
						bool successor_pending = false;
						for (auto succ : successors_of[j])
						{
							successor_pending = true;
							Time max_susp = succ.second.max();
							earliest_certain_successor_jobs_ready_time = 
								std::min(earliest_certain_successor_jobs_ready_time, 
									std::max(succ.first->latest_arrival(), lft + max_susp));
						}
						if (successor_pending)
							job_finish_times.push_back(std::make_pair(j, finish_times));
						added_j = true;
					}

					bool successor_pending = false;
					for (auto succ : successors_of[job]) {
						auto to_job = succ.first->get_job_index();
						if (!scheduled_jobs.contains(to_job))
						{
							successor_pending = true;
							Time max_susp = succ.second.max();
							earliest_certain_successor_jobs_ready_time = 
								std::min(earliest_certain_successor_jobs_ready_time, 
									std::max(succ.first->latest_arrival(), job_lft + max_susp));
						}
					}
					if (successor_pending)
						job_finish_times.push_back(std::make_pair(job, std::make_pair(job_eft, job_lft)));
				}

				if (!added_j)
				{
					bool successor_pending = false;
					for (auto succ : successors_of[j])
					{
						successor_pending = true;
						Time max_susp = succ.second.max();
						earliest_certain_successor_jobs_ready_time = 
							std::min(earliest_certain_successor_jobs_ready_time, 
								std::max(succ.first->latest_arrival(), lft + max_susp));
					}
					if (successor_pending)
						job_finish_times.push_back(std::make_pair(j, finish_times));
				}

				DM("*** new state: constructed " << *this << std::endl);
			}

			// Return the core availability resulting from scheduling job j 
			// in the current state
			void next_core_avail(Job_index j, const Job_precedence_set& predecessors,
				Interval<Time> start_times, Interval<Time> finish_times,
				CoreAvailability& result) const
			{
				int n_cores = core_avail.size();
				auto est = start_times.min();
				auto lst = start_times.max();
				auto eft = finish_times.min();
				auto lft = finish_times.max();

				DM("est: " << est << std::endl
					<< "lst: " << lst << std::endl
					<< "eft: " << eft << std::endl
					<< "lft: " << lft << std::endl);

				int n_prec = 0;

				// Compute the number of active predecessors.
				// certain_jobs is sorted. If predecessors is sorted, some improvement is possible.
				for (const auto& rj : certain_jobs) {
					auto x = rj.first;
					if (contains(predecessors, x)) {
						n_prec++;
					}
				}

				// compute the cores availability intervals
				std::vector<Time> ca, pa;
				ca.reserve(n_cores);
				pa.reserve(n_cores);

				// Keep pa and ca sorted, by adding the value at the correct place.
				bool eft_added_to_pa = false;
				bool lft_added_to_ca = false;

				// note, we must skip first element in from.core_avail
				if (n_prec > 1) {
					// if there are n_prec predecessors running, n_prec cores must be available when j starts
					for (int i = 1; i < n_prec; i++) {
						pa.push_back(est); // TODO: GN: check whether we can replace by est all the time since predecessors must possibly be finished by est to let j start
						ca.push_back(std::min(lst, std::max(est, core_avail[i].max())));
					}
				}
				else {
					n_prec = 1;
				}
				for (int i = n_prec; i < core_avail.size(); i++) {
					if (!eft_added_to_pa && eft < core_avail[i].min()) 
					{
						pa.push_back(eft);
						eft_added_to_pa = true;
					}
					pa.push_back(std::max(est, core_avail[i].min()));
					if (!lft_added_to_ca && lft < core_avail[i].max()) {
						ca.push_back(lft);
						lft_added_to_ca = true;
					}
					ca.push_back(std::max(est, core_avail[i].max()));
				}
				if (!eft_added_to_pa) pa.push_back(eft);
				if (!lft_added_to_ca) ca.push_back(lft);

				for (int i = 0; i < core_avail.size(); i++) 
				{
					DM(i << " -> " << pa[i] << ":" << ca[i] << std::endl);
					result.emplace_back(pa[i], ca[i]);
				}
			}

			Interval<Time> core_availability() const
			{
				assert(core_avail.size() > 0);
				return core_avail[0];
			}

			Time earliest_finish_time() const
			{
				return core_avail[0].min();
			}

			bool get_finish_times(Job_index j, Interval<Time>& ftimes) const
			{
				int offset = jft_find(j);
				if (offset < job_finish_times.size() && job_finish_times[offset].first == j)
				{
					ftimes = job_finish_times[offset].second;
					return true;
				}
				else {
					return false;
				}
			}

			Time get_earliest_certain_successor_jobs_ready_time() const
			{
				return earliest_certain_successor_jobs_ready_time;
			}

			bool core_avail_overlap(const CoreAvailability& other) const
			{
				assert(core_avail.size() == other.size());
				for (int i = 0; i < core_avail.size(); i++)
					if (!core_avail[i].intersects(other[i]))
						return false;
				return true;
			}

			bool can_merge_with(const Schedule_state<Time>& other, bool useJobFinishTimes = false) const
			{
				if (core_avail_overlap(other.core_avail))
				{
					if (useJobFinishTimes)
						return check_finish_times_overlap(other.job_finish_times);
					else
						return true;
				}
				else
					return false;
			}

			bool can_merge_with(const CoreAvailability& cav, const JobFinishTimes& jft, bool useJobFinishTimes = false) const
			{
				if (core_avail_overlap(cav))
				{
					if (useJobFinishTimes)
						return check_finish_times_overlap(jft);
					else
						return true;
				}
				else
					return false;
			}

			bool try_to_merge(const Schedule_state<Time>& other, bool useJobFinishTimes = false)
			{
				if (!can_merge_with(other, useJobFinishTimes))
					return false;

				merge(other.core_avail, other.job_finish_times, other.certain_jobs, other.earliest_certain_successor_jobs_ready_time);

				DM("+++ merged " << other << " into " << *this << std::endl);
				return true;
			}

			void merge(
				const CoreAvailability& cav, 
				const JobFinishTimes& jft,
				const JobFinishTimes& cert_j,
				Time ecsj_ready_time)
			{
				for (int i = 0; i < core_avail.size(); i++)
					core_avail[i] |= cav[i];

				// vector to collect joint certain jobs
				std::vector<std::pair<Job_index, Interval<Time>>> new_cj;

				// walk both sorted job lists to see if we find matches
				auto it = certain_jobs.begin();
				auto jt = cert_j.begin();
				while (it != certain_jobs.end() &&
					jt != cert_j.end()) {
					if (it->first == jt->first) {
						// same job
						new_cj.emplace_back(it->first, it->second | jt->second);
						it++;
						jt++;
					}
					else if (it->first < jt->first)
						it++;
					else
						jt++;
				}
				// move new certain jobs into the state
				certain_jobs.swap(new_cj);

				// merge job_finish_times
				widen_finish_times(jft);

				// update certain ready time of jobs with predecessors
				earliest_certain_successor_jobs_ready_time = std::max(earliest_certain_successor_jobs_ready_time, ecsj_ready_time);

				DM("+++ merged (cav,jft,cert_t) into " << *this << std::endl);
			}

			friend std::ostream& operator<< (std::ostream& stream,
				const Schedule_state<Time>& s)
			{
				stream << "Global::State(";
				for (const auto& a : s.core_avail)
					stream << "[" << a.from() << ", " << a.until() << "] ";
				stream << "(";
				for (const auto& rj : s.certain_jobs)
					stream << rj.first << "";
				stream << ") " << ")";
				stream << " @ " << &s;
				return stream;
			}

			void print_vertex_label(std::ostream& out,
				const typename Job<Time>::Job_set& jobs) const
			{
				for (const auto& a : core_avail)
					out << "[" << a.from() << ", " << a.until() << "] ";
				out << "\\n";
				bool first = true;
				out << "{";
				for (const auto& rj : certain_jobs) {
					if (!first)
						out << ", ";
					out << "T" << jobs[rj.first].get_task_id()
						<< "J" << jobs[rj.first].get_job_id() << ":"
						<< rj.second.min() << "-" << rj.second.max();
					first = false;
				}
				out << "}";
			}

		private:

			// Check whether the job_finish_times overlap.
			bool check_finish_times_overlap(const JobFinishTimes& from_pwj) const
			{
				bool allJobsIntersect = true;
				// The JobFinishTimes vectors are sorted.
				// Check intersect for matching jobs.
				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();
				while (from_it != from_pwj.end() &&
					state_it != job_finish_times.end())
				{
					if (from_it->first == state_it->first)
					{
						if (!from_it->second.intersects(state_it->second))
						{
							allJobsIntersect = false;
							break;
						}
						from_it++;
						state_it++;
					}
					else if (from_it->first < state_it->first)
						from_it++;
					else
						state_it++;
				}
				return allJobsIntersect;
			}

			void widen_finish_times(const JobFinishTimes& from_pwj)
			{
				// The JobFinishTimes vectors are sorted.
				// Assume check_overlap() is true.
				auto from_it = from_pwj.begin();
				auto state_it = job_finish_times.begin();
				while (from_it != from_pwj.end() &&
					state_it != job_finish_times.end())
				{
					if (from_it->first == state_it->first)
					{
						state_it->second.widen(from_it->second);
						from_it++;
						state_it++;
					}
					else if (from_it->first < state_it->first)
						from_it++;
					else
						state_it++;
				}
			}
		};

		template<class Time> class Schedule_node
		{
		private:

			typedef typename std::vector<Interval<Time>> CoreAvailability;
			Time earliest_pending_release;
			Time next_certain_successor_job_ready_time;
			Time next_certain_source_job_release;

			Job_set scheduled_jobs;
			hash_value_t lookup_key;
			Interval<Time> finish_time;
			unsigned int num_cpus;
			unsigned int num_jobs_scheduled;

			// no accidental copies
			Schedule_node(const Schedule_node& origin) = delete;

			typedef Schedule_state<Time> State;

			struct eft_compare
			{
				bool operator() (State* x, State* y) const
				{
					return x->earliest_finish_time() < y->earliest_finish_time();
				}
			};

			typedef typename std::multiset<State*, eft_compare> State_ref_queue;
			State_ref_queue states;

		public:

			// initial node
			Schedule_node(
				unsigned int num_cores,
				const Time next_earliest_release = 0,
				const Time next_certain_source_job_release = Time_model::constants<Time>::infinity() // the next time a job without predecessor is certainly released
			)
				: lookup_key{ 0 }
				, num_cpus(num_cores)
				, finish_time{ 0,0 }
				, num_jobs_scheduled(0)
				, earliest_pending_release{ next_earliest_release }
				, next_certain_successor_job_ready_time{ Time_model::constants<Time>::infinity() }
				, next_certain_source_job_release{ next_certain_source_job_release }
			{
			}

			// transition: new node by scheduling a job in an existing node
			Schedule_node(
				const Schedule_node& from,
				const Job<Time>& j,
				std::size_t idx,
				const Time next_earliest_release,
				const Time next_certain_source_job_release // the next time a job without predecessor is certainly released
			)
				: scheduled_jobs{ from.scheduled_jobs, idx }
				, lookup_key{ from.next_key(j) }
				, num_cpus(from.num_cpus)
				, num_jobs_scheduled(from.num_jobs_scheduled + 1)
				, finish_time{ 0, Time_model::constants<Time>::infinity() }
				, earliest_pending_release{ next_earliest_release }
				, next_certain_source_job_release{ next_certain_source_job_release }
				, next_certain_successor_job_ready_time{ Time_model::constants<Time>::infinity() }
			{
			}

			~Schedule_node()
			{
				for (State* s : states)
					delete s;
			}

			const unsigned int number_of_scheduled_jobs() const
			{
				return num_jobs_scheduled;
			}

			Time earliest_job_release() const
			{
				return earliest_pending_release;
			}

			Time get_next_certain_source_job_release() const
			{
				return next_certain_source_job_release;
			}

			Time next_certain_job_ready_time() const
			{
				return std::min(next_certain_successor_job_ready_time, next_certain_source_job_release);
			}

			hash_value_t get_key() const
			{
				return lookup_key;
			}

			const Job_set& get_scheduled_jobs() const
			{
				return scheduled_jobs;
			}

			const bool job_incomplete(Job_index j) const
			{
				return !scheduled_jobs.contains(j);
			}

			const bool job_ready(const Job_precedence_set& predecessors) const
			{
				for (auto j : predecessors)
					if (!scheduled_jobs.contains(j))
						return false;
				return true;
			}

			bool matches(const Schedule_node& other) const
			{
				return lookup_key == other.lookup_key &&
					scheduled_jobs == other.scheduled_jobs;
			}

			hash_value_t next_key(const Job<Time>& j) const
			{
				return get_key() ^ j.get_key();
			}

			// RV: finish_range / finish_time contains information about the
			//     earliest and latest core availability for core 0.
			//     whenever a state is changed (through merge) or added,
			//     that interval should be adjusted.
			const Interval<Time>& finish_range() const
			{
				return finish_time;
			}

			/*void update_finish_range(const Interval<Time>& update)
			{
				finish_time.widen(update);
			}*/

			void add_state(State* s)
			{
				// Update finish_time
				Interval<Time> ft = s->core_availability();
				if (states.empty()) {
					finish_time = ft;
					next_certain_successor_job_ready_time = s->get_earliest_certain_successor_jobs_ready_time();
				}
				else {
					finish_time.widen(ft);
					next_certain_successor_job_ready_time = std::max(next_certain_successor_job_ready_time, s->get_earliest_certain_successor_jobs_ready_time());
				}
				states.insert(s);
			}

			friend std::ostream& operator<< (std::ostream& stream,
				const Schedule_node<Time>& n)
			{
				stream << "Node(" << n.states.size() << ")";
				return stream;
			}

			//return the number of states in the node
			int states_size() const
			{
				return states.size();
			}

			const State* get_first_state() const
			{
				auto first = states.begin();
				return *first;
			}

			const State* get_last_state() const
			{
				auto last = --(states.end());
				return *last;
			}

			const State_ref_queue* get_states() const
			{
				return &states;
			}

			bool merge_states(const Schedule_state<Time>& s, bool useJobFinishTimes = false)
			{
				// RV: instead of merging with only one state, try to merge with more states if possible.
				int merge_budget = 1;
				static StatCollect stats = StatCollect("merge");
				stats.tick(merge_budget);

				State* last_state_merged;
				bool result = false;
				for (auto& state : states)
				{
					if (result == false)
					{
						if (state->try_to_merge(s, useJobFinishTimes))
						{
							// Update the node finish_time
							finish_time.widen(s.core_availability());
							//update the certain next job ready time
							next_certain_successor_job_ready_time = std::max(next_certain_successor_job_ready_time, s.get_earliest_certain_successor_jobs_ready_time());

							result = true;

							// Try to merge with a few more states.
							// std::cerr << "Merged with " << merge_budget << " of " << states.size() << " states left.\n";
							merge_budget--;
							if (merge_budget == 0)
								break;

							last_state_merged = state;
						}
					}
					else // if we already merged with one state at least
					{
						if (state->try_to_merge(*last_state_merged, useJobFinishTimes))
						{
							// the state was merged => we can thus remove the old one from the list of states
							states.erase(last_state_merged);
							delete last_state_merged;

							// Try to merge with a few more states.
							// std::cerr << "Merged with " << merge_budget << " of " << states.size() << " states left.\n";
							merge_budget--;
							if (merge_budget == 0)
								break;
							
							last_state_merged = state;
						}
					}
				}

				stats.tick(result);
				stats.print();

				return result;
			}
		};

	}
}

#endif
