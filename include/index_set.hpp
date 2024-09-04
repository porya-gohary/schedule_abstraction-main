#ifndef INDEX_SET_H
#define INDEX_SET_H

namespace NP {

	class Index_set
	{
	public:

		// Using uint8_t to store 8 bits in each byte
		typedef std::vector<uint8_t> Set_type;

		// new empty job set
		Index_set() : the_set() {}

		// derive a new set by "cloning" an existing set and adding an index
		Index_set(const Index_set& from, std::size_t idx)
				: the_set(std::max(from.the_set.size(), (idx / 8) + 1), 0)
		{
			// Copy the existing set
			std::copy(from.the_set.begin(), from.the_set.end(), the_set.begin());
			set_bit(idx, true);
		}

		// create the diff of two job sets (intended for debugging only)
		Index_set(const Index_set &a, const Index_set &b)
				: the_set(std::max(a.the_set.size(), b.the_set.size()), 0)
		{
			auto limit = std::min(a.the_set.size(), b.the_set.size());
			for (std::size_t i = 0; i < limit * 8; ++i) {
				set_bit(i, a.contains(i) ^ b.contains(i));
			}
		}

		bool operator==(const Index_set &other) const
		{
			return the_set == other.the_set;
		}

		bool operator!=(const Index_set &other) const
		{
			return the_set != other.the_set;
		}

		bool contains(std::size_t idx) const
		{
			if (idx / 8 >= the_set.size()) {
				return false;
			}
			return get_bit(idx);
		}

		bool includes(std::vector<std::size_t> indices) const
		{
			for (auto i : indices)
				if (!contains(i))
					return false;
			return true;
		}

		bool is_subset_of(const Index_set& other) const
		{
			for (std::size_t i = 0; i < the_set.size() * 8; i++) {
				if (contains(i) && !other.contains(i)) {
					return false;
				}
			}
			return true;
		}

		std::size_t size() const
		{
			std::size_t count = 0;
			for (std::size_t i = 0; i < the_set.size() * 8; ++i) {
				if (contains(i)) {
					count++;
				}
			}
			return count;
		}

		void add(std::size_t idx)
		{
			if (idx / 8 >= the_set.size()) {
				the_set.resize((idx / 8) + 1, 0);
			}
			set_bit(idx, true);
		}

		friend std::ostream& operator<< (std::ostream& stream,
										 const Index_set& s)
		{
			bool first = true;
			stream << "{";
			for (std::size_t i = 0; i < s.the_set.size() * 8; ++i) {
				if (s.contains(i)) {
					if (!first)
						stream << ", ";
					first = false;
					stream << i;
				}
			}
			stream << "}";

			return stream;
		}

	private:

		Set_type the_set;

		// Helper functions to set and get individual bits
		void set_bit(std::size_t idx, bool value)
		{
			std::size_t byte_index = idx / 8;
			std::size_t bit_index = idx % 8;
			if (value) {
				the_set[byte_index] |= (1 << bit_index);
			} else {
				the_set[byte_index] &= ~(1 << bit_index);
			}
		}

		bool get_bit(std::size_t idx) const
		{
			std::size_t byte_index = idx / 8;
			std::size_t bit_index = idx % 8;
			return the_set[byte_index] & (1 << bit_index);
		}

		// no accidental copies
		Index_set(const Index_set& origin) = delete;
	};
}

#endif
