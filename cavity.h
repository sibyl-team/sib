// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein

#ifndef CAVITY_H
#define CAVITY_H

#include <numeric>
#include <iterator>


// this function computes dest_i = op_{j\neq i} src_j for all i
// given an associative binary operator op with neutral element init in O(n)
// Note: if op has inverse inv() AND is commutative, this is equivalent to
// op((op_j a_j), inv(a_i))
template<typename T, typename iterator_t, typename C>
T cavity(iterator_t const & beg, iterator_t const &end, iterator_t const & dest, T const & init, C const & op)
{
	if (beg == end)
		return init;
	std::partial_sum(beg, end, dest, op);
	iterator_t dest2 = std::next(dest, std::distance(beg, end) - 1);
	T full = *dest2;
	T right = init;
	for (iterator_t it = std::prev(end); it != beg; --dest2, --it) {
		*dest2 = op(*std::prev(dest2), right);
		right = op(*it, right);
	}
	*dest = right;
	return full;
}



#endif
