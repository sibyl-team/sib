// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein

#ifndef CAVITY_H
#define CAVITY_H

#include <vector>
#include <iostream>

template<typename T>
class Cavity : public std::vector<T> {
public:


	template<class container_t, typename C>
	Cavity(container_t const & c, T const & init, C const & convolute)
	{
		initialize(c.begin(), c.end(), init, convolute);
	}

	template<class iterator_t, typename C>
	Cavity(iterator_t const & beg, iterator_t const & end, T const & init, C const & convolute)
	{
		initialize(beg, end, init, convolute);
	}
	template<class iterator_t, typename C>
	void initialize(iterator_t const & beg, iterator_t const & end, T const & init, C const & convolute)
	{
		Cavity & me = *this;
		std::vector<T>::resize(distance(beg, end), init);
		int const n = std::vector<T>::size();
		if (n == 0) {
			full_ = init;
			return;
		}

		if (n == 1) {
			me[0] = init;
			full_ = *beg;
			return;
		}

		iterator_t beg2 = beg, end2 = end;
		--end2;
		me[n - 1] = *end2;
		//--end2;
		for (int i = n - 2; i >= 0; --i)
			me[i] = convolute(*(--end2), me[i + 1]);
		full_ = me[0];

		me[0] = me[1];
		T left = *beg2;
		for (int i = 1; i < n - 1; ++i) {
			me[i] = convolute(left, me[i + 1]);
			//std::cout << left.size() << " " << (beg2 + 1)->size() << "  " << i << " " << n << std::endl;
			//assert(left.size() == (beg2 + 1)->size());
			left = convolute(left, *(++beg2));
		}
		me[n - 1] = left;
	}

	T const & full() const { return full_; }

private:
	T full_;
};



#endif
