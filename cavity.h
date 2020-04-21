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
		std::vector<T>::resize(distance(beg, end), init);
		int const n = std::vector<T>::size();
		if (n == 0) {
			full_ = init;
			return;
		}

		if (n == 1) {
			std::vector<T>::operator[](0) = init;
			full_ = *beg;
			return;
		}

		std::vector<T> right(n);
		iterator_t beg2 = beg, end2 = end;
		--end2;
		right[n - 1] = *end2;
		//--end2;
		for (int i = n - 2; i >= 0; --i)
			right[i] = convolute(*(--end2), right[i + 1]);
		full_ = right[0];

		std::vector<T>::operator[](0) = right[1];
		T left = *beg2;
		for (int i = 1; i < n - 1; ++i) {
			std::vector<T>::operator[](i) = convolute(left, right[i + 1]);
			//std::cout << left.size() << " " << (beg2 + 1)->size() << "  " << i << " " << n << std::endl;
			//assert(left.size() == (beg2 + 1)->size());
			left = convolute(left, *(++beg2));
		}
		std::vector<T>::operator[](n - 1) = left;
	}

	T const & full() const { return full_; }

private:
	T full_;
};



#endif
