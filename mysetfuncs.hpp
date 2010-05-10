#ifndef MYSETFUNCS_HPP
#define MYSETFUNCS_HPP
#include <set>
#include <algorithm>


template<typename T> bool set_update(std::set<T>&  target, const std::set<T>&  source) {
	typename std::set<T>::iterator it1;
	for (it1 = source.begin(); it1 != source.end(); ++it1) {
		target.insert(*it1);
	}
	return 1;
};


template<typename T> bool set_union_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target,target.begin()));
};


template<typename T> bool set_difference_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target,target.begin()));
};
	


template<typename T> bool set_symmetric_difference_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_symmetric_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target,target.begin()));
};




#endif

