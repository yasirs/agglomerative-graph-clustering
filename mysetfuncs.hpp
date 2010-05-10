#ifndef MYSETFUNCS_HPP
#define MYSETFUNCS_HPP
#include <set>

template<typename T> bool set_update(set<T>& target, const set<T>& source) {
	set<T>::iterator it;
	for (it = source.begin(); it != source.end(); ++it) {
		target.insert(*it);
	}
	return 1;
};

template<typename T> bool set_union_update(set<T>& target, const set<T>& s1, const set<T>& s2) {
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), insert_iterator<set<T> >(target,target.begin()));
};


template<typename T> bool set_difference_update(set<T>& target, const set<T>& s1, const set<T>& s2) {
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), insert_iterator<set<T> >(target,target.begin()));
};
	


template<typename T> bool set_symmetric_difference_update(set<T>& target, const set<T>& s1, const set<T>& s2) {
	set_symmetric_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), insert_iterator<set<T> >(target,target.begin()));
};




#endif

