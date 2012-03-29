#ifndef MYSETFUNCS_HPP
#define MYSETFUNCS_HPP
#include <set>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "nodetree.hpp"
#include <iterator>
#include <sstream>
#if ISVC
#include <unordered_set>
#include <unordered_map>
#else
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#endif


/*
#include <float.h>
namespace std {
	int isnan(double x);
}

int std::isnan(double x) {
	return _isnan(x);
}*/

int atoi(const char* a) {
	std::stringstream ss;
	int i;
	ss << a;
	ss >> i;
	return i;
}

int atoi(const std::string& a) {
	std::stringstream ss;
	int i;
	ss << a;
	ss >> i;
	return i;
}


template<typename T> int num_set_common(const std::set<T>& s1, const std::set<T>& s2) {
	std::set<T> target;
	set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target, target.begin()));
	return target.size();
};




template<class K, class V> int num_common_keys(std::tr1::unordered_map<K,V>& m1, std::tr1::unordered_map<K,V>& m2) {
	std::set<K> s;
	get_common_keys(m1,m2,s);
	return s.size();
};


template<class K, class V> void get_common_keys(std::tr1::unordered_map<K,V>& m1, std::tr1::unordered_map<K,V>& m2, std::set<K>& s) {
	typename std::tr1::unordered_map<K,V>::iterator it1;
	typename std::tr1::unordered_map<K,V>::iterator it2;
	K k1, k2;
	for (it1 = m1.begin(); it1 != m1.end(); it1++) {
		k1 = it1->first;
		for (it2 = m2.begin(); it2 != m2.end(); it2++) {
			k2 = it2->first;
			if (k1==k2) s.insert(k1);
		}
	}
};




template<class K, class V> int num_common_keys(std::map<K,V>& m1, std::map<K,V>& m2) {
	std::set<K> s;
	get_common_keys(m1,m2,s);
	return s.size();
};


template<class K, class V> void get_common_keys(std::map<K,V>& m1, std::map<K,V>& m2, std::set<K>& s) {
	typename std::map<K,V>::iterator it1;
	typename std::map<K,V>::iterator it2;
	K k1, k2;
	for (it1 = m1.begin(); it1 != m1.end(); it1++) {
		k1 = it1->first;
		for (it2 = m2.begin(); it2 != m2.end(); it2++) {
			k2 = it2->first;
			if (k1==k2) s.insert(k1);
		}
	}
};

template<typename T> bool set_update(std::set<T>&  target, const std::set<T>&  source) {
	typename std::set<T>::iterator it1;
	for (it1 = source.begin(); it1 != source.end(); ++it1) {
		target.insert(*it1);
	}
	return 1;
};


template<typename T> bool set_union_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	std::insert_iterator<std::set<T> > ii(target,target.begin());
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), ii);
	return 1;
};


template<typename T> class set_union_Enumerator{
	public:
		set_union_Enumerator(std::set<T>& s1, std::set<T>& s2);
		bool next(T& z);
	private:
		typename std::set<T>::iterator x1, x2, last1, last2;
};

template<typename T> set_union_Enumerator<T>::set_union_Enumerator(std::set<T>& s1, std::set<T>& s2) {
	this->x1 = s1.begin(); this->last1 = s1.end();
	this->x2 = s2.begin(); this->last2 = s2.end();
}


template<typename T> bool set_union_Enumerator<T>::next(T& result) {
	if ((x1 != last1)&&(x2 != last2)) {
		if (*x1 < *x2) {
			result = *x1;
			x1++;
			return 1;
		} else if (*x2 < *x1) {
			result = *x2;
			x2++;
			return 1;
		} else {
			result = *x1;
			x1++;
			x2++;
			return 1;
		}
	}
	if (x1 != last1) {
		result = *x1;
		x1++;
		return 1;
	}
	if (x2 != last2) {
		result = *x2;
		x2++;
		return 1;
	}
	return 0;
}

template<typename T> bool set_union_update(std::tr1::unordered_set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_union(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::tr1::unordered_set<T> >(target,target.begin()));
	return 1;
};




template<typename T> bool set_difference_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target,target.begin()));
	return 1;
};
	


template<typename T> bool set_symmetric_difference_update(std::set<T>& target, const std::set<T>& s1, const std::set<T>& s2) {
	set_symmetric_difference(s1.begin(), s1.end(), s2.begin(), s2.end(), std::insert_iterator<std::set<T> >(target,target.begin()));
	return 1;
};



// string functions

std::vector<std::string>* splitspaces(const std::string& str, const std::string& delimiters = " ")
{
    std::vector<std::string>* tokens = new std::vector<std::string>;
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens->push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
};






//TODO comment out this debug function


Node* NfromMap(std::map<int, Node*>& nm, int i) {
	if (nm.find(i)==nm.end()) {
		std::cout << "key " << i << "not present\n";
		throw 1;
	}
	return nm[i];
};


void listintset(std::set<int>& a) {
	std::set<int>::iterator it;
	for (it=a.begin(); it!=a.end(); ++it) {
		std::cout << (*it)<<", ";
	}
	std::cout << "\n";
};


void listintmap(std::map<int, std::set<int> >& a, int i) {
	std::set<int>::iterator it;
	for (it=a[i].begin(); it!=a[i].end(); ++it) {
		std::cout << (*it)<<", ";
	}
	std::cout << "\n";
};

void i2smap(std::map<int, std::string >& a, int i) {
	std::cout << a[i] << "\n";
};

void i2fmap(std::map<int, float>& a, int i) {
	std::cout << a[i] << "\n";
};





#endif

