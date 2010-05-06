#ifndef DATAMAP_HPP

#include <map>
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"
#include "scoremap.hpp"


class dataMap{
	public:
		std::map<int, std::map<int,float> > dat;
		std::vector<int> degrees;
		bool has_uv(int u, int v); //done
		bool AddPair(int u, int v, float d); //done
		bool Addto(int u, int v, float d); //done
		bool eraseAll(); //done
		bool erase(int u, int v); //done
		std::set<int>* neighbors(int i); //done
};

std::set<int>* dataMap::neighbors(int i) {
	std::set<int>* pans = new std::set<int>;
	if (dat.find(i)!=dat.end()) {
		return pans;
	}
	std::map<int,float>::iterator it;
	for(it = dat[i].begin();it != dat[i].end(); ++it) {
		(*pans).insert( (*it).first);
	}
	return pans;
};

bool dataMap::has_uv(int u, int v) {
	if (dat.find(u)!= dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			return 1;
		}
	}
	return 0;
};

bool dataMap::AddPair(int u, int v, float d) {
	if (degrees.size()<(u+1)) degrees.resize(u+1);
	if (dat.find(u) != dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			degrees[u] = degrees[u] - dat[u][v] + d;
			dat[u][v] = d;
			return 0;
		} else {
			degrees[u] = degrees[u] + d;
			dat[u][v]=d;
		}
	} else {
		std::map<int,float> mnew;
		mnew[v] = d;
		dat[u] = mnew;
		degrees[u] = degrees[u]+ d;
	}
	return 1;
};


bool dataMap::Addto(int u, int v, float d) {
	if (degrees.size()<(u+1)) degrees.resize(u+1);
	if (dat.find(u) != dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			dat[u][v] = dat[u][v]+d;
			degrees[u] = degrees[u] + d;
			return 1;
		} else {
			degrees[u] = degrees[u] + d;
			dat[u][v]=d;
		}
	} else {
		std::map<int,float> mnew;
		mnew[v] = d;
		dat[u] = mnew;
		degrees[u] = degrees[u] + d;
	}
	return 0;
};

bool dataMap::eraseAll() {
	if (dat.empty()) {
		return 0;
	} else {
		dat.erase(dat.begin(), dat.end());
		degrees.resize(0);
		return 1;
	}
};

bool dataMap::erase(int u, int v) {
	if (dat.find(u)!=dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			degrees[u] = degrees[u] - dat[u][v];
			dat[u].erase(v);
			
		} else {
			return 0;
		}
	} else {
		return 0;
	}
	if (dat[u].empty()) {
		dat.erase(u);
	}
};

#endif
