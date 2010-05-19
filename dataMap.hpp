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
		std::map<int, float> degrees; // for weighted
		std::map<int, float> nV; // for binomial model
		std::map<int, float> selfMissing; // for weighted networks 
		bool has_uv(int u, int v); //done
		bool AddPair(int u, int v, float d); //done
		bool Addto(int u, int v, float d); //done
		bool eraseAll(); //done
		bool erase(int u, int v); //done
		std::set<int>* neighbors(int i); //done
		float getuv_ifhas(int u, int v);
		float get_uv(int u, int v);
};



float dataMap::get_uv(int u, int v) {
	//NOTE: returns zero if not present
	if (dat.find(u)!=dat.end()) {
		if (dat[u].find(v)!=dat[u].end()) {
			return dat[u][v];
		} else {
			return 0;
		}

	} else {
		return 0;
	}
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
	if (dat.find(u) != dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			dat[u][v] = d;
			return 0;
		} else {
			dat[u][v]=d;
		}
	} else {
		std::map<int,float> mnew;
		mnew[v] = d;
		dat[u] = mnew;
	}
	return 1;
};


bool dataMap::Addto(int u, int v, float d) {
	if (dat.find(u) != dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			dat[u][v] = dat[u][v]+d;
			return 1;
		} else {
			dat[u][v]=d;
		}
	} else {
		std::map<int,float> mnew;
		mnew[v] = d;
		dat[u] = mnew;
	}
	return 0;
};

bool dataMap::eraseAll() {
	if (dat.empty()) {
		return 0;
	} else {
		dat.erase(dat.begin(), dat.end());
		degrees.erase(degrees.begin(), degrees.end());
		return 1;
	}
};

bool dataMap::erase(int u, int v) {
	if (dat.find(u)!=dat.end()) {
		if (dat[u].find(v) != dat[u].end()) {
			dat[u].erase(v);
			if (dat[u].empty()) dat.erase(u);
			return 1;			
		} else {
			if (dat[u].empty()) dat.erase(u);
			return 0;
		}
	} else {
		return 0;
	}
};

#endif
