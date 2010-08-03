#ifndef DATAMAP_HPP
#define DATAMAP_HPP
#if ISVC
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
#include <map>
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"
#include "scoremap.hpp"


class dataMap{
	public:
		std::tr1::unordered_map<int, std::tr1::unordered_map<int,float> > dat;
		std::tr1::unordered_map<int, float> degrees; // for weighted
		std::tr1::unordered_map<int, float> nV; // for binomial model
		std::tr1::unordered_map<int, float> selfMissing; // for weighted networks 
		bool has_uv(int u, int v); //done
		bool AddPair(int u, int v, float d); //done
		bool Addto(int u, int v, float d); //done
		bool eraseAll(); //done
		bool erase(int u, int v); //done
		std::set<int>* neighbors(int i); //done
		float getuv_ifhas(int u, int v);
		float get_uv(int u, int v);
		float getDegree(int i);
		bool allErase(int a, int b, int numV);
		void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbors);
		void addMergedData(int a, int b, int c, std::set<int>& fNeighbours);
};

void dataMap::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPair(c,x,this->get_uv(a,x) + this->get_uv(b,x));
		this->AddPair(x,c,this->get_uv(x,a) + this->get_uv(b,x));
	}
}
	


void dataMap::initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
	int u, v;
	std::set<int> emptySet;
	for (std::map<int, graphData::destList*>::iterator it1 (D.edgeList.begin()); it1 != D.edgeList.end(); ++it1) {
		u = (*it1).first;
		this->nV[u]=1;
		if (fNeighbours.find(u)==fNeighbours.end()) fNeighbours[u] = emptySet;
		if (this->degrees.find(u)==this->degrees.end()) this->degrees[u]=0;
		if (this->selfMissing.find(u)==this->selfMissing.end()) this->selfMissing[u]=0;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			this->nV[v]=1;
			this->AddPair(u,v,(*it2).second);
			this->AddPair(v,u,(*it2).second);
			if (this->degrees.find(v)==this->degrees.end()) this->degrees[v]=0;
			if (this->selfMissing.find(v)==this->selfMissing.end()) this->selfMissing[v]=0;
			this->degrees[u] += (*it2).second;
			this->degrees[v] += (*it2).second;
			if (fNeighbours.find(v)==fNeighbours.end()) fNeighbours[v] = emptySet;
			fNeighbours[u].insert(v);
			fNeighbours[v].insert(u);
		}
	}
}


bool dataMap::allErase(int a, int b, int numV) {
	int x;
	std::tr1::unordered_map<int, std::tr1::unordered_map<int,float> >::iterator outIt(this->dat.find(a));
	if (outIt == this->dat.end()) {
		return 0;
	}
	std::tr1::unordered_map<int,float>::iterator innerIt( (*outIt).second.begin() );
	for (; innerIt != (*outIt).second.end(); ++innerIt) {
		x = (*innerIt).first;
		if ((x != b) and (x != a)) {
			this->dat[x].erase(a);
		}
	}
	this->dat.erase(a);
	outIt = this->dat.find(b);
	if (outIt == this->dat.end()) {
		return 0;
	}
	innerIt=  (*outIt).second.begin();
	for (; innerIt != (*outIt).second.end(); ++innerIt) {
		x = (*innerIt).first;
		if ((x != b) and (x != a)) {
			this->dat[x].erase(b);
		}
	}
	this->dat.erase(b);
	if (a>=numV) { // why check for this? I forgot
		this->degrees.erase(a);
		this->nV.erase(a);
		this->selfMissing.erase(a);
	}
	if (b>=numV) {
		this->degrees.erase(b);
		this->nV.erase(b);
		this->selfMissing.erase(b);
	}
	/*this->degrees.erase(a); this->degrees.erase(b);
	this->nV.erase(a); this->nV.erase(b);
	this->selfMissing.erase(a); this->selfMissing.erase(b);*/
	return 1;
}



float dataMap::getDegree(int i) {
	return degrees[i];
};

float dataMap::get_uv(int u, int v) {
	//NOTE: returns zero if not present
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(dat.find(u));
	if (outIt!=dat.end()) {
		std::tr1::unordered_map<int,float>::iterator inIt ( (*outIt).second.find(v));
		if ( inIt !=(*outIt).second.end()) {
			return (*inIt).second;
		} else {
			return 0;
		}

	} else {
		return 0;
	}
};

std::set<int>* dataMap::neighbors(int i) {
	std::set<int>* pans = new std::set<int>;
	std::tr1::unordered_map<int, std::tr1::unordered_map<int,float> >::iterator outIt (dat.find(i));
	if (outIt==dat.end()) {
		return pans;
	}
	std::tr1::unordered_map<int,float>::iterator it ((*outIt).second.begin());
	for(it = (*outIt).second.begin();it != (*outIt).second.end(); ++it) {
		(*pans).insert( (*it).first);
	}
	return pans;
};

bool dataMap::has_uv(int u, int v) {
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(dat.find(u));
	if (outIt != dat.end()) {
		std::tr1::unordered_map<int,float>::iterator inIt ( (*outIt).second.find(v));
		if ( inIt != (*outIt).second.end()) {
			return 1;
		}
	}
	return 0;
};

bool dataMap::AddPair(int u, int v, float d) {
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(dat.find(u));
	if ( outIt != dat.end()) {
		std::tr1::unordered_map<int,float>::iterator inIt ( (*outIt).second.find(v));
		if ( inIt != (*outIt).second.end()) {
			(*inIt).second = d;
			return 0;
		} else {
			(*outIt).second[v]=d;
		}
	} else {
		std::tr1::unordered_map<int,float> mnew;
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
		std::tr1::unordered_map<int,float> mnew;
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
