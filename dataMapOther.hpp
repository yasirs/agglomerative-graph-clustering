#ifndef DATAMAPOTHER_HPP
#define DATAMAPOTHER_HPP

#include "dataMap.hpp"


class dataMapOther: public dataMap{
	public:
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> > oDat;
		std::tr1::unordered_map<int, float> oDegrees; // for weighted
		std::tr1::unordered_map<int, float> oNV; // for binomial model
		std::tr1::unordered_map<int, float> oSelfMissing; // for weighted networks 
		bool AddPairOther(int u, int v, float d); // done
		bool AddtoOther(int u, int v, float d); // done
		float get_uvOther(int u, int v); // done
		virtual bool allErase(int a, int b, int numV); // done
		virtual void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
			std::cerr << "Wrong initializer function called for the version of dataMap that tracks a second graph (residual)\n";
			throw (1);
		}
		virtual void initialize(graphData& D, graphData& Dother, std::map<int,std::set<int> >& fNeighbours); // done
		virtual void addMergedData(int a, int b, int c, std::set<int>& fNeighnours); // done
};


void dataMapOther::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPairOther(c,x,this->get_uvOther(a,x) + this->get_uvOther(b,x));
		this->AddPairOther(x,c,this->get_uvOther(x,a) + this->get_uvOther(b,x));
		this->AddPair(c,x,this->get_uv(a,x) + this->get_uv(b,x));
		this->AddPair(x,c,this->get_uv(x,a) + this->get_uv(b,x));
	}
}

bool dataMapOther::AddtoOther(int u, int v, float d) {
	if (oDat.find(u) != oDat.end()) {
		if (oDat[u].find(v) != oDat[u].end()) {
			oDat[u][v] = oDat[u][v]+d;
			return 1;
		} else {
			oDat[u][v]=d;
		}
	} else {
		std::tr1::unordered_map<int,float> mnew;
		mnew[v] = d;
		oDat[u] = mnew;
	}
	return 0;
};


bool dataMapOther::AddPairOther(int u, int v, float d) {
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(oDat.find(u));
	if ( outIt != oDat.end()) {
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
		oDat[u] = mnew;
	}
	return 1;
};


float dataMapOther::get_uvOther(int u, int v) {
	//NOTE: returns zero if not present
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(dat.find(u));
	if (outIt!=oDat.end()) {
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

bool dataMapOther::allErase(int a, int b, int numV) {
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
			this->oDat[x].erase(a);
		}
	}
	this->dat.erase(a);
	this->oDat.erase(a);
	outIt = this->dat.find(b);
	if (outIt == this->dat.end()) {
		return 0;
	}
	innerIt=  (*outIt).second.begin();
	for (; innerIt != (*outIt).second.end(); ++innerIt) {
		x = (*innerIt).first;
		if ((x != b) and (x != a)) {
			this->dat[x].erase(b);
			this->oDat[x].erase(b);
		}
	}
	this->dat.erase(b);
	this->oDat.erase(b);
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


void dataMapOther::initialize(graphData& D, graphData& Dother, std::map<int,std::set<int> >& fNeighbours) {
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
	for (std::map<int, graphData::destList*>::iterator it1 (Dother.edgeList.begin()); it1 != Dother.edgeList.end(); ++it1) {
		u = (*it1).first;
		this->nV[u]=1;
		if (this->oDegrees.find(u)==this->oDegrees.end()) this->oDegrees[u]=0;
		if (this->oSelfMissing.find(u)==this->oSelfMissing.end()) this->oSelfMissing[u]=0;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			this->oNV[v]=1;
			this->AddPairOther(u,v,(*it2).second);
			this->AddPairOther(v,u,(*it2).second);
			if (this->oDegrees.find(v)==this->oDegrees.end()) this->oDegrees[v]=0;
			if (this->oSelfMissing.find(v)==this->oSelfMissing.end()) this->oSelfMissing[v]=0;
			this->oDegrees[u] += (*it2).second;
			this->oDegrees[v] += (*it2).second;
			if (fNeighbours.find(v)==fNeighbours.end()) fNeighbours[v] = emptySet;
			fNeighbours[u].insert(v);
			fNeighbours[v].insert(u);
		}
	}
}



#endif
