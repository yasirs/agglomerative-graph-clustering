#ifndef DATAMAPOTHER_HPP
#define DATAMAPOTHER_HPP

#include "dataMap.hpp"


class dataMapOther: public dataMap{
	public:
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> > oDat;
		std::tr1::unordered_map<int, float> oDegrees; // for weighted
		std::tr1::unordered_map<int, float> oNV; // for binomial model
		std::tr1::unordered_map<int, float> oSelfMissing; // for weighted networks 
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> > sDat;
		std::tr1::unordered_map<int, float> sDegrees; // fsr weighted
		std::tr1::unordered_map<int, float> sNV; // fsr binsmial msdel
		std::tr1::unordered_map<int, float> sSelfMissing; // fsr weighted netwsrks 
		bool AddPairSoFar(int u, int v, float d); // done
		bool AddtoSoFar(int u, int v, float d); // done
		float get_uvSoFar(int u, int v); // done
		bool AddPairOriginal(int u, int v, float d); // done
		bool AddtoOriginal(int u, int v, float d); // done
		float get_uvOriginal(int u, int v); // done
		virtual bool allErase(int a, int b, int numV); // done
		virtual void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
			std::cerr << "Wrong initializer function called for the version of dataMap that tracks a second graph (residual)\n";
			throw (1);
		}
		virtual void initialize(graphData& D, graphData& Doriginal, graphData& DsoFar, std::map<int,std::set<int> >& fNeighbours); // done
		virtual void addMergedData(int a, int b, int c, std::set<int>& fNeighbours); // done
		virtual ~dataMapOther();
};


dataMapOther::~dataMapOther() {
	oDat.clear();
	oDegrees.clear();
	oNV.clear();
	oSelfMissing.clear();
	sDat.clear();
	sDegrees.clear();
	sNV.clear();
	sSelfMissing.clear();
};
	

void dataMapOther::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPairOriginal(c,x,this->get_uvOriginal(a,x) + this->get_uvOriginal(b,x));
		this->AddPairOriginal(x,c,this->get_uvOriginal(x,a) + this->get_uvOriginal(b,x));
		this->AddPairSoFar(c,x,this->get_uvSoFar(a,x) + this->get_uvSoFar(b,x));
		this->AddPairSoFar(x,c,this->get_uvSoFar(x,a) + this->get_uvSoFar(b,x));
		this->AddPair(c,x,this->get_uv(a,x) + this->get_uv(b,x));
		this->AddPair(x,c,this->get_uv(x,a) + this->get_uv(b,x));
	}
}


bool dataMapOther::AddtoOriginal(int u, int v, float d) {
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



bool dataMapOther::AddtoSoFar(int u, int v, float d) {
	if (sDat.find(u) != sDat.end()) {
		if (sDat[u].find(v) != sDat[u].end()) {
			sDat[u][v] = sDat[u][v]+d;
			return 1;
		} else {
			sDat[u][v]=d;
		}
	} else {
		std::tr1::unordered_map<int,float> mnew;
		mnew[v] = d;
		sDat[u] = mnew;
	}
	return 0;
};


bool dataMapOther::AddPairOriginal(int u, int v, float d) {
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


bool dataMapOther::AddPairSoFar(int u, int v, float d) {
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(sDat.find(u));
	if ( outIt != sDat.end()) {
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
		sDat[u] = mnew;
	}
	return 1;
};




float dataMapOther::get_uvOriginal(int u, int v) {
	//NOTE: returns zero if not present
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(oDat.find(u));
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


float dataMapOther::get_uvSoFar(int u, int v) {
	//NOTE: returns zero if not present
	std::tr1::unordered_map<int, std::tr1::unordered_map<int, float> >::iterator outIt(sDat.find(u));
	if (outIt!=sDat.end()) {
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
			this->sDat[x].erase(a);
		}
	}
	this->dat.erase(a);
	this->oDat.erase(a);
	this->sDat.erase(a);
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
			this->sDat[x].erase(b);
		}
	}
	this->dat.erase(b);
	this->oDat.erase(b);
	this->sDat.erase(b);
	if (a>=numV) { // why check for this? I forgot
		this->oDegrees.erase(a);
		this->oNV.erase(a);
		this->oSelfMissing.erase(a);
		this->sDegrees.erase(a);
		this->sNV.erase(a);
		this->sSelfMissing.erase(a);
		this->degrees.erase(a);
		//this->nV.erase(a);
		this->selfMissing.erase(a);
	}
	if (b>=numV) {
		this->oDegrees.erase(b);
		this->oNV.erase(b);
		this->oSelfMissing.erase(b);
		this->sDegrees.erase(b);
		this->sNV.erase(b);
		this->sSelfMissing.erase(b);
		this->degrees.erase(b);
		//this->nV.erase(b);
		this->selfMissing.erase(b);
	}
	/*this->degrees.erase(a); this->degrees.erase(b);
	this->nV.erase(a); this->nV.erase(b);
	this->selfMissing.erase(a); this->selfMissing.erase(b);*/
	return 1;
}


void dataMapOther::initialize(graphData& D, graphData& Doriginal, graphData& DsoFar, std::map<int,std::set<int> >& fNeighbours) {
	int u, v;
	std::set<int> emptySet;
	this->nV.resize(D.numV,0);
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
	for (std::map<int, graphData::destList*>::iterator it1 (Doriginal.edgeList.begin()); it1 != Doriginal.edgeList.end(); ++it1) {
		u = (*it1).first;
		this->nV[u]=1;
		if (this->oDegrees.find(u)==this->oDegrees.end()) this->oDegrees[u]=0;
		if (this->oSelfMissing.find(u)==this->oSelfMissing.end()) this->oSelfMissing[u]=0;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			this->oNV[v]=1;
			this->AddPairOriginal(u,v,(*it2).second);
			this->AddPairOriginal(v,u,(*it2).second);
			if (this->oDegrees.find(v)==this->oDegrees.end()) this->oDegrees[v]=0;
			if (this->oSelfMissing.find(v)==this->oSelfMissing.end()) this->oSelfMissing[v]=0;
			this->oDegrees[u] += (*it2).second;
			this->oDegrees[v] += (*it2).second;
		}
	}
	for (std::map<int, graphData::destList*>::iterator it1 (DsoFar.edgeList.begin()); it1 != DsoFar.edgeList.end(); ++it1) {
		u = (*it1).first;
		this->nV[u]=1;
		if (this->sDegrees.find(u)==this->sDegrees.end()) this->sDegrees[u]=0;
		if (this->sSelfMissing.find(u)==this->sSelfMissing.end()) this->sSelfMissing[u]=0;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			this->sNV[v]=1;
			this->AddPairSoFar(u,v,(*it2).second);
			this->AddPairSoFar(v,u,(*it2).second);
			if (this->sDegrees.find(v)==this->sDegrees.end()) this->sDegrees[v]=0;
			if (this->sSelfMissing.find(v)==this->sSelfMissing.end()) this->sSelfMissing[v]=0;
			this->sDegrees[u] += (*it2).second;
			this->sDegrees[v] += (*it2).second;
		}
	}
}



#endif
