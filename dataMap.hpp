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
#include "modelStats.hpp"


class dataMap{
	public:
		std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> > datpair;
		//std::tr1::unordered_map<int, ModelSelfStatsBase*> datvert;
		std::vector<ModelSelfStatsBase*> datvert;
		char gtype;

		std::tr1::unordered_map<int, ModelPairStatsBase*>::iterator temp_inIt;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> >::iterator temp_outIt; 
		bool AdEdge(int u, int v, float x); //done
		bool has_uv(int u, int v); //done
		bool AddPair(int u, int v, float d);
		bool Addto(int u, int v, float d);
		bool eraseAll();
		std::set<int>* neighbors(int i);
		ModelPairStatsBase* getuv_ifhas(int u, int v);
		ModelPairStatsBase* get_uv(int u, int v); //done
		float get_uvSimple(int u, int v); //done
		//virtual bool allErase(int a, int b, int numV); //dont even go there yet, we dont need to erase
		virtual void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbors); //done
		virtual void addMergedData(int a, int b, int c, std::set<int>& fNeighbours);
		virtual ~dataMap();
};

bool dataMap::AddEdge(int u, int v, float x) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = datpair.find(u);
	if (temp_outIt == datpair.end()) {
		datpair[u] = std::tr1::unordered_map<int,ModelPairStatsBase*>;
		temp_outIt = datpair.find(u);
	}
	temp_inIt = *temp_outIt.second.find(v);
	if (temp_inIt == temp_outIt.second.end()) {
		// add a pointer
		ModelPairStatsBase *p;
		if (gtype=='b') p = new BinomialPairStats;
		else if (gtype=='p') p = new PoissonPairStats;
		else if (gtype=='w') p = new WPairStats;
		else {cout << "bad graph type\n";
		temp_outIt.second[v] = p;
		temp_inIt = *temp_outIt.second.find(v);
		ans = 1;
	}
	temp_inIt.second->AddEdge(x);
	datvert[u]->AddOutGoing(x);
	return ans;
}

void dataMap::initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
	this->gtype = D.gtype;
	int u, v;
	for (u=0; u != D.numV; u++) {
		assert(datvert.size()==u);
			if (this->gtype=='b') {
				datvert.push_back(new BinomialSelfStats);
			} else if (this->gtype=='p') {
				datvert.push_back(new PoissonSelfStats);
			} else if (this->gtype=='w') {
				datvert.push_back(new WSelfStats);
			} else {
				std::cerr << "Error! "<<this->gtype<<" graph not recognised while initializing datamap\n";
				throw(0);
			}
	}
	std::set<int> emptySet;
	for (std::map<int, graphData::destList*>::iterator it1 (D.edgeList.begin()); it1 != D.edgeList.end(); ++it1) {
		u = (*it1).first;
		if (fNeighbours.find(u)==fNeighbours.end()) fNeighbours[u] = emptySet;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			assert(this->AddEdge(u,v,(*it2).second));
			if (fNeighbours.find(v)==fNeighbours.end()) fNeighbours[v] = emptySet;
			fNeighbours[u].insert(v);
			fNeighbours[v].insert(u);
		}
	}
}



dataMap::~dataMap() {
	for (temp_outIt= datpair.begin(); temp_outIt != datpair.end(); ++tempoutIt) {
		for (temp_inIt = *temp_outIt.second.begin(); temp_inIt != *temp_outIt.second.end(); ++temp_inIt) {
			delete *temp_inIt.second;
		}
	}
	datpair.clear();
	std::tr1::unordered_map<int, ModelSelfStatsBase*>::iterator vit(datvert.begin());
	for (; vit != datvert.end(); ++vit) {
		delete *vit.second;
	}
	datvert.clear();

};

void dataMap::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	assert(this->AddPair(c,c,this->get_uv(a,a) + this->get_uv(b,b) + this->get_uv(a,b)));
	this->nV.push_back(this->nV[a]+this->nV[b]);
	this->degrees[c] = this->degrees[a] + this->degrees[b];
	this->selfMissing[c] = this->selfMissing[a] + this->selfMissing[b] + (this->degrees[a] * this->degrees[b]) - this->get_uv(a,b);
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPair(c,x,this->get_uv(a,x) + this->get_uv(b,x));
		this->AddPair(x,c,this->get_uv(x,a) + this->get_uv(b,x));
	}
}
	




/*bool dataMap::allErase(int a, int b, int numV) {
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
	if (a>=numV) { // why check for this? because we need nv etc. to compute scores for the vertex edges
		this->degrees.erase(a);
		//this->nV.erase(a);
		this->selfMissing.erase(a);
	}
	if (b>=numV) {
		this->degrees.erase(b);
		//this->nV.erase(b);
		this->selfMissing.erase(b);
	}
	//this->degrees.erase(a); this->degrees.erase(b);
	//this->nV.erase(a); this->nV.erase(b);
	//this->selfMissing.erase(a); this->selfMissing.erase(b);
	return 1;
}*/



float dataMap::getDegree(int i) {
	return degrees[i];
};



float dataMap::get_uvSimple(int u, int v) {
	//NOTE: returns NULL if not present
	temp_outIt = dat.find(u);
	if (temp_outIt!=dat.end()) {
		temp_inIt = (*temp_outIt).second.find(v);
		if ( temp_inIt !=(*temp_outIt).second.end()) {
			return (*temp_inIt).second->simple();
		} else {
			return 0;
		}

	} else {
		return 0;
	}
};

ModelPairStatsBase* dataMap::get_uv(int u, int v) {
	//NOTE: returns NULL if not present
	temp_outIt = dat.find(u);
	if (temp_outIt!=dat.end()) {
		temp_inIt = (*temp_outIt).second.find(v);
		if ( temp_inIt !=(*temp_outIt).second.end()) {
			return (*temp_inIt).second;
		} else {
			return NULL;
		}

	} else {
		return NULL;
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
	tmep_outIt = dat.find(u);
	if (temp_outIt != dat.end()) {
		temp_inIt = (*outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			return 1;
		}
	}
	return 0;
};

bool dataMap::AddPair(int u, int v, float d) {
	temp_outIt = dat.find(u);
	if ( temp_outIt != dat.end()) {
		temp_inIt = (*outIt).second.find(v);
		if ( temp_inIt != (*outIt).second.end()) {
			(*temp_inIt).second = d;
			return 0;
		} else {
			(*temp_outIt).second[v]=d;
		}
	} else {
		std::tr1::unordered_map<int,ModelStatsPairBase*> mnew;
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
/*
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
};*/

#endif
