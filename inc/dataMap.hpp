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
#include <cassert>
#include "graphData.hpp"
#include "modelStats.hpp"


class dataMap{
	public:
		std::map<int,std::set<int> > fNeighbors;
		std::map<int,std::set<int> > secondNeighbors;
		int dim;
                virtual void dbgcheck() {
                        std::cout << "Inside dataMap object!\n";
                        std::cout << "datpair[0] length = "<< datpair[0].size() << "\n";
                }

		std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> >* datpair;
		std::vector<ModelSelfStatsBase*>* datvert;
		char* gtype;
		ModelPairStatsBase* MyNullPairStat;
		float MLcenterscore(int a, int b, int d);
		float MLdeltascore(int a, int b, int x, int d);
		float FBdeltascore(int a, int b, int x, int d);
		float FBcenterscore(int a, int b, int d);
		std::tr1::unordered_map<int, ModelPairStatsBase*>::iterator temp_inIt;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> >::iterator temp_outIt; 
		bool AddEdge(int u, int v, float x, int d);
		bool has_uv(int u, int v, int d);
		bool AddPair(int u, int v, ModelPairStatsBase* p, int d);
		std::set<int>* neighbors(int i);
		ModelPairStatsBase* getuv_ifhas(int u, int v, int d);
		ModelPairStatsBase* get_uv(int u, int v, int d);
		float get_uvSimple(int u, int v, int d);
		virtual void initialize(graphData* D, int _dim);
		virtual void addMergedData(int a, int b, int c);
		virtual ~dataMap();
		virtual void initVert(unsigned int u);
		virtual std::string DerivedType() { return std::string("dataMap"); }
		bool deleteNeighbors(int a);
};

bool dataMap::deleteNeighbors(int a) {
		int x;
		std::set<int>::iterator intit;
		for (intit = fNeighbors[a].begin(); intit != fNeighbors[a].end(); ++intit) {
			x = *intit;
			fNeighbors[x].erase(a);
			fNeighbors[a].erase(x);
		}
		fNeighbors.erase(a);
		for (intit = secondNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
			x = *intit;
			secondNeighbors[x].erase(a);
			secondNeighbors[a].erase(x);
		}
		secondNeighbors.erase(a);
		return true;
}
	


void dataMap::initVert(unsigned int u) {
	for (int d=0;d<this->dim;d++) {
		assert(datvert[d].size()==u);
		if (gtype[d]=='b') {
			datvert[d].push_back(new BinomialSelfStats);
			this->AddPair(u,u,new BinomialPairStats,d);
		} else if (gtype[d]=='p') {
			datvert[d].push_back(new PoissonSelfStats);
			this->AddPair(u,u,new PoissonPairStats,d);
		} else if (gtype[d]=='w') {
			datvert[d].push_back(new WSelfStats);
			this->AddPair(u,u,new WPairStats,d);
		} else if (gtype[d]=='d') {
			datvert[d].push_back(new DcorrSelfStats);
			this->AddPair(u,u,new DcorrPairStats,d);
		} else if (gtype[d]=='g') {
			datvert[d].push_back(new GaussianSelfStats);
			this->AddPair(u,u,new GaussianPairStats,d);
		} else {
			std::cerr << "Error! "<<gtype[d]<<" graph not recognised while initializing datamap\n";
			throw(0);
		}
	}
}


float dataMap::FBcenterscore(int a, int b, int d) {
	assert(d<this->dim);
	ModelPairStatsBase* p = this->get_uv(a,b,d);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->FBcenterscore(datvert[d][a],datvert[d][b],this->get_uv(a,a,d),this->get_uv(b,b,d));
}

float dataMap::MLcenterscore(int a, int b, int d) {
	assert(d<this->dim);
	ModelPairStatsBase* p = this->get_uv(a,b,d);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->MLcenterscore(datvert[d][a],datvert[d][b],this->get_uv(a,a,d),this->get_uv(b,b,d));
}

float dataMap::FBdeltascore(int a, int b, int x, int d) {
	assert(d<this->dim);
	ModelPairStatsBase* p = this->get_uv(a,b,d);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->FBdeltascore(datvert[d][a],datvert[d][b],datvert[d][x],this->get_uv(a,x,d),this->get_uv(b,x,d));
}


float dataMap::MLdeltascore(int a, int b, int x, int d) {
	assert(d<this->dim);
	ModelPairStatsBase* p = this->get_uv(a,b,d);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->MLdeltascore(datvert[d][a],datvert[d][b],datvert[d][x],this->get_uv(a,x,d),this->get_uv(b,x,d));
}


bool dataMap::AddEdge(int u, int v, float x, int d) {
	// returns True if a new edge was added, false otherwise
	assert(d<this->dim);
	bool ans = 0;
	temp_outIt = datpair[d].find(u);
	if (temp_outIt == datpair[d].end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		datpair[d][u] = nd;
		temp_outIt = datpair[d].find(u);
	}
	temp_inIt = temp_outIt->second.find(v);
	if (temp_inIt == temp_outIt->second.end()) {
		// add a pointer
		ModelPairStatsBase *p;
		if (gtype[d]=='b') p = new BinomialPairStats;
		else if (gtype[d]=='p') p = new PoissonPairStats;
		else if (gtype[d]=='w') p = new WPairStats;
		else if (gtype[d]=='d') p = new DcorrPairStats;	
		else if (gtype[d]=='g') p = new GaussianPairStats;	
		else {std::cout << "bad graph type\n"; throw(0);}
		temp_outIt->second[v] = p;
		temp_inIt = temp_outIt->second.find(v);
		ans = 1;
	}
	temp_inIt->second->AddEdge(x);
	datvert[d][u]->AddOutGoing(x);
	return ans;
}

void dataMap::initialize(graphData* D,int _dim) {
	std::set<int> emptySet;
	this->dim = _dim;
	this->datvert = new std::vector<ModelSelfStatsBase*>[this->dim];
	this->datpair = new std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> >[this->dim];
	this->gtype = new char[this->dim];
	unsigned int u, v;
	for (int d=0; d<this->dim;d++) {
		this->gtype[d] = D[d].gtype;
		if (gtype[d]=='b') MyNullPairStat = new BinomialPairStats;
		else if (gtype[d]=='p') MyNullPairStat = new PoissonPairStats;
		else if (gtype[d]=='w') MyNullPairStat = new WPairStats;
		else if (gtype[d]=='d') MyNullPairStat = new DcorrPairStats;
		else if (gtype[d]=='g') MyNullPairStat = new GaussianPairStats;
		else {std::cout << "bad graph type\n";}
		for (u=0; u != D[d].numV; u++) {
			this->initVert(u);
		}
		std::set<int> emptySet;
		for (std::map<int, graphData::destList*>::iterator it1 (D[d].edgeList.begin()); it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			if (fNeighbors.find(u)==fNeighbors.end()) fNeighbors[u] = emptySet;
			for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
				v = (*it2).first;
				assert(this->AddEdge(u,v,(*it2).second,d));
				if (fNeighbors.find(v)==fNeighbors.end()) fNeighbors[v] = emptySet;
				fNeighbors[u].insert(v);
				fNeighbors[v].insert(u);
			}
		}
	}
	std::set<int>::iterator neighbit;
	for (int u=0; u<D[0].numV; u++) {
		if (secondNeighbors.find(u)==secondNeighbors.end()) secondNeighbors[u] = emptySet;
		for (neighbit = fNeighbors[u].begin(); neighbit != fNeighbors[u].end(); ++neighbit) {
			v = *neighbit;
			set_difference_update(secondNeighbors[u],fNeighbors[v],fNeighbors[u]);
		}
		secondNeighbors[u].erase(u);
	}
}



dataMap::~dataMap() {
	for (int d=0; d<this->dim; d++) {
		for (temp_outIt= datpair[d].begin(); temp_outIt != datpair[d].end(); ++temp_outIt) {
			for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
				delete temp_inIt->second;
			}
		}
		datpair[d].clear();
		std::vector<ModelSelfStatsBase*>::iterator vit(datvert[d].begin());
		for (; vit != datvert[d].end(); ++vit) {
			delete *vit;
		}
		datvert[d].clear();
	}
	delete MyNullPairStat;
	delete[] datvert;
	delete[] datpair;
};

void dataMap::addMergedData(int a, int b, int c) {
	std::set<int>::iterator intit;
	int x;
	for (int d=0; d<this->dim; d++) {
		assert(this->AddPair(c,c, this->get_uv(a,a,d)->Add3(this->get_uv(a,b,d), this->get_uv(b,b,d)),d ) );
		this->datvert[d].push_back(this->datvert[d][a]->Add2(this->datvert[d][b],this->get_uv(a,b,d)));
		for (std::set<int>::iterator intit (fNeighbors[c].begin()) ; intit != fNeighbors[c].end(); ++intit) {
			x = (*intit);
			this->AddPair(c,x,this->MyNullPairStat->Add3(this->get_uv(a,x,d),this->get_uv(b,x,d)),d);
			this->AddPair(x,c,this->MyNullPairStat->Add3(this->get_uv(x,a,d),this->get_uv(x,b,d)),d);
		}
	}
	// also need to update neighbors
	std::set<int> emptySet, tempSet;
	fNeighbors[c] = emptySet;
	secondNeighbors[c] = emptySet;
	set_union_update(fNeighbors[c],fNeighbors[a],fNeighbors[b]);
	fNeighbors[c].erase(a); fNeighbors[c].erase(b); 
	tempSet = emptySet;
	set_union_update(tempSet,secondNeighbors[a],secondNeighbors[b]);
	set_difference_update(secondNeighbors[c],tempSet,fNeighbors[c]);
	secondNeighbors[c].erase(a); secondNeighbors[c].erase(b);
	// add it c to the list of its neighbors as well
	for (intit = fNeighbors[c].begin(); intit != fNeighbors[c].end(); ++intit) {
		x = (*intit);
		fNeighbors[x].insert(c);
	}
	for (intit = secondNeighbors[c].begin(); intit != secondNeighbors[c].end(); ++intit) {
		x = (*intit);
		secondNeighbors[x].insert(c);
	}
	
}
	







float dataMap::get_uvSimple(int u, int v, int d) {
	//NOTE: returns 0 if not present
	temp_outIt = datpair[d].find(u);
	if (temp_outIt!=datpair[d].end()) {
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

ModelPairStatsBase* dataMap::get_uv(int u, int v, int d) {
	//NOTE: returns NULL if not present
	temp_outIt = datpair[d].find(u);
	if (temp_outIt!=datpair[d].end()) {
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
	for (int d=0;d<this->dim;d++) {
		temp_outIt  = datpair[d].find(i);
		temp_inIt = (*temp_outIt).second.begin();
		for(;temp_inIt != (*temp_outIt).second.end(); ++temp_inIt) {
			(*pans).insert( (*temp_inIt).first);
		}
	}
	return pans;
};

bool dataMap::has_uv(int u, int v, int d) {
	temp_outIt = datpair[d].find(u);
	if (temp_outIt != datpair[d].end()) {
		temp_inIt = (*temp_outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			return 1;
		}
	}
	return 0;
};

bool dataMap::AddPair(int u, int v, ModelPairStatsBase* p,int d) {
	temp_outIt = datpair[d].find(u);
	if ( temp_outIt != datpair[d].end()) {
		temp_inIt = (*temp_outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			(*temp_inIt).second = p;
			return 0;
		} else {
			(*temp_outIt).second[v]=p;
		}
	} else {
		std::tr1::unordered_map<int,ModelPairStatsBase*> mnew;
		mnew[v] = p;
		datpair[d][u] = mnew;
	}
	return 1;
};


#endif
