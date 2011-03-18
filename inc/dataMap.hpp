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
//#include "nodetree.hpp"
//#include "scoremap.hpp"
#include "modelStats.hpp"


class dataMap{
	public:
                virtual void dbgcheck() {
                        std::cout << "Inside dataMap object!\n";
                        std::cout << "datpair length = "<< datpair.size() << "\n";
                }

		std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> > datpair;
		//std::tr1::unordered_map<int, ModelSelfStatsBase*> datvert;
		std::vector<ModelSelfStatsBase*> datvert;
		char gtype;
		ModelPairStatsBase* MyNullPairStat;
		float MLcenterscore(int a, int b);
		float MLdeltascore(int a, int b, int x);
		float FBdeltascore(int a, int b, int x);
		float FBcenterscore(int a, int b);
		std::tr1::unordered_map<int, ModelPairStatsBase*>::iterator temp_inIt;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> >::iterator temp_outIt; 
		bool AddEdge(int u, int v, float x); //done
		bool has_uv(int u, int v); //done
		bool AddPair(int u, int v, ModelPairStatsBase* p);
		//bool Addto(int u, int v, float d);
		//bool eraseAll();
		std::set<int>* neighbors(int i);
		ModelPairStatsBase* getuv_ifhas(int u, int v);
		ModelPairStatsBase* get_uv(int u, int v); //done
		float get_uvSimple(int u, int v); //done
		//virtual bool allErase(int a, int b, int numV); //dont even go there yet, we dont need to erase
		virtual void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbors); //done
		virtual void addMergedData(int a, int b, int c, std::set<int>& fNeighbours);
		virtual ~dataMap();
		virtual void initVert(unsigned int u);
		virtual std::string DerivedType() { return std::string("dataMap"); }
};


void dataMap::initVert(unsigned int u) {
	assert(datvert.size()==u);
	if (gtype=='b') {
		datvert.push_back(new BinomialSelfStats);
		this->AddPair(u,u,new BinomialPairStats);
	} else if (gtype=='p') {
		datvert.push_back(new PoissonSelfStats);
		this->AddPair(u,u,new PoissonPairStats);
	} else if (gtype=='w') {
		datvert.push_back(new WSelfStats);
		this->AddPair(u,u,new WPairStats);
	} else if (gtype=='d') {
		datvert.push_back(new DcorrSelfStats);
		this->AddPair(u,u,new DcorrPairStats);
	} else if (gtype=='g') {
		datvert.push_back(new GaussianSelfStats);
		this->AddPair(u,u,new GaussianPairStats);
	} else {
		std::cerr << "Error! "<<gtype<<" graph not recognised while initializing datamap\n";
		throw(0);
	}
}


float dataMap::FBcenterscore(int a, int b) {
	ModelPairStatsBase* p = this->get_uv(a,b);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->FBcenterscore(datvert[a],datvert[b],this->get_uv(a,a),this->get_uv(b,b));
}

float dataMap::MLcenterscore(int a, int b) {
	ModelPairStatsBase* p = this->get_uv(a,b);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->MLcenterscore(datvert[a],datvert[b],this->get_uv(a,a),this->get_uv(b,b));
}

float dataMap::FBdeltascore(int a, int b, int x) {
	ModelPairStatsBase* p = this->get_uv(a,b);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->FBdeltascore(datvert[a],datvert[b],datvert[x],this->get_uv(a,x),this->get_uv(b,x));
}


float dataMap::MLdeltascore(int a, int b, int x) {
	ModelPairStatsBase* p = this->get_uv(a,b);
	if (p==NULL) {
		p = this->MyNullPairStat;
	}
	return p->MLdeltascore(datvert[a],datvert[b],datvert[x],this->get_uv(a,x),this->get_uv(b,x));
}


bool dataMap::AddEdge(int u, int v, float x) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = datpair.find(u);
	if (temp_outIt == datpair.end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		datpair[u] = nd;
		temp_outIt = datpair.find(u);
	}
	temp_inIt = temp_outIt->second.find(v);
	if (temp_inIt == temp_outIt->second.end()) {
		// add a pointer
		ModelPairStatsBase *p;
		if (gtype=='b') p = new BinomialPairStats;
		else if (gtype=='p') p = new PoissonPairStats;
		else if (gtype=='w') p = new WPairStats;
		else if (gtype=='d') p = new DcorrPairStats;	
		else if (gtype=='g') p = new GaussianPairStats;	
		else {std::cout << "bad graph type\n";}
		temp_outIt->second[v] = p;
		temp_inIt = temp_outIt->second.find(v);
		ans = 1;
	}
	temp_inIt->second->AddEdge(x);
	datvert[u]->AddOutGoing(x);
	return ans;
}

void dataMap::initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
	this->gtype = D.gtype;
	if (gtype=='b') MyNullPairStat = new BinomialPairStats;
	else if (gtype=='p') MyNullPairStat = new PoissonPairStats;
	else if (gtype=='w') MyNullPairStat = new WPairStats;
	else if (gtype=='d') MyNullPairStat = new DcorrPairStats;
	else if (gtype=='g') MyNullPairStat = new GaussianPairStats;
	else {std::cout << "bad graph type\n";}
	int u, v;
	for (u=0; u != D.numV; u++) {
		this->initVert(u);
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
	for (temp_outIt= datpair.begin(); temp_outIt != datpair.end(); ++temp_outIt) {
		for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
			delete temp_inIt->second;
		}
	}
	datpair.clear();
	std::vector<ModelSelfStatsBase*>::iterator vit(datvert.begin());
	for (; vit != datvert.end(); ++vit) {
		delete *vit;
	}
	datvert.clear();
	delete MyNullPairStat;

};

void dataMap::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	assert(this->AddPair(c,c, this->get_uv(a,a)->Add3(this->get_uv(a,b), this->get_uv(b,b)) ) );
	this->datvert.push_back(this->datvert[a]->Add2(this->datvert[b],this->get_uv(a,b)));
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPair(c,x,this->MyNullPairStat->Add3(this->get_uv(a,x),this->get_uv(b,x)));
		this->AddPair(x,c,this->MyNullPairStat->Add3(this->get_uv(x,a),this->get_uv(x,b)));
	}
}
	







float dataMap::get_uvSimple(int u, int v) {
	//NOTE: returns 0 if not present
	temp_outIt = datpair.find(u);
	if (temp_outIt!=datpair.end()) {
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
	temp_outIt = datpair.find(u);
	if (temp_outIt!=datpair.end()) {
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
	temp_outIt  = datpair.find(i);
	if (temp_outIt==datpair.end()) {
		return pans;
	}
	temp_inIt = (*temp_outIt).second.begin();
	for(;temp_inIt != (*temp_outIt).second.end(); ++temp_inIt) {
		(*pans).insert( (*temp_inIt).first);
	}
	return pans;
};

bool dataMap::has_uv(int u, int v) {
	temp_outIt = datpair.find(u);
	if (temp_outIt != datpair.end()) {
		temp_inIt = (*temp_outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			return 1;
		}
	}
	return 0;
};

bool dataMap::AddPair(int u, int v, ModelPairStatsBase* p) {
	temp_outIt = datpair.find(u);
	if ( temp_outIt != datpair.end()) {
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
		datpair[u] = mnew;
	}
	return 1;
};

/*
bool dataMap::Addto(int u, int v, float d) {
	if (datpair.find(u) != datpair.end()) {
		if (dat[u]pair.find(v) != dat[u]pair.end()) {
			datpair[u][v] = datpair[u][v]+d;
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
*/
/*bool dataMap::eraseAll() {
	if (datpair.empty()) {
		return 0;
	} else {
		datpair.erase(datpair.begin(), datpair.end());
		degrees.erase(degrees.begin(), degrees.end());
		return 1;
	}
};*/
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
