#ifndef DATAMAPOTHER_HPP
#define DATAMAPOTHER_HPP

#include "dataMap.hpp"


class dataMapOther: public dataMap{
	public:
		virtual void dbgcheck() {
			std::cout << "Inside dataMapOther object!\n";
			std::cout << "oDatpair length = " << oDatpair.size() << ", datpair length = "  << datpair.size() << "\n";
		}
		dataMapOther() { 
			std::cout << "inside dataMapOther constructor!\n";
		}
		
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> > oDatpair;
		std::vector<ModelSelfStatsBase*> oDatvert;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> > sDatpair;
		std::vector<ModelSelfStatsBase*> sDatvert;
		bool AddPairSoFar(int u, int v, ModelPairStatsBase* p);
		bool AddEdgeSoFar(int u, int v, float x);
		ModelPairStatsBase* get_uvSoFar(int u, int v); // done
		bool AddPairOriginal(int u, int v, ModelPairStatsBase* p); // done
		bool AddEdgeOriginal(int u, int v, float x);
		ModelPairStatsBase* get_uvOriginal(int u, int v); // done
		//virtual bool allErase(int a, int b, int numV); // done
		virtual void initialize(graphData& D, std::map<int,std::set<int> >& fNeighbours) {
			std::cerr << "Wrong initializer function called for the version of dataMap that tracks a second graph (residual)\n";
			throw (1);
		}
		virtual void initialize(graphData& D, graphData& Doriginal, graphData& DsoFar, std::map<int,std::set<int> >& fNeighbours); // done
		virtual void addMergedData(int a, int b, int c, std::set<int>& fNeighbours); // done
		virtual void initVert(unsigned int u);
		virtual ~dataMapOther();
		virtual std::string DerivedType() { return std::string("dataMapOther"); }
};


void dataMapOther::initVert(unsigned int u) {
	dataMap::initVert(u); // call the base class function

	assert(oDatvert.size()==u);
	if (gtype=='b') {
		oDatvert.push_back(new BinomialSelfStats);
		this->AddPairOriginal(u,u,new BinomialPairStats);
	} else if (gtype=='p') {
		oDatvert.push_back(new PoissonSelfStats);
		this->AddPairOriginal(u,u,new PoissonPairStats);
	} else if (gtype=='w') {
		oDatvert.push_back(new WSelfStats);
		this->AddPairOriginal(u,u,new WPairStats);
	} else if (gtype=='d') {
		oDatvert.push_back(new DcorrSelfStats);
		this->AddPairOriginal(u,u,new DcorrPairStats);
	} else if (gtype=='g') {
		oDatvert.push_back(new GaussianSelfStats);
		this->AddPairOriginal(u,u,new GaussianPairStats);
	} else {
		std::cerr << "Error! "<<gtype<<" graph not recognised while initializing oDatamap\n";
		throw(0);
	}

	assert(sDatvert.size()==u);
	if (gtype=='b') {
		sDatvert.push_back(new BinomialSelfStats);
		this->AddPairSoFar(u,u,new BinomialPairStats);
	} else if (gtype=='p') {
		sDatvert.push_back(new PoissonSelfStats);
		this->AddPairSoFar(u,u,new WPairStats);
	} else if (gtype=='w') {
		sDatvert.push_back(new WSelfStats);
		this->AddPairSoFar(u,u,new WPairStats);
	} else if (gtype=='d') {
		sDatvert.push_back(new DcorrSelfStats);
		this->AddPairSoFar(u,u,new DcorrPairStats);
	} else if (gtype=='g') {
		sDatvert.push_back(new GaussianSelfStats);
		this->AddPairSoFar(u,u,new GaussianPairStats);
	} else {
		std::cerr << "Error! "<<gtype<<" graph not recognised while initializing sDatamap\n";
		throw(0);
	}
}

dataMapOther::~dataMapOther() {
	for (temp_outIt= oDatpair.begin(); temp_outIt != oDatpair.end(); ++temp_outIt) {
		for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
			delete temp_inIt->second;
		}
	}
	oDatpair.clear();
	std::vector<ModelSelfStatsBase*>::iterator vit(oDatvert.begin());
	for (; vit != oDatvert.end(); ++vit) {
		delete *vit;
	}
	oDatvert.clear();

	for (temp_outIt= sDatpair.begin(); temp_outIt != sDatpair.end(); ++temp_outIt) {
		for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
			delete temp_inIt->second;
		}
	}
	sDatpair.clear();
	for (vit = sDatvert.begin(); vit != sDatvert.end(); ++vit) {
		delete *vit;
	}
	sDatvert.clear();
};

void dataMapOther::addMergedData(int a, int b, int c, std::set<int>& fNeighbours) {
	assert(this->AddPair(c,c, this->get_uv(a,a)->Add3(this->get_uv(a,b), this->get_uv(b,b)) ) );
	assert(this->AddPairOriginal(c,c, this->get_uvOriginal(a,a)->Add3(this->get_uvOriginal(a,b), this->get_uvOriginal(b,b)) ) );
	assert(this->AddPairSoFar(c,c, this->get_uvSoFar(a,a)->Add3(this->get_uvSoFar(a,b), this->get_uvSoFar(b,b)) ) );
	this->datvert.push_back(this->datvert[a]->Add2(this->datvert[b],this->get_uv(a,b)));
	this->oDatvert.push_back(this->oDatvert[a]->Add2(this->oDatvert[b],this->get_uv(a,b)));
	this->sDatvert.push_back(this->sDatvert[a]->Add2(this->sDatvert[b],this->get_uv(a,b)));
	int x;
	for (std::set<int>::iterator intit (fNeighbours.begin()) ; intit != fNeighbours.end(); ++intit) {
		x = (*intit);
		this->AddPair(c,x,this->MyNullPairStat->Add3(this->get_uv(a,x),this->get_uv(b,x)));
		this->AddPair(x,c,this->MyNullPairStat->Add3(this->get_uv(x,a),this->get_uv(x,b)));

		this->AddPairOriginal(c,x,this->MyNullPairStat->Add3(this->get_uvOriginal(a,x),this->get_uvOriginal(b,x)));
		this->AddPairOriginal(x,c,this->MyNullPairStat->Add3(this->get_uvOriginal(x,a),this->get_uvOriginal(x,b)));

		this->AddPairSoFar(c,x,this->MyNullPairStat->Add3(this->get_uvSoFar(a,x),this->get_uvSoFar(b,x)));
		this->AddPairSoFar(x,c,this->MyNullPairStat->Add3(this->get_uvSoFar(x,a),this->get_uvSoFar(x,b)));
	}
}


bool dataMapOther::AddPairOriginal(int u, int v, ModelPairStatsBase* p) {
	temp_outIt = oDatpair.find(u);
	if ( temp_outIt != oDatpair.end()) {
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
		oDatpair[u] = mnew;
	}
	return 1;
};


bool dataMapOther::AddPairSoFar(int u, int v, ModelPairStatsBase* p) {
	temp_outIt = sDatpair.find(u);
	if ( temp_outIt != sDatpair.end()) {
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
		sDatpair[u] = mnew;
	}
	return 1;
};


bool dataMapOther::AddEdgeOriginal(int u, int v, float x) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = oDatpair.find(u);
	if (temp_outIt == oDatpair.end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		oDatpair[u] = nd;
		temp_outIt = oDatpair.find(u);
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
		else {std::cout << "bad graph type\n"; throw(0);}
		temp_outIt->second[v] = p;
		temp_inIt = temp_outIt->second.find(v);
		ans = 1;
	}
	temp_inIt->second->AddEdge(x);
	oDatvert[u]->AddOutGoing(x);
	return ans;
}

bool dataMapOther::AddEdgeSoFar(int u, int v, float x) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = sDatpair.find(u);
	if (temp_outIt == sDatpair.end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		sDatpair[u] = nd;
		temp_outIt = sDatpair.find(u);
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
		else {std::cout << "bad graph type\n"; throw(0);}
		temp_outIt->second[v] = p;
		temp_inIt = temp_outIt->second.find(v);
		ans = 1;
	}
	temp_inIt->second->AddEdge(x);
	sDatvert[u]->AddOutGoing(x);
	return ans;
}





ModelPairStatsBase* dataMapOther::get_uvOriginal(int u, int v) {
	//NOTE: returns zero if not present
	temp_outIt = oDatpair.find(u);
	if (temp_outIt!=oDatpair.end()) {
		temp_inIt =  (*temp_outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			return (*temp_inIt).second;
		} else {
			return NULL;
		}

	} else {
		return NULL;
	}
};


ModelPairStatsBase* dataMapOther::get_uvSoFar(int u, int v) {
	//NOTE: returns zero if not present
	temp_outIt = sDatpair.find(u);
	if (temp_outIt != sDatpair.end()) {
		temp_inIt =  (*temp_outIt).second.find(v);
		if ( temp_inIt != (*temp_outIt).second.end()) {
			return (*temp_inIt).second;
		} else {
			return NULL;
		}

	} else {
		return NULL;
	}
};


void dataMapOther::initialize(graphData& D, graphData& Doriginal, graphData& DsoFar, std::map<int,std::set<int> >& fNeighbours) {
	std::cout << "inside dataMapOther initialize\n";
	assert(D.gtype==Doriginal.gtype);
	assert(D.gtype==DsoFar.gtype);
	this->gtype = D.gtype;

	if (gtype=='b') MyNullPairStat = new BinomialPairStats;
	else if (gtype=='p') MyNullPairStat = new PoissonPairStats;
	else if (gtype=='w') MyNullPairStat = new WPairStats;
	else if (gtype=='d') MyNullPairStat = new DcorrPairStats;
	else if (gtype=='g') MyNullPairStat = new GaussianPairStats;
	else {std::cout << "bad graph type\n";}

	unsigned int u, v;
	std::set<int> emptySet;
	for (u=0; u != D.numV; u++) {
		this->initVert(u);
	}

	std::map<int, graphData::destList*>::iterator it1 (D.edgeList.begin());
	for (; it1 != D.edgeList.end(); ++it1) {
		u = (*it1).first;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			assert(this->AddEdge(u,v,(*it2).second));
		}
	}

	// Use the Original graph for making first neighbors
	for (it1 = Doriginal.edgeList.begin(); it1 != Doriginal.edgeList.end(); ++it1) {
		u = (*it1).first;
		if (fNeighbours.find(u)==fNeighbours.end()) fNeighbours[u] = emptySet;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			assert(this->AddEdgeOriginal(u,v,(*it2).second));
			if (fNeighbours.find(v)==fNeighbours.end()) fNeighbours[v] = emptySet;
			fNeighbours[u].insert(v);
			fNeighbours[v].insert(u);
		}
	}

	for (it1 = DsoFar.edgeList.begin(); it1 != DsoFar.edgeList.end(); ++it1) {
		u = (*it1).first;
		for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
			v = (*it2).first;
			assert(this->AddEdgeSoFar(u,v,(*it2).second));
		}
	}

}



#endif
