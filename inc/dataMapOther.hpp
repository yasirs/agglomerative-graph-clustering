#ifndef DATAMAPOTHER_HPP
#define DATAMAPOTHER_HPP

#include "dataMap.hpp"


class dataMapOther: public dataMap{
	public:
		virtual void dbgcheck() {
			std::cout << "Inside dataMapOther object!\n";
			std::cout << "oDatpair[0] length = " << oDatpair[0].size() << ", datpair[0] length = "  << datpair[0].size() << "\n";
		}
		dataMapOther() { 
			std::cout << "inside dataMapOther constructor!\n";
		}

		std::map<int,std::set<int> > oFNeighbors;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> >* oDatpair;
		std::vector<ModelSelfStatsBase*>* oDatvert;
		std::tr1::unordered_map<int, std::tr1::unordered_map<int, ModelPairStatsBase*> >* sDatpair;
		std::vector<ModelSelfStatsBase*>* sDatvert;
		bool AddPairSoFar(int u, int v, ModelPairStatsBase* p, int d);
		bool AddEdgeSoFar(int u, int v, float x, int d);
		ModelPairStatsBase* get_uvSoFar(int u, int v, int d); // done
		bool AddPairOriginal(int u, int v, ModelPairStatsBase* p, int d); // done
		bool AddEdgeOriginal(int u, int v, float x, int d);
		ModelPairStatsBase* get_uvOriginal(int u, int v, int d); // done
		virtual void initialize(graphData* D, int _dim) {
			std::cerr << "Wrong initializer function called for the version of dataMap that tracks a second graph (residual)\n";
			throw (1);
		}
		virtual void initialize(graphData* D, graphData* Doriginal, graphData* DsoFar, int _dim); // done
		virtual void addMergedData(int a, int b, int c); // done
		virtual void initVert(unsigned int u);
		virtual ~dataMapOther();
		virtual std::string DerivedType() { return std::string("dataMapOther"); }
};


void dataMapOther::initVert(unsigned int u) {
	dataMap::initVert(u); // call the base class function
	for(int d=0;d<this->dim;d++) {
		assert(oDatvert[d].size()==u);
		if (gtype[d]=='b') {
			oDatvert[d].push_back(new BinomialSelfStats);
			this->AddPairOriginal(u,u,new BinomialPairStats,d);
		} else if (gtype[d]=='p') {
			oDatvert[d].push_back(new PoissonSelfStats);
			this->AddPairOriginal(u,u,new PoissonPairStats,d);
		} else if (gtype[d]=='w') {
			oDatvert[d].push_back(new WSelfStats);
			this->AddPairOriginal(u,u,new WPairStats,d);
		} else if (gtype[d]=='d') {
			oDatvert[d].push_back(new DcorrSelfStats);
			this->AddPairOriginal(u,u,new DcorrPairStats,d);
		} else if (gtype[d]=='g') {
			oDatvert[d].push_back(new GaussianSelfStats);
			this->AddPairOriginal(u,u,new GaussianPairStats,d);
		} else {
			std::cerr << "Error! "<<gtype[d]<<" graph not recognised while initializing oDatamap\n";
			throw(0);
		}
	
		assert(sDatvert[d].size()==u);
		if (gtype[d]=='b') {
			sDatvert[d].push_back(new BinomialSelfStats);
			this->AddPairSoFar(u,u,new BinomialPairStats,d);
		} else if (gtype[d]=='p') {
			sDatvert[d].push_back(new PoissonSelfStats);
			this->AddPairSoFar(u,u,new WPairStats,d);
		} else if (gtype[d]=='w') {
			sDatvert[d].push_back(new WSelfStats);
			this->AddPairSoFar(u,u,new WPairStats,d);
		} else if (gtype[d]=='d') {
			sDatvert[d].push_back(new DcorrSelfStats);
			this->AddPairSoFar(u,u,new DcorrPairStats,d);
		} else if (gtype[d]=='g') {
			sDatvert[d].push_back(new GaussianSelfStats);
			this->AddPairSoFar(u,u,new GaussianPairStats,d);
		} else {
			std::cerr << "Error! "<<gtype[d]<<" graph not recognised while initializing sDatamap\n";
			throw(0);
		}
	}
}

dataMapOther::~dataMapOther() {
	for (int d=0;d<this->dim;d++) {
		for (temp_outIt= oDatpair[d].begin(); temp_outIt != oDatpair[d].end(); ++temp_outIt) {
			for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
				delete temp_inIt->second;
			}
		}
		oDatpair[d].clear();
		std::vector<ModelSelfStatsBase*>::iterator vit(oDatvert[d].begin());
		for (; vit != oDatvert[d].end(); ++vit) {
			delete *vit;
		}
		oDatvert[d].clear();
	
		for (temp_outIt= sDatpair[d].begin(); temp_outIt != sDatpair[d].end(); ++temp_outIt) {
			for (temp_inIt = temp_outIt->second.begin(); temp_inIt != temp_outIt->second.end(); ++temp_inIt) {
				delete temp_inIt->second;
			}
		}
		sDatpair[d].clear();
		for (vit = sDatvert[d].begin(); vit != sDatvert[d].end(); ++vit) {
			delete *vit;
		}
		sDatvert[d].clear();
	}
	delete[] sDatvert;
	delete[] oDatvert;
	delete[] sDatpair;
	delete[] oDatpair;
};

void dataMapOther::addMergedData(int a, int b, int c) {
	std::set<int>::iterator intit;
	int x;
	int d;
	for(d=0; d< this->dim;d++) {
		assert(this->AddPair(c,c, this->get_uv(a,a,d)->Add3(this->get_uv(a,b,d), this->get_uv(b,b,d)),d ) );
		assert(this->AddPairOriginal(c,c, this->get_uvOriginal(a,a,d)->Add3(this->get_uvOriginal(a,b,d), this->get_uvOriginal(b,b,d)),d ) );
		assert(this->AddPairSoFar(c,c, this->get_uvSoFar(a,a,d)->Add3(this->get_uvSoFar(a,b,d), this->get_uvSoFar(b,b,d)),d ) );
		this->datvert[d].push_back(this->datvert[d][a]->Add2(this->datvert[d][b],this->get_uv(a,b,d)));
		this->oDatvert[d].push_back(this->oDatvert[d][a]->Add2(this->oDatvert[d][b],this->get_uv(a,b,d)));
		this->sDatvert[d].push_back(this->sDatvert[d][a]->Add2(this->sDatvert[d][b],this->get_uv(a,b,d)));
	}
	// also need to update neighbors
	std::set<int> emptySet, tempSet;
	fNeighbors[c] = emptySet; oFNeighbors[c] = emptySet;
	secondNeighbors[c] = emptySet;
	set_union_update(fNeighbors[c],fNeighbors[a],fNeighbors[b]); set_union_update(oFNeighbors[c],oFNeighbors[a],oFNeighbors[b]);
	fNeighbors[c].erase(a); fNeighbors[c].erase(b); oFNeighbors[c].erase(a); oFNeighbors[c].erase(b); 
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
	for (intit = oFNeighbors[c].begin(); intit != oFNeighbors[c].end(); ++intit) {
		x = (*intit);
		oFNeighbors[x].insert(c);
	}
	for (d = 0;d<this->dim;d++) {
		for (std::set<int>::iterator intit (fNeighbors[c].begin()) ; intit != fNeighbors[c].end(); ++intit) {
			x = (*intit);
			this->AddPair(c,x,this->MyNullPairStat->Add3(this->get_uv(a,x,d),this->get_uv(b,x,d)),d);
			this->AddPair(x,c,this->MyNullPairStat->Add3(this->get_uv(x,a,d),this->get_uv(x,b,d)),d);
	
			this->AddPairSoFar(c,x,this->MyNullPairStat->Add3(this->get_uvSoFar(a,x,d),this->get_uvSoFar(b,x,d)),d);
			this->AddPairSoFar(x,c,this->MyNullPairStat->Add3(this->get_uvSoFar(x,a,d),this->get_uvSoFar(x,b,d)),d);
		}
		for (std::set<int>::iterator intit (oFNeighbors[c].begin()) ; intit != oFNeighbors[c].end(); ++intit) {
			x = (*intit);
			this->AddPairOriginal(c,x,this->MyNullPairStat->Add3(this->get_uvOriginal(a,x,d),this->get_uvOriginal(b,x,d)),d);
			this->AddPairOriginal(x,c,this->MyNullPairStat->Add3(this->get_uvOriginal(x,a,d),this->get_uvOriginal(x,b,d)),d);	
		}
	}

}


bool dataMapOther::AddPairOriginal(int u, int v, ModelPairStatsBase* p,int d) {
	temp_outIt = oDatpair[d].find(u);
	if ( temp_outIt != oDatpair[d].end()) {
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
		oDatpair[d][u] = mnew;
	}
	return 1;
};


bool dataMapOther::AddPairSoFar(int u, int v, ModelPairStatsBase* p,int d) {
	temp_outIt = sDatpair[d].find(u);
	if ( temp_outIt != sDatpair[d].end()) {
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
		sDatpair[d][u] = mnew;
	}
	return 1;
};


bool dataMapOther::AddEdgeOriginal(int u, int v, float x,int d) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = oDatpair[d].find(u);
	if (temp_outIt == oDatpair[d].end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		oDatpair[d][u] = nd;
		temp_outIt = oDatpair[d].find(u);
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
	oDatvert[d][u]->AddOutGoing(x);
	return ans;
}

bool dataMapOther::AddEdgeSoFar(int u, int v, float x,int d) {
	// returns True if a new edge was added, false otherwise
	bool ans = 0;
	temp_outIt = sDatpair[d].find(u);
	if (temp_outIt == sDatpair[d].end()) {
		std::tr1::unordered_map<int,ModelPairStatsBase*> nd;
		sDatpair[d][u] = nd;
		temp_outIt = sDatpair[d].find(u);
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
	sDatvert[d][u]->AddOutGoing(x);
	return ans;
}





ModelPairStatsBase* dataMapOther::get_uvOriginal(int u, int v, int d) {
	//NOTE: returns zero if not present
	temp_outIt = oDatpair[d].find(u);
	if (temp_outIt!=oDatpair[d].end()) {
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


ModelPairStatsBase* dataMapOther::get_uvSoFar(int u, int v, int d) {
	//NOTE: returns zero if not present
	temp_outIt = sDatpair[d].find(u);
	if (temp_outIt != sDatpair[d].end()) {
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


void dataMapOther::initialize(graphData* D, graphData* Doriginal, graphData* DsoFar, int _dim) {
	std::set<int> emptySet;
	std::cout << "inside dataMapOther initialize\n";
	unsigned int u, v;
	this->dim = _dim;
	this->datpair = new std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> >[this->dim];
	this->datvert = new std::vector<ModelSelfStatsBase*>[this->dim];
	this->oDatpair = new std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> >[this->dim];
	this->oDatvert = new std::vector<ModelSelfStatsBase*>[this->dim];
	this->sDatpair = new std::tr1::unordered_map<int, std::tr1::unordered_map<int,ModelPairStatsBase*> >[this->dim];
	this->sDatvert = new std::vector<ModelSelfStatsBase*>[this->dim];
	this->gtype = new char[this->dim];
	for(int d=0;d<this->dim;d++) {
		assert(D[d].gtype==Doriginal[d].gtype);
		assert(D[d].gtype==DsoFar[d].gtype);
		this->gtype[d] = D[d].gtype;
	
		if (gtype[d]=='b') MyNullPairStat = new BinomialPairStats;
		else if (gtype[d]=='p') MyNullPairStat = new PoissonPairStats;
		else if (gtype[d]=='w') MyNullPairStat = new WPairStats;
		else if (gtype[d]=='d') MyNullPairStat = new DcorrPairStats;
		else if (gtype[d]=='g') MyNullPairStat = new GaussianPairStats;
		else {std::cout << "bad graph type\n";}
	
		std::set<int> emptySet;
		for (u=0; u != D[d].numV; u++) {
			this->initVert(u);
		}

		// use current graph to make current first neighbors
		std::map<int, graphData::destList*>::iterator it1 (D[d].edgeList.begin());
		for (; it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
				v = (*it2).first;
				assert(this->AddEdge(u,v,(*it2).second,d));
				if (fNeighbors.find(v)==fNeighbors.end()) fNeighbors[v] = emptySet;
				fNeighbors[u].insert(v);
				fNeighbors[v].insert(u);
			}
		}
	
		// Use the Original graph for making Original first neighbors
		for (it1 = Doriginal[d].edgeList.begin(); it1 != Doriginal[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			if (fNeighbors.find(u)==fNeighbors.end()) fNeighbors[u] = emptySet;
			for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
				v = (*it2).first;
				assert(this->AddEdgeOriginal(u,v,(*it2).second,d));
				if (oFNeighbors.find(v)==oFNeighbors.end()) oFNeighbors[v] = emptySet;
				oFNeighbors[u].insert(v);
				oFNeighbors[v].insert(u);
			}
		}
	
		for (it1 = DsoFar[d].edgeList.begin(); it1 != DsoFar[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			for (graphData::destList::iterator it2 ((*it1).second->begin()); it2 != (*it1).second->end(); ++it2) {
				v = (*it2).first;
				assert(this->AddEdgeSoFar(u,v,(*it2).second,d));
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



#endif
