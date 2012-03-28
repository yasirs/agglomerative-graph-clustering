#ifndef LINKPREDICTOR_HPP
#define LINKPREDICTOR_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "Engine.hpp"
#include <cassert>
#include <cmath>
#include <math.h>

class linkPredictor {
	public:
		std::map<int,std::map<int, ModelParamBase**> > topParams;
		int dim;
		dataMap* w;
		TreeClass* tree;
		graphData* D;
		bool attached;
		linkPredictor() {
			dim = 0;
			attached = 0;
		}
		virtual void attach(Engine* e);
		float predictEdge(unsigned int u, unsigned int v, int d);
		graphData* makeNonEdgePred(graphData* Dref);
		graphData* makeEdgePred(graphData* Dref);
		std::vector<float>* returnNonEdgePred(graphData* Dref, std::vector<std::pair<int,int> >& nonEdges);
		graphData* makeCompleteEdgePred();
		graphData* copyNoEdges(graphData* Dold);
		void addPredstoGraph(graphData* PD);
		~linkPredictor();
};

void linkPredictor::addPredstoGraph(graphData* PD) {
	// this adds predictions for existing edges in PD
	int d,u,v;
	float weight;
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	for (d=0; d<dim; d++) {
		for (dit = PD[d].edgeList.begin(); dit != PD[d].edgeList.end(); ++dit) {
			u = (*dit).first;
			for (eit = PD[d].edgeList[u]->begin(); eit != PD[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				weight = this->predictEdge(u,v,d);
				assert ( PD[d].Add_uv(u,v,weight));
			}
		}
	}
};


graphData* linkPredictor::copyNoEdges(graphData* Dold) {
	int d;
	graphData* Dnew = new graphData[dim];
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = Dold[d].int2Name;
		Dnew[d].name2Int = Dold[d].name2Int;
		Dnew[d].gtype = 'b'; // TODO:: come back to this: why 'b', and not Dold.gtype?
		Dnew[d].numV = Dold[d].numV;
	}
	return Dnew;
};


graphData* linkPredictor::makeEdgePred(graphData* Dref) {
	assert (attached);
	graphData* PD;
	float NP, weight;
	int d,u,v;
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	PD = new graphData[dim];
	for (d=0; d<dim;d++) {
		NP = 0;
		PD[d].gtype = 'w';
		PD[d].Etot = 0;
		PD[d].numV = D[d].numV;
		PD[d].int2Name = D[d].int2Name;
		PD[d].name2Int = D[d].name2Int;
		for (dit = Dref[d].edgeList.begin(); dit != Dref[d].edgeList.end(); ++dit) {
			u = (*dit).first;
			for (eit = (*dit).second->begin(); eit != (*dit).second->end(); eit++) {
				v = (*eit).first;
				weight = this->predictEdge(u,v,d);
				PD[d].Add_uv(u,v,weight);
				PD[d].Etot += weight;
				NP += 1;
			}
		}
	}
	return PD;
}

std::vector<float>* linkPredictor::returnNonEdgePred(graphData* Dref, std::vector<std::pair<int,int> >& nonEdges) {
	assert(attached);
	std::vector<float>* nePreds;
	nePreds = new std::vector<float>[dim];
	int u,v,d,i;
	float w;
	for (d=1;d<dim;d++) assert(D[d].numV == D[0].numV);
	for (d=0;d<dim;d++) {
		i = 0;
		for (u=0; u<D[d].numV; u++) {
			for (v=0;v<u; v++) {
				if (! Dref[d].has_uv(u,v)) {
					assert(nonEdges[i].first==u);
					assert(nonEdges[i].second==v);
					nePreds.append(this->predictEdge(u,v,d));
				}
			}
		}
	}
	return nePreds;
}


graphData* linkPredictor::makeCompleteEdgePred() {
	assert (attached);
	graphData* PD;
	float w, NP;
	int d;
	unsigned int u,v;
	PD = new graphData[dim];
	for (d=0; d<dim;d++) {
		NP = 0;
		PD[d].gtype = 'w';
		PD[d].Etot = 0;
		PD[d].numV = D[d].numV;
		PD[d].int2Name = D[d].int2Name;
		PD[d].name2Int = D[d].name2Int;
		for (u=0; u<D[d].numV; u++) {
			for (v=0; v<D[d].numV; v++) {
				if (u!=v) {
					w = this->predictEdge(u,v,d);
					assert(! PD[d].Add_uv(u,v,w));
					PD[d].Etot += w;
					NP += 1;
				}
			}
		}
		PD[d].aveP = NP/(D[d].numV * D[d].numV);
	}
	return PD;
};




graphData* linkPredictor::makeNonEdgePred(graphData* Dref) {
	assert (attached);
	graphData* PD;
	float w, NP;
	int d;
	unsigned u,v;
	PD = new graphData[dim];
	for (d=0; d<dim;d++) {
		NP = 0;
		PD[d].gtype = 'w';
		PD[d].Etot = 0;
		PD[d].numV = D[d].numV;
		PD[d].int2Name = D[d].int2Name;
		PD[d].name2Int = D[d].name2Int;
		for (u=0; u<D[d].numV; u++) {
			for (v=0; v<D[d].numV; v++) {
				if (u!=v) {
					if (! Dref[d].has_uv(u,v)) {
						w = this->predictEdge(u,v,d);
						assert(! PD[d].Add_uv(u,v,w));
						PD[d].Etot += w;
						NP += 1;
					}
				}
			}
		}
		PD[d].aveP = NP/(D[d].numV * D[d].numV);
	}
	return PD;
};






void linkPredictor::attach(Engine* e) {
	// delete old topParams, if any
	int d;
	if (attached) {
		std::map<int, std::map<int, ModelParamBase**> >::iterator outit;
		std::map<int, ModelParamBase**>::iterator init;
		for (outit = topParams.begin(); outit != topParams.end(); ++outit) {
			for (init = (*outit).second.begin(); init != (*outit).second.end(); ++init) {
				for (d=0;d<dim;d++) {
					delete init->second[d];
				}
				delete[] (*init).second;
			}
			topParams[(*outit).first].clear();
		}
		topParams.clear();
	}
	// assign new tree, w
	tree = e->tree;
	w = e->w;
	dim = e->dim;
	D = e->D;
	// compute new topThetas
	int n1, n2;
	std::set<int>::iterator intit1, intit2, intit3;
	for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); ++intit1) {
		n1 = (*intit1);
		topParams[n1] = std::map<int, ModelParamBase**>();
		for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); ++intit2) {
			n2 = (*intit2);
			if (n1 != n2) {
				topParams[n1][n2] = new ModelParamBase* [dim];
				for (d=0;d<dim;d++) {
					if (D[d].gtype=='w') {
						topParams[n1][n2][d] = new WParam;
					} else if (D[d].gtype=='b') {
						topParams[n1][n2][d] = new BinomialParam;
					} else if (D[d].gtype=='p') {
						topParams[n1][n2][d] = new PoissonParam;
					} else if (D[d].gtype=='d') {
						topParams[n1][n2][d] = new DcorrParam;
					} else if (D[d].gtype=='g') {
						topParams[n1][n2][d] = new GaussianParam;
					} else {
						std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Params).\n";
						throw 1;
					}
					
					topParams[n1][n2][d]->calculate(w->get_uv(n1,n2,d),w->datvert[d][n1],w->datvert[d][n2]);
					topParams[n1][n2][d]->cleanup();
				}
			}
		}
	}
	attached = 1;
};



float linkPredictor::predictEdge(unsigned int u, unsigned int v, int d) {
	if (!attached) {
		std::cerr << "can't do link prediction, not attached\n";
		throw 1;
	}
	if (u>=D[d].numV) {
		std::cerr << "Only "<<D[d].numV<<" vertices in the graph, can't predict edges from vertex "<<u<<"\n";
		throw 1;
	}
	if (v>=D[d].numV) {
		std::cerr << "Only "<<D[d].numV<<" vertices in the graph, can't predict edges to vertex "<<v<<"\n";
		throw 1;
	}
	std::pair<int,int> ipair = tree->getLCA(u,v);
	ModelParamBase* param;
	if (ipair.first==ipair.second) {
		param = tree->nodeMap[ipair.first]->params[d];
	} else {
		param = topParams[ipair.first][ipair.second][d];
	}
	return param->predict(w->datvert[d][u],w->datvert[d][v]);
};

linkPredictor::~linkPredictor() {
	std::map<int, std::map<int, ModelParamBase**> >::iterator outit;
	std::map<int, ModelParamBase**>::iterator init;
	for (outit = topParams.begin(); outit != topParams.end(); ++outit) {
		for (init = (*outit).second.begin(); init != (*outit).second.end(); ++init) {
			for (int d=0;d<dim;d++) {
				delete (*init).second[d];
			}
			delete[] (*init).second;
		}
		topParams[(*outit).first].clear();
	}
	topParams.clear();
}



#endif
