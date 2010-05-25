#ifndef LINKPREDICTOR_HPP
#define LINKPREDICTOR_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "Engine.hpp"
#include <cassert>

class linkPredictor {
	std::map<int,std::map<int, float*> > topThetas;
	int dim;
	dataMap* w;
	TreeClass* tree;
	graphData* D;
	bool attached;
	linkPredictor() {attach = 0;}
	void attach(Engine* e);
	float predictEdge(int u, int v, int d);
	graphData* makeNonEdgePred();
	graphData* makeCompleteEdgePred();
	void addPredstoGraph(graphData* PD);
	
};

graphData* linkPredictor::makeNonEdgePred() {
	assert (attached);
	graphData* PD;
	float w, NP;
	PD = new graphData[dim];
	for (d=0; d<dim;d++) {
		NP = 0;
		PD[d].Etot = 0;
		PD[d].numV = D[d].numV;
		PD[d].int2Name = D[d].int2Name;
		PD[d].name2Int = D[d].name2Int;
		for (u=0; u<D[d].numV; u++) {
			for (v=0; v<D[d].numV; v++) {
				if (! D[d].has_uv(u,v)) {
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

void linkPredictor::addPredstoGraph(graphData* PD) {
	int d;
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	for (d=0; d<dim; d++) {
		for (dit = D[d].edgeList.begin(); dit != D[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = D[d].edgeList[u]->begin(); eit != D[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				w = this->predictEdge(u,v,d);
				assert ( PD[d].Add_uv(u,v,w));
			}
		}
	}
};


		
};



void linkPredictor::attach(Engine* e) {
	// delete old topThetas, if any
	std::map<int, std::map<int, float*> >::iterator outit;
	std::map<int,float*>::iterator init;
	for (outit = topThetas.begin(); outit != topThetas.end(); outit++) {
		for (init = (*outit).second.begin(); init != (*outit).second.end(); ++init) {
			delete[] (*init).second;
		}
	}
	// assign new tree, w
	tree = e->tree;
	w = e->w;
	dim = e->dim;
	D = e->D;
	// compute new topThetas
	int n1, n2, d;
	float theta;
	std::set<int>::iterator intit1, intit2, intit3;
	for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); intit1++) {
		n1 = (*intit1);
		topThetas[n1] = std::map<int, float*>();
		for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); intit2++) {
			n2 = (*intit2);
			topThetas[n1][n2] = new float[dim];
			for (d=0;d<dim;d++) {
				if (D[d].gtype=='w') {
					theta = w[d].get_uv(n1,n2)/(w[d].degrees[n1] * w[d].degrees[n2]);
					topThetas[n1][n2][d] = theta;
				} else if (D[d].gtype=='b') {
					theta = w[d].get_uv(n1,n2)/(w[d].nV[n1] * w[d].nV[n2]);
					topThetas[n1][n2][d] = theta;
				} else {
					std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Thetas).\n";
					throw 1;
				}
			}
		}
	}
	attached = 1;
};

float linkPredictor::predictEdge(int u, int v, int d) {
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
	float wpredicted, theta;
	if (ipair.first==ipair.second) {
		theta = tree->nodeMap[ipair.first]->theta[d];
	} else {
		theta = topThetas[ipair.first][ipair.second][d];
	}
	if (D[d].gtype=='w') {
		wpredicted = w[d].degrees[u] * w[d].degrees[v] * theta;
	} else if (D[d].gtype=='b') {
		wpredicted = w[d].nV[u] * w[d].nV[v] * theta;
	} else {
		std::cerr << "graph type "<<D[d].gtype<<" not supported yet for predicting edges.\n";
		throw 1;
	}
	return wpredicted;
};


#endif
