#ifndef LINKPREDICTOROTHER_HPP
#define LINKPREDICTOROTHER_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "Engine.hpp"
#include <cassert>
#include <cmath>
#include <math.h>

class linkPredictorOther {
	public:
		std::map<int,std::map<int, float*> > topThetas;
		int dim;
		dataMapOther* w;
		TreeClassOther* tree;
		graphData* D;
		bool attached;
		linkPredictorOther() {
			attached = 0;
		}
		void attach(Engine* e);
		float predictEdge(int u, int v, int d);
		graphData* makeNonEdgePred(graphData* Dref);
		graphData* makeEdgePred(graphData* Dref);
		graphData* makeCompleteEdgePred();
		graphData* copyNoEdges(graphData* Dold);
		void addPredstoGraph(graphData* PD);
		void updateSoFar(graphData* GsoFar);
		~linkPredictorOther();
};


graphData* linkPredictorOther::copyNoEdges(graphData* Dold) {
	int d;
	graphData* Dnew = new graphData[dim];
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = Dold[d].int2Name;
		Dnew[d].name2Int = Dold[d].name2Int;
		Dnew[d].gtype = 'b';
		Dnew[d].numV = Dold[d].numV;
	}
	return Dnew;
};


graphData* linkPredictorOther::makeEdgePred(graphData* Dref) {
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

graphData* linkPredictorOther::makeNonEdgePred(graphData* Dref) {
	assert (attached);
	graphData* PD;
	float w, NP;
	int d,u,v;
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

void linkPredictorOther::addPredstoGraph(graphData* PD) {
	// this add predictions for existing edges in PD
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


void linkPredictorOther::updateSoFar(graphData* GsoFar) {
	// TODO:: lazy V^2 computation right now, need to update it to go over only the nonzero thetas
	float wpredicted, wnew;
	for (int d=0; d< this->dim;d++) {
		for (int u= 0; u<GsoFar[d].numV; u++) {
			for (int v= 0; v<GsoFar[d].numV; v++) {
				wpredicted = this->predictEdge(u,v,d);
				wnew = 1 - (1 - GsoFar[d].get_uv(u,v))*(1 - wpredicted);
				if (wnew>EPS) {
					GsoFar[d].set_uv(u,v,wnew);
				} else {
					GsoFar[d].delete_uv(u,v);
				}
			}
		}
	}
}


linkPredictorOther::~linkPredictorOther() {
	std::map<int, std::map<int, float*> >::iterator outit;
	std::map<int,float*>::iterator init;
	for (outit = topThetas.begin(); outit != topThetas.end(); ++outit) {
		for (init = (*outit).second.begin(); 
			init != (*outit).second.end();
			++init) {
			delete[] (*init).second;
		}
		topThetas[(*outit).first].clear();
	}
	topThetas.clear();
}



void linkPredictorOther::attach(Engine* e) {
	// delete old topThetas, if any
	std::map<int, std::map<int, float*> >::iterator outit;
	std::map<int,float*>::iterator init;
	for (outit = topThetas.begin(); outit != topThetas.end(); ++outit) {
		for (init = (*outit).second.begin(); 
			init != (*outit).second.end();
			++init) {
			delete[] (*init).second;
		}
		topThetas[(*outit).first].clear();
	}
	topThetas.clear();

	// assign new tree, w
	tree = (TreeClassOther*) e->tree;
	w = (dataMapOther*) e->w;
	dim = e->dim;
	D = e->D;
	// compute new topThetas
	int n1, n2, d;
	float thSoFar, thOriginal;
	std::set<int>::iterator intit1, intit2, intit3;
	for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); ++intit1) {
		n1 = (*intit1);
		topThetas[n1] = std::map<int, float*>();
		for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); ++intit2) {
			n2 = (*intit2);
			topThetas[n1][n2] = new float[dim];
			for (d=0;d<dim;d++) {
				if (D[d].gtype=='w') {
					thSoFar = w[d].get_uvSoFar(n1,n2) / (w[d].sDegrees[n1] * w[d].sDegrees[n2]);
					if (std::isnan(thSoFar)) thSoFar=0;
					thOriginal = w[d].get_uvOriginal(n1,n2) / (w[d].oDegrees[n1] * w[d].oDegrees[n2]);
					if (std::isnan(thOriginal)) thOriginal=0;
					topThetas[n1][n2][d] = std::max(0.0f,(thOriginal - thSoFar)/(1 - thSoFar));
				} else if (D[d].gtype=='b') {
					thOriginal = w[d].get_uvOriginal(n1,n2)/(w[d].oNV[n1] * w[d].oNV[n2]);
					if (std::isnan(thOriginal)) thOriginal=0;
					thSoFar = w[d].get_uvSoFar(n1,n2)/(w[d].oNV[n1] * w[d].oNV[n2]);
					if (std::isnan(thSoFar)) thSoFar=0;
					topThetas[n1][n2][d] = std::max(0.0f,(thOriginal - thSoFar)/(1 - thSoFar));
				} else {
					std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Thetas).\n";
					throw 1;
				}
			}
		}
	}
	attached = 1;


}






float linkPredictorOther::predictEdge(int u, int v, int d) {
	float ma, mb;
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
		ma = w[d].degrees[u];
		mb = w[d].degrees[v];
	} else if (D[d].gtype=='b') {
		ma = w[d].nV[u];
		mb = w[d].nV[v];
	} else {
		std::cerr << "graph type "<<D[d].gtype<<" not supported yet for predicting edges.\n";
		throw 1;
	}
	wpredicted = ma * theta * mb;
	if ((wpredicted>1)||(wpredicted<0)) {
		std::cerr << "bad weight predicted?\n";
		// TODO:: erase this! only for debug
	}
	return wpredicted;
};


#endif
