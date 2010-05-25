#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

#include "nodetree.hpp"
#include "graphData.hpp"


graphData* getResidual(graphData* D, TreeClass* tree, dataMap* w, const int dim) {
	int d,u,v,n1,n2;
	std::stack<int> st;
	graphData *Dnew = new graphData[dim];
	std::set<int>::iterator intit1, intit2, intit3;
	std::pair<int,int> ipair;
	float theta, wthis, wpredicted, NE;
	std::map<int,std::map<int, float*> > topThetas;
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
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
					std::cerr << "graph type "<<D[d].gtype<<" not yet supported for residuals.\n";
					throw 1;
				}
			}
		}
	}
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = D[d].int2Name;
		Dnew[d].name2Int = D[d].name2Int;
		Dnew[d].gtype = 'w';
		Dnew[d].numV = D[d].numV;
	}
	for (d=0;d<dim;d++) {
		Dnew[d].Etot = 0;
		NE = 0;
		for (dit = D[d].edgeList.begin(); dit != D[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = D[d].edgeList[u]->begin(); eit != D[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wthis = (*eit).second;
				ipair = tree->getLCA(u,v);
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
					std::cerr << "graph type "<< D[d].gtype <<" not supported yet for computing residuals\n";
					throw 1;
				}
				if ((wthis-wpredicted)>EPS) {
					Dnew[d].set_uv(u,v,wthis-wpredicted);
					Dnew[d].Etot += wthis-wpredicted;
					NE += 1;
				}
			}
		}
		Dnew[d].aveP = NE/(Dnew[d].numV * Dnew[d].numV);
	}
	return Dnew;
};

#endif
