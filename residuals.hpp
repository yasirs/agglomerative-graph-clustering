#include "nodetree.hpp"
#include "graphData.hpp"


graphData* getResidual(graphData* D, TreeClass* tree, dataMap* w, const int dim) {
	int d,i,u,v,n1,n2;
	std::stack<int> st;
	graphData *Dnew = new graphData[dim];
	std::set<int>::iterator intit1, intit2, intit3;
	std::pair<int,int> ipair;
	float theta, wthis, wpredicted;
	std::map<int,std::map<int, float*> > topThetas;
	graphData::destList::iterator eit;
	for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); intit1++) {
		n1 = (*intit1);
		topThetas[n1] = std::map<int, float*>();
		for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); intit2++) {
			n2 = (*intit2);
			topThetas[n1][n2] = new float[dim];
			for (d=0;d<dim;d++) {
				if (D[d].gtype=='w') {
					theta = (w[d].get_uv(n1,n1) + w[d].get_uv(n2,n2) + w[d].get_uv(n1,n2))/
						(w[d].selfMissing[n1] + w[d].selfMissing[n2] + (w[d].degrees[n1] * w[d].degrees[n2]) - w[d].get_uv(n1,n2));
					theta = theta/(theta+1.0f);
					topThetas[n1][n2][d] = theta;
				} else if (D[d].gtype=='b') {
					theta = (w[d].get_uv(n1,n1) + w[d].get_uv(n2,n2) + w[d].get_uv(n1,n2))/
						(w[d].nV[n1] + w[d].nV[n2]);
					topThetas[n1][n2][d] = theta;
				} else {
					std::cout << "graph type "<<D[d].gtype<<" not yet supported for residuals.\n";
					throw 1;
				}
			}
		}
	}
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = D[d].int2Name;
		Dnew[d].name2Int = D[d].name2Int;
		Dnew[d].gtype = 'w';
	}
	for (d=0;d<dim;d++) {
		for (u=0;u<D[0].numV; u++) {
			for (eit = D[d].edgeList[u]->begin(); eit != D[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wthis = (*eit).second;
				ipair = tree->getLCA(u,v);
				if (ipair.first==ipair.second) {
					theta = tree->nodeMap[ipair.first]->theta[d];
				} else {
					theta = topThetas[ipair.first][ipair.second][d];
				}
				wpredicted = w[d].degrees[u] * w[d].degrees[v] * theta;
				if ((wthis-wpredicted)>EPS) {
					Dnew[d].set_uv(u,v,wthis-wpredicted);
				}
			}
		}
	}
	return Dnew;
	
};
