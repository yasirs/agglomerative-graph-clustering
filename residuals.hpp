#include "nodetree.hpp"
#include "graphData.hpp"


graphData* getResidual(graphData* D, TreeClass* tree, dataMap* w, const int dim) {
	int d;
	int t,p,n,i,u,v;
	std::stack<int> st;
	graphData *Dnew = new graphData[dim];
	std::set<int>::iterator intit, intit2, intit3;
	std::pair<int,int> ipair;
	float theta;
	std::map<int,std::map<int, float*> > topThetas;
	graphData::destList::iterator eit;
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = D[d].int2Name;
		Dnew[d].name2Int = D[d].name2Int;
	}
	for (d=0;d<dim;d++) {
		for (u=0;u<D[0].numV; i++) {
			for (eit = edgeList[u]->begin(); eit = edgeList[u]->end(); eit++) {
				v = (*eit).first;
				w = (*eit).second;
				ipair tree->getLCA(u,v);
				if (ipair.first==ipair.second) theta = tree->nodeMap[ipair.first]->theta[d];
				else theta = topThetas[ipair.first][ipair.second][d];
				wpredicted = w[d].degrees[u] * w[d].degrees[v] * theta;
				if ((w-wpredicted)>EPS) Dnew[d].set_uv(u,v,w-wpredicted);
			}
		}
	}
	
};
