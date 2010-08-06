#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "linkPredictor.hpp"


/*void updateSoFar(graphData* Doriginal, linkPredictor& lp, graphData* GsoFar) {
	float NE = 0;
	int u,v,d;
	float wpredicted, w_old;
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	for (d=0;d<lp.dim;d++) {
		NE = 0;
		for (dit = Doriginal[d].edgeList.begin(); dit != Doriginal[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = Doriginal[d].edgeList[u]->begin(); eit != Doriginal[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wpredicted = lp.predictEdge(u,v,d);
				if (wpredicted > EPS) {
					w_old = GsoFar[d].get_uv(u,v);
					wpredicted = 1.0 - (1.0 - w_old)*(1.0 - wpredicted);
					GsoFar[d].set_uv(u, v, wpredicted);
					GsoFar[d].Etot += wpredicted;
					GsoFar[d].Etot -= w_old;
					NE += 1;
				}
			}
		}
		GsoFar[d].aveP = NE/(GsoFar[d].numV * (GsoFar[d].numV-1));
	}
};
*/
			
graphData* residualDiff(graphData* Doriginal, graphData* GsoFar, int dim) {
	graphData *Dnew = new graphData[dim];
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	float NE = 0;
	int u,v,d;
	float wpredicted, wthis, thisweight;
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = Doriginal[d].int2Name;
		Dnew[d].name2Int = Doriginal[d].name2Int;
		Dnew[d].gtype = 'b';
		Dnew[d].numV = Doriginal[d].numV;
	}
	for (d=0;d<dim;d++) {
		Dnew[d].Etot = 0;
		NE = 0;
		for (dit = Doriginal[d].edgeList.begin(); dit != Doriginal[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = Doriginal[d].edgeList[u]->begin(); eit != Doriginal[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wthis = (*eit).second;
				wpredicted = GsoFar[d].get_uv(u,v);
				thisweight = wthis - wpredicted;
				if ((thisweight)>EPS) {
					Dnew[d].set_uv(u,v,thisweight);
					Dnew[d].Etot += thisweight;
					NE += 1;
				}
			}
		}
		Dnew[d].aveP = NE/(Dnew[d].numV * (Dnew[d].numV-1));
	}
	return Dnew;
};



/*
graphData* getORResidual(graphData* Doriginal, linkPredictor& lp, graphData* GsoFar) {
	// for now, this is necessary
	assert(Doriginal->gtype=='b');
	graphData *Dnew = new graphData[lp.dim];
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	float NE = 0;
	int u,v,d;
	float wpredicted, wthis, thisweight, w_old;
	for (d=0;d<lp.dim;d++) {
		Dnew[d].int2Name = Doriginal[d].int2Name;
		Dnew[d].name2Int = Doriginal[d].name2Int;
		Dnew[d].gtype = 'b';
		Dnew[d].numV = Doriginal[d].numV;
	}
	for (d=0;d<lp.dim;d++) {
		Dnew[d].Etot = 0;
		NE = 0;
		for (dit = Doriginal[d].edgeList.begin(); dit != Doriginal[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = Doriginal[d].edgeList[u]->begin(); eit != Doriginal[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wthis = (*eit).second;
				wpredicted = lp.predictEdge(u,v,d);
				// correct wpredicted
				w_old = GsoFar[d].get_uv(u,v);
				wpredicted = 1.0 - (1.0 - w_old)*(1.0 - wpredicted);
				if ((wthis-wpredicted)>EPS) {
					thisweight = (wthis-wpredicted)/(1-wpredicted);
					Dnew[d].set_uv(u,v,thisweight);
					Dnew[d].Etot += thisweight;
					NE += 1;
				}
			}
		}
		Dnew[d].aveP = NE/(Dnew[d].numV * (Dnew[d].numV-1));
	}
	return Dnew;
};
*/	


graphData* getResidual(graphData* D, linkPredictor& lp) {
	graphData *Dnew = new graphData[lp.dim];
	graphData::destList::iterator eit;
	std::map<int, graphData::destList*>::iterator dit;
	float NE = 0;
	int u,v,d;
	float wpredicted, wthis;
	for (d=0;d<lp.dim;d++) {
		Dnew[d].int2Name = D[d].int2Name;
		Dnew[d].name2Int = D[d].name2Int;
		Dnew[d].gtype = 'w';
		Dnew[d].numV = D[d].numV;
	}
	for (d=0;d<lp.dim;d++) {
		Dnew[d].Etot = 0;
		NE = 0;
		for (dit = D[d].edgeList.begin(); dit != D[d].edgeList.end(); dit++) {
			u = (*dit).first;
			for (eit = D[d].edgeList[u]->begin(); eit != D[d].edgeList[u]->end(); eit++) {
				v = (*eit).first;
				wthis = (*eit).second;
				wpredicted = lp.predictEdge(u,v,d);
				if ((wthis-wpredicted)>EPS) {
					Dnew[d].set_uv(u,v,wthis-wpredicted);
					Dnew[d].Etot += wthis-wpredicted;
					NE += 1;
				}
			}
		}
		Dnew[d].aveP = NE/(Dnew[d].numV * (Dnew[d].numV-1));
	}
	return Dnew;
};


#endif
