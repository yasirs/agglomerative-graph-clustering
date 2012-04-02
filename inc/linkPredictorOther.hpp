#ifndef LINKPREDICTOROTHER_HPP
#define LINKPREDICTOROTHER_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "Engine.hpp"
#include "mygenerators.hpp"
#include "linkPredictor.hpp"
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <math.h>

class linkPredictorOther: public linkPredictor {
	public:
		//std::map<int,std::map<int, float*> > topThetas;
		//int dim;
		dataMapOther* ww;
		TreeClassOther* ttree;
		//graphData* D;
		//bool attached;
		//linkPredictorOther() {
		//	attached = 0;
		//}
		virtual void attach(Engine* e);
		//float predictEdge(int u, int v, int d);
		//graphData* makeNonEdgePred(graphData* Dref);
		//graphData* makeEdgePred(graphData* Dref);
		//graphData* makeCompleteEdgePred();
		//graphData* copyNoEdges(graphData* Dold);
		//void addPredstoGraph(graphData* PD);
		void updateSoFarLazy(graphData* GsoFar);
		void updateTrackEdgeHoles(graphData* soFar);
		void integrateEdgeHoles(graphData* soFarEdgeHoles, graphData* soFarIntegrated);
		std::vector<float>* integrateEdgeHoleNonEdgePred(graphData *Dref, std::vector<std::pair<int, int> >& nonEdges, graphData* soFarEdgeHoles);

		virtual ~linkPredictorOther();
};


void linkPredictorOther::updateTrackEdgeHoles(graphData* soFar) {
	if (dim%2 != 0) {
		throw(std::runtime_error("To track edges and holes, need 2*d graphs"));
	}
	for (int d=0; d<dim/2; d++) {
		for (int u=0; u < soFar[d].numV; u++) {
			for (int v=0; v < soFar[d].numV; v++) {
				float e_hat = this->predictEdge(u,v,2*d);
				float h_hat = this->predictEdge(u,v,2*d+1);
				float del_e = std::max(ww->get_uvSimple(u,v,2*d) - e_hat,   0.0f);
				float del_h = std::max(ww->get_uvSimple(u,v,2*d+1) - h_hat, 0.0f);
				float e_new = soFar[2*d].get_uv(u,v) + e_hat + del_e - del_h;
				float h_new = soFar[2*d+1].get_uv(u,v) + h_hat + del_h - del_e;
				soFar[2*d].set_uv(u,v,e_new);
				soFar[2*d+1].set_uv(u,v,h_new);
			}
		}
	}
}


std::vector<float>* linkPredictorOther::integrateEdgeHoleNonEdgePred(graphData* Dref, std::vector<std::pair<int,int> >& nonEdges, graphData* soFarEdgeHoles) {
	if (dim%2 != 0) {
		throw(std::runtime_error("To track edges and holes, need 2*d graphs"));
	}
	std::vector<float>* nePreds;
	nePreds = new std::vector<float>[dim/2];
	int u,v,d,i;
	float w;
	for (d=1;d<dim/2;d++) assert(D[d].numV == D[0].numV);
	for (d=0;d<dim/2;d++) {
		i = 0;
		for (u=0; u<D[d].numV; u++) {
			for (v=0;v<u; v++) {
				if (! Dref[d].has_uv(u,v)) {
					assert(nonEdges[i].first==u);
					assert(nonEdges[i].second==v);
					nePreds[d].push_back(.5 +.5*(soFarEdgeHoles[2*d].get_uv(u,v) - soFarEdgeHoles[2*d+1].get_uv(u,v)));
					i++;
				}
			}
		}
	}
	return nePreds;
}



void linkPredictorOther::integrateEdgeHoles(graphData* soFarEdgeHoles, graphData* soFarIntegrated) {
	if (dim%2 != 0) {
		throw(std::runtime_error("To track edges and holes, need 2*d graphs"));
	}
	for (int d=0; d<dim/2; d++) {
		for (int u=0; u < soFarEdgeHoles[d].numV; u++) {
			for (int v=0; v < soFarEdgeHoles[d].numV; v++) {
				soFarIntegrated[d].set_uv(u,v, 0.5+0.5*(soFarEdgeHoles[2*d].get_uv(u,v) - soFarEdgeHoles[2*d+1].get_uv(u,v)));
			}
		}
	}
}
				


linkPredictorOther::~linkPredictorOther() {
	/*if (attached) {
		std::map<int, std::map<int, ModelParamBase**> >::iterator outit;
		std::map<int, ModelParamBase**>::iterator init;
		for (outit = topParams.begin(); outit != topParams.end(); ++outit) {
			for (init = (*outit).second.begin(); init != (*outit).second.end(); ++init) {
				for (int d=0;d<dim;d++) {
					delete init->second[d];
				}
				delete[] (*init).second;
			}
			topParams[(*outit).first].clear();
		}
		topParams.clear();
	}*/
	// do nothing, the base class function does it for us
}





void linkPredictorOther::attach(Engine* e) {
	// delete old topParams, if any
	int d;
	ModelParamBase **tp;
	if (attached) {
		//delete ww;
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
	ModelParamBase **oriP, **sofP;
		
	// assign new tree, w
	tree = e->tree;
	assert(tree->isOther());
	ttree = (TreeClassOther* ) tree;
	w = e->w;
	dim = e->dim;
	assert(w->DerivedType().compare("dataMapOther")==0);
	ww = (dataMapOther*) w;
	D = e->D;
	
	oriP = new ModelParamBase*[dim];
	sofP = new ModelParamBase*[dim];
	for (d=0;d<dim;d++) {
		if (D[d].gtype=='w') {
			oriP[d] = new WParam;
			sofP[d] = new WParam;
		} else if (D[d].gtype=='b') {
			oriP[d] = new BinomialParam;
			sofP[d] = new BinomialParam;
		} else if (D[d].gtype=='p') {
			oriP[d] = new PoissonParam;
			sofP[d] = new PoissonParam;
		} else if (D[d].gtype=='d') {
			oriP[d] = new DcorrParam;
			sofP[d] = new DcorrParam;
		} else if (D[d].gtype=='g') {
			oriP[d] = new GaussianParam;
			sofP[d] = new GaussianParam;
		} else {
			std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Params).\n";
			throw 1;
		}
	}
	// compute new topThetas
	int n1, n2;
	std::set<int>::iterator intit1, intit2, intit3;
	for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); ++intit1) {
		n1 = (*intit1);
		topParams[n1] = std::map<int, ModelParamBase**>();
		for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); ++intit2) {
			n2 = (*intit2);
			if (n1 != n2) {
				tp = new ModelParamBase*[dim];
				topParams[n1][n2] = tp; 
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
					topParams[n1][n2][d]->calculate(ww->get_uv(n1,n2,d),ww->datvert[d][n1],ww->datvert[d][n2]);
					topParams[n1][n2][d]->cleanup();					
				}
			}
		}
	}
	attached = 1;
	for (d=0;d<dim;d++) {
		delete oriP[d];
		delete sofP[d];
	}
	delete[] sofP;
	delete[] oriP;
};

void linkPredictorOther::updateSoFarLazy(graphData* GsoFar) {
	float wnew;
	std::map<int, std::map<int, ModelParamBase**> >::iterator outit;
	std::map<int, ModelParamBase**>::iterator init;
	for (int d=0; d<this->dim;d++) {
		for (outit = topParams.begin(); outit != topParams.end(); ++outit) {
			int n1 = (*outit).first;
			for (init = (*outit).second.begin(); init != (*outit).second.end(); ++init) {
				// need to process this
				int n2 = (*init).first;
				ModelParamBase* par = (*init).second[d];
				if (not par->isZero()) {
					AllChildVertGenerator G1(tree, n1); 
					int u;
					ModelSelfStatsBase* su;
					for (; (! G1.isDone()); ) {
						u = G1.goNext();
						su = ww->datvert[d][u];
						int v;
						ModelSelfStatsBase* sv;
						AllChildVertGenerator G2(tree, n2);
						for (; (! G2.isDone()); ) {
							v = G2.goNext();
							if (u != v) {
								sv = ww->datvert[d][v];
								//oldpred = par->predict(su,sv);
								wnew = par->updatedSoFar(su,sv,GsoFar[d].get_uv(u,v));
								//wnew = 1 - (1 - GsoFar[d].get_uv(u,v))*(1 - wpredicted);
								if (wnew>EPS) {
									GsoFar[d].set_uv(u,v,wnew);
								} else {
									GsoFar[d].delete_uv(u,v);
								}
							}
						}
					}
				}
			}
		}
		NCNZGenerator NodeGen(tree, d);
		int n;
		for (; (! NodeGen.isDone()); ) {
			n = NodeGen.goNext();
			Node* pnode = tree->nodeMap[n];
			assert(pnode->isOther());
			ModelParamBase* par = ((NodeOther*)pnode )->params[d];
			if (not par->isZero()) {
				if (pnode->collapsed) {
					// need to go over the vertices
					assert(pnode->vertsComputed);
					AllChildVertGenerator G1(tree, n); 
					int u;
					ModelSelfStatsBase* su;
					for (; (! G1.isDone()); ) {
						u = G1.goNext();
						su = ww->datvert[d][u];
						int v;
						ModelSelfStatsBase* sv;
						AllChildVertGenerator G2(tree, n);
						for (; (! G2.isDone()); ) {
							v = G2.goNext();
							if (u != v) {
								sv = ww->datvert[d][v];
								wnew = par->updatedSoFar(su,sv,GsoFar[d].get_uv(u,v));
								if (wnew>EPS) {
									GsoFar[d].set_uv(u,v,wnew);
								} else {
									GsoFar[d].delete_uv(u,v);
								}
							}
						}
					}
				} else {
					// need to go over child - vertices (pair-wise)
					for (std::set<int>::iterator cit1(pnode->childSet.begin()); cit1 != pnode->childSet.end(); cit1++) {
						for (std::set<int>::iterator cit2(pnode->childSet.begin()); cit2 != pnode->childSet.end(); cit2++) {
							if ((*cit1) != (*cit2)) {
								AllChildVertGenerator G1(tree, *cit1); 
								int u;
								ModelSelfStatsBase* su;
								for (; (! G1.isDone()); ) {
									u = G1.goNext();
									su = ww->datvert[d][u];
									int v;
									ModelSelfStatsBase* sv;
									AllChildVertGenerator G2(tree, *cit2);
									for (; (! G2.isDone()); ) {
										v = G2.goNext();
										if (u != v) {
											sv = ww->datvert[d][v];
											wnew = par->updatedSoFar(su,sv,GsoFar[d].get_uv(u,v));
											if (wnew>EPS) {
												GsoFar[d].set_uv(u,v,wnew);
											} else {
												GsoFar[d].delete_uv(u,v);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

#endif			

