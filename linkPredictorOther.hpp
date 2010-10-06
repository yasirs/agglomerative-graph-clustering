#ifndef LINKPREDICTOROTHER_HPP
#define LINKPREDICTOROTHER_HPP

#include "nodetree.hpp"
#include "graphData.hpp"
#include "Engine.hpp"
#include "mygenerators.hpp"
#include <cassert>
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
		void updateSoFar(graphData* GsoFar);
		void updateSoFarLazy(graphData* GsoFar);
		virtual ~linkPredictorOther();
};


linkPredictorOther::~linkPredictorOther() {
	// do nothing, the base class function does it for us
}




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

void linkPredictorOther::attach(Engine* e) {
	// delete old topParams, if any
	int d;
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
	ModelParamBase *oriP, *sofP;
	// assign new tree, w
	tree = e->tree;
	w = e->w;
	ww = (dataMapOther*) w;
	dim = e->dim;
	D = e->D;
	// compute new topThetas
	int n1, n2;
	std::set<int>::iterator intit1, intit2, intit3;
	for (d=0;d<dim;d++) {
		if (D[d].gtype=='w') {
			oriP = new WParam;
			sofP = new WParam;
		} else if (D[d].gtype=='b') {
			oriP = new BinomialParam;
			sofP = new BinomialParam;
		} else if (D[d].gtype=='p') {
			oriP = new PoissonParam;
			sofP = new PoissonParam;
		} else if (D[d].gtype=='d') {
			oriP = new DcorrParam;
			sofP = new DcorrParam;
		} else {
			std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Params).\n";
			throw 1;
		}
		
		for (intit1 = tree->topLevel.begin(); intit1 != tree->topLevel.end(); ++intit1) {
			n1 = (*intit1);
			topParams[n1] = std::map<int, ModelParamBase**>();
			for (intit2 = tree->topLevel.begin(); intit2 != tree->topLevel.end(); ++intit2) {
				n2 = (*intit2);
				if (n1 != n2) {
					topParams[n1][n2] = new ModelParamBase* [dim];
					if (D[d].gtype=='w') {
						topParams[n1][n2][d] = new WParam;
					} else if (D[d].gtype=='b') {
						topParams[n1][n2][d] = new BinomialParam;
					} else if (D[d].gtype=='p') {
						topParams[n1][n2][d] = new PoissonParam;
					} else if (D[d].gtype=='d') {
						topParams[n1][n2][d] = new DcorrParam;
					} else {
						std::cerr << "graph type "<<D[d].gtype<<" not yet supported for link prediction (top Params).\n";
						throw 1;
					}
					oriP->calculate(ww[d].get_uvOriginal(n1,n2),ww[d].oDatvert[n1],ww[d].oDatvert[n2]);
					oriP->cleanup();
					sofP->calculate(ww[d].get_uvSoFar(n1,n2),ww[d].sDatvert[n1],ww[d].sDatvert[n2]);
					sofP->cleanup();
					topParams[n1][n2][d]->bestfromSoFar(oriP,sofP);
					//topParams[n1][n2][d]->cleanup();
				}
			}
		}
		delete oriP;
		delete sofP;
	}
	attached = 1;
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
				AllChildVertGenerator G1(tree, n1); 
				int u;
				ModelSelfStatsBase* su;
				for (; (! G1.isDone()); ) {
					u = G1.goNext();
					su = w[d].datvert[u];
					int v;
					ModelSelfStatsBase* sv;
					AllChildVertGenerator G2(tree, n2);
					for (; (! G2.isDone()); ) {
						v = G2.goNext();
						if (u != v) {
							sv = w[d].datvert[v];
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
		NCNZGenerator NodeGen(tree, d);
		int n;
		for (; (! NodeGen.isDone()); ) {
			n = NodeGen.goNext();
			Node* pnode = tree->nodeMap[n];
			ModelParamBase* par = pnode->params[d];
			if (pnode->collapsed) {
				// need to go over the vertices
				assert(pnode->vertsComputed);
				AllChildVertGenerator G1(tree, n); 
				int u;
				ModelSelfStatsBase* su;
				for (; (! G1.isDone()); ) {
					u = G1.goNext();
					su = w[d].datvert[u];
					int v;
					ModelSelfStatsBase* sv;
					AllChildVertGenerator G2(tree, n);
					for (; (! G2.isDone()); ) {
						v = G2.goNext();
						if (u != v) {
							sv = w[d].datvert[v];
							//wpredicted = par->predict(su,sv);
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
								su = w[d].datvert[u];
								int v;
								ModelSelfStatsBase* sv;
								AllChildVertGenerator G2(tree, *cit2);
								for (; (! G2.isDone()); ) {
									v = G2.goNext();
									if (u != v) {
										sv = w[d].datvert[v];
										//wpredicted = par->predict(su,sv);
										//wnew = 1 - (1 - GsoFar[d].get_uv(u,v))*(1 - wpredicted);
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

#endif			

