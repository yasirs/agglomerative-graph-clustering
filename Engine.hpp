#ifndef ENGINE_HPP
#define ENGINE_HPP
#include <map>
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"
#include "scoremap.hpp"
#include "dataMap.hpp"
#include "mysetfuncs.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>



#if (! NOGSL)
#include "gsl/gsl_sf.h"
#else
double lgamma(double x) {
	double ans;
	ans = x*log(x)-x;
	return ans;
};
double gsl_sf_lnbeta(double a, double b) {
	double ans;
	ans = lgamma(a)+lgamma(b) -lgamma(a+b);
	return ans;
};
double gsl_sf_gamma(double a) {
	return exp(lgamma(a));
};

#endif

#ifndef DEBUGMODE
#define DEBUGMODE 0
#endif


#ifndef NOREFERENCE
#define NOREFERENCE 0
#endif

class Engine{
	public:
		graphData *D;
		int dim;
		int curLev;
		dataMap* w;
		scoremap sm;
		int* NodeMembership;
		TreeClass *tree;
		std::map<int,std::set<int> > firstNeighbors;
		std::map<int,std::set<int> > secondNeighbors;
		//std::map<int, float> groupDegrees;
		~Engine();
		Engine(graphData* G, int d);
		bool initializeFirstLev();
		int run();
		bool doStage();
		bool cleanUp();
		double deltascore(int d, int a, int b, int x);
		double centerscore(int d, int a, int b);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
		void printJaccardFile(const char* filename, int d, bool edges);
		void printHyperGeomFile(const char* filename, int d, bool edges);
		void printDegreeProdFile(const char* filename, int d, bool edges);
		void printCommonNeighbFile(const char* filename, int d, bool edges);
};


void Engine::printCommonNeighbFile(const char* fn, int d, bool edges) {
	int u,v, c;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		for (v=0; v<u; v++) {
			if ((edges)or(not D[d].has_uv(u,v))) {
				c = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
				file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << c << '\n';
			}
		}
	}
	file.close();
};


void Engine::printJaccardFile(const char* fn, int d, bool edges) {
	int u,v, aub, aib, ad;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		ad = D[d].degree(u);
		for (v=0; v<u; v++) {
			if ((edges)or(not D[d].has_uv(u,v))) {
				aib = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
				aub = (ad * D[d].degree(v)) - aib;
				file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << (aib+0.0f)/aub << '\n';
			}
		}
	}
	file.close();
};


void Engine::printHyperGeomFile(const char* fn, int d, bool edges) {
	int u, v, m, n, c, t, x, dmin;
	double s;
	std::ofstream file;
	file.open(fn,std::ios::out);
	t = D[d].numV-2;
	for (u=0; u<D[d].numV; u++) {
		m = D[d].degree(u);
		for (v=0; v<u; v++) {
			if ((edges)or(not D[d].has_uv(u,v))) {
				m = D[d].degree(v);
				dmin = std::min(m,n);
				c = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
				s = 0;
				for (x = c; x<= dmin; x++) {
					s = s + 1.0 / (gsl_sf_gamma(x)*gsl_sf_gamma(m-x)*gsl_sf_gamma(n-x)*gsl_sf_gamma(t-m-n+x));
				}
				s = s * gsl_sf_gamma(m) * gsl_sf_gamma(n) * gsl_sf_gamma(t-m) * gsl_sf_gamma(t-n) / gsl_sf_gamma(t);
				s = log10(s);
				file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << s << '\n';
			}
		}
	}
	file.close();
};


void Engine::printDegreeProdFile(const char* fn, int d, bool edges) {
	int u,v, ad, bd;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		ad = D[d].degree(u);
		for (v=0; v<u; v++) {
			if ((edges)or(not D[d].has_uv(u,v))) {
				bd = D[d].degree(v);
				file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << ad*bd << '\n';
			}
		}
	}
	file.close();
};













Engine::Engine(graphData* G, int d) {
	D = G;
	dim = d;
	w = new dataMap[dim];
	tree = new TreeClass(G);
};


Engine::~Engine() {
	delete[] w;
};


double Engine::centerscore(int d, int a, int b) {
	assert(a!=b);
	double ans;
	if (D[d].gtype=='b') {
		float Tab, Eab, Eaa, Ebb, Haa, Hbb, Taa, Tbb, Tcc;
		Tab = w[d].nV[a] * w[d].nV[b];
		Taa = w[d].nV[a] * (w[d].nV[a] - 1)/2.0f;
		Tbb = w[d].nV[b] * (w[d].nV[b] - 1)/2.0f;
		Tcc = Taa + Tbb + Tab;
		Eab = w[d].get_uv(a,b);
		Eaa = w[d].get_uv(a,a);
		Ebb = w[d].get_uv(b,b);
		Haa = Taa - Eaa;
		Hbb = Tbb - Ebb;
		ans =   gsl_sf_lnbeta(Eaa + Eab + Ebb + 1,Haa + Hbb + Tab - Eab + 1)
			- gsl_sf_lnbeta(Eaa + 1, Haa +1)
			- gsl_sf_lnbeta(Ebb + 1, Hbb +1)
			- gsl_sf_lnbeta(Eab+1,Tab - Eab +1);
		//for the random reference model
#if (! NOREFERENCE)
		ans -=  gsl_sf_lnbeta(Tcc*D[d].aveP + 1, Tcc*(1 - D[d].aveP) + 1)
			- gsl_sf_lnbeta(Taa*D[d].aveP + 1, Taa*(1 - D[d].aveP) + 1)
			- gsl_sf_lnbeta(Tbb*D[d].aveP + 1, Tbb*(1 - D[d].aveP) + 1)
			- gsl_sf_lnbeta(Tab*D[d].aveP + 1, Tab*(1 - D[d].aveP) + 1);
#endif
	} else if (D[d].gtype=='w') {
		float Tab, Eab, Eaa, Ebb, Haa, Hbb, Taa, Tbb, Tcc;
		Eab = w[d].get_uv(a,b);
		Eaa = w[d].get_uv(a,a);
		Ebb = w[d].get_uv(b,b);
		Haa = w[d].selfMissing[a];
		Hbb = w[d].selfMissing[a];
		Tab = w[d].degrees[a] * w[d].degrees[b];
		Taa = Eaa + Haa;
		Tbb = Ebb + Hbb;
		Tcc = Taa + Tbb + Tab;
		ans =   gsl_sf_lnbeta(Eaa + Eab + Ebb + 1,Haa + Hbb + Tab - Eab + 1)
			- gsl_sf_lnbeta(Eaa + 1, Haa +1)
			- gsl_sf_lnbeta(Ebb + 1, Hbb +1)
			- gsl_sf_lnbeta(Eab+1,Tab - Eab +1);
		//for the random reference model
#if (! NOREFERENCE)
		float Ebaraa, Ebarbb, Ebarab, Ebarcc;
		Ebarab = w[d].degrees[a] * w[d].degrees[b];
		Ebaraa = Taa/(2.0f*D[d].Etot);
		Ebarbb = Tbb/(2.0f*D[d].Etot);
		Ebarcc = Ebaraa + Ebarbb + Ebarbb;
		ans -=  gsl_sf_lnbeta(Ebarcc + 1, Tcc - Ebarcc + 1)
			- gsl_sf_lnbeta(Ebaraa + 1, Taa - Ebaraa + 1)
			- gsl_sf_lnbeta(Ebarbb + 1, Tbb - Ebarbb + 1)
			- gsl_sf_lnbeta(Ebarab + 1, Tab - Ebarab + 1);	
#endif
	}
	return ans;
};


double Engine::deltascore(int d, int a, int b, int x) {
	assert((a!=b)&&(a!=x)&&(b!=x));
	double ans;
	float da, db, dx;
	if (D[d].gtype=='b') {
		float Eax, Ebx, Hax, Hbx, Tax, Tbx;
		float Ebarax, Ebarbx, Hbarax, Hbarbx;
		Eax = w[d].get_uv(a,x);
		Ebx = w[d].get_uv(b,x);
		da = w[d].nV[a];
		db = w[d].nV[b];
		dx = w[d].nV[x];
		Tax = da * dx;
		Tbx = db * dx;
		Hax = Tax - Eax;
		Hbx = Tbx - Ebx;
		Ebarax = Tax * D[d].aveP;
		Ebarbx = Tbx * D[d].aveP;
		Hbarax = Tax - Ebarax;
		Hbarbx = Tbx - Ebarbx;
		ans =    (	gsl_sf_lnbeta(Eax+Ebx+1.0f,Hax+Hbx+1.0f)
				-gsl_sf_lnbeta(Eax+1.0f,Hax+1.0f)
				-gsl_sf_lnbeta(Ebx+1.0f,Hbx+1.0f)
			 );
#if (! NOREFERENCE)
		ans -=   (	gsl_sf_lnbeta(Ebarax+Ebarbx+1.0f,Hbarax+Hbarbx+1.0f)
				-gsl_sf_lnbeta(Ebarax+1.0f,Hbarax+1.0f)
				-gsl_sf_lnbeta(Ebarbx+1.0f,Hbarbx+1.0f)
			 );
#endif
	} else if (D[d].gtype=='w') {
		float Eax, Ebx, Hax, Hbx;
		float Ebarax, Ebarbx, Hbarax, Hbarbx;
		Eax = w[d].get_uv(a,x);
		Ebx = w[d].get_uv(b,x);
		da = w[d].degrees[a];
		db = w[d].degrees[b];
		dx = w[d].degrees[x];
		Hax = da * dx - Eax;
		Hbx = db * dx - Ebx;
		Ebarax = da * dx/(2.0f*D[d].Etot);
		Ebarbx = db * dx/(2.0f*D[d].Etot);
		Hbarax = da * dx - Ebarax;
		Hbarbx = db * dx - Ebarbx;
		ans =    (	gsl_sf_lnbeta(Eax+Ebx+1.0f,Hax+Hbx+1.0f)
				-gsl_sf_lnbeta(Eax+1.0f,Hax+1.0f)
				-gsl_sf_lnbeta(Ebx+1.0f,Hbx+1.0f)
			 );
#if (! NOREFERENCE)
		ans -=	 (	gsl_sf_lnbeta(Ebarax+Ebarbx+1.0f,Hbarax+Hbarbx+1.0f)
				-gsl_sf_lnbeta(Ebarax+1.0f,Hbarax+1.0f)
				-gsl_sf_lnbeta(Ebarbx+1.0f,Hbarbx+1.0f)
			 );
#endif
	}
	return ans;
};


int Engine::run() {
	std::set<int> emptySet, tempSet;
	std::set<int> possSet1, possSet2;
	std::set<int>::iterator intit, intit2, intit3;
	int numJoins = 0;
	int a,b,c,x,y,z,d;
	float wc;
	double cscore, jscore;
	Node* pnode;
	scoremap::pairScore pscore;
	std::map<int, scoremap::smap>::iterator smOut;
	std::map<int, scoremap::twoScores>::iterator smIn;
	std::map<int, std::map<int,float> >::iterator datOut;
	std::map<int, float>::iterator datIn;
	while (sm.hasPos()) {
		// join and do stuff
		pscore = sm.popBestScore();
		a = pscore.u; b = pscore.v;
		sm.erase(b,a);
		// let us create new group c and create heirarchical relations
		c = tree->numNodes;
		c++;
		tree->numNodes = c;
		pnode = new Node(c,-1,0);
		pnode->theta = new float[dim];
		tree->nodeMap[c] = pnode;
		tree->nodeMap[a]->parent = c;
		tree->nodeMap[b]->parent = c;
		tree->nodeMap[c]->childSet.insert(a);
		tree->nodeMap[c]->childSet.insert(b);
		tree->topLevel.erase(a);
		tree->topLevel.erase(b);
		tree->topLevel.insert(c);
		if ((tree->nodeMap[a]->collapsed)&&(tree->nodeMap[b]->collapsed)&&(pscore.s.centerMscore>0)) {
			pnode->collapsed = 1;
			for(intit= tree->nodeMap[a]->vertexSet.begin(); intit!= tree->nodeMap[a]->vertexSet.end(); ++intit) {
				pnode->vertexSet.insert(*intit);
			}
			for(intit= tree->nodeMap[b]->vertexSet.begin(); intit!= tree->nodeMap[b]->vertexSet.end(); ++intit) {
				pnode->vertexSet.insert(*intit);
			}
			pnode->vertsComputed = 1;
		} else {
			pnode->collapsed = 0;
			pnode->vertsComputed = 0;
		}
		// compute d,w,n,m for c
		for (d=0;d<dim;d++) {
			wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
			assert(w[d].AddPair(c,c,wc));
			w[d].degrees[c] = w[d].degrees[a] + w[d].degrees[b];
			w[d].selfMissing[c] = w[d].selfMissing[a] + w[d].selfMissing[b] + (w[d].degrees[a] * w[d].degrees[b]) - w[d].get_uv(a,b);
			w[d].nV[c] = w[d].nV[a] + w[d].nV[b];
			if (D[d].gtype=='w') {
				tree->nodeMap[c]->theta[d] = w[d].get_uv(a,b)/(w[d].degrees[a] * w[d].degrees[b]);
			} else if (D[d].gtype=='b') {
				tree->nodeMap[c]->theta[d] = w[d].get_uv(a,b)/(w[d].nV[a] * w[d].nV[b]);
			} else {
				std::cerr << "graph type "<<D[d].gtype<<" not yet supported (during theta calculation).\n";
				throw 1;
			}

		}
		firstNeighbors[c] = emptySet;
		secondNeighbors[c] = emptySet;
		set_union_update(firstNeighbors[c],firstNeighbors[a],firstNeighbors[b]);
		firstNeighbors[c].erase(a); firstNeighbors[c].erase(b); 
		tempSet = emptySet;
		set_union_update(tempSet,secondNeighbors[a],secondNeighbors[b]);
		set_difference_update(secondNeighbors[c],tempSet,firstNeighbors[c]);
		secondNeighbors[c].erase(a); secondNeighbors[c].erase(b);
		for (d=0;d<dim;d++) {
			for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
				x = (*intit);
				w[d].AddPair(c,x,w[d].get_uv(a,x) + w[d].get_uv(b,x));
				w[d].AddPair(x,c,w[d].get_uv(x,a) + w[d].get_uv(b,x));
			}
		}
		// delete a,b scores
		for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
			x = *intit;
			if (b!=x) {
				sm.erase(a,x);
				sm.erase(x,a);
			}
		}
		for (intit = secondNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
			x = *intit;
			if (b!=x) {
				sm.erase(a,x);
				sm.erase(x,a);
			}
		}
		for (intit = firstNeighbors[b].begin(); intit != firstNeighbors[b].end(); ++intit) {
			x = *intit;
			if (a!=x) {
				sm.erase(b,x);
				sm.erase(x,b);
			}
		}
		for (intit = secondNeighbors[b].begin(); intit != secondNeighbors[b].end(); ++intit) {
			x = *intit;
			if (a!=x) {
				sm.erase(b,x);
				sm.erase(x,b);
			}
		}

		// update all existing old scores
		for (smOut = sm.scores.begin(); smOut != sm.scores.end(); ++smOut) {
			x = (*smOut).first;
			for (smIn = (*smOut).second.scoreDest.begin(); smIn != (*smOut).second.scoreDest.end(); ++smIn) {
				y = (*smIn).first;
				jscore =0;
				for (d=0;d<dim;d++) {
					jscore = jscore + deltascore(d,x,y,c)-deltascore(d,x,y,a)-deltascore(d,x,y,b);
				}
				sm.AddTo(x,y,jscore);
			}
		}
		
		//TODO:: debug check
		assert(sm.has_u(c)<0);
		// create new c,x scores
		for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
			x = *intit;

			firstNeighbors[x].insert(c);
			if ((x!=a)&&(x!=b)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					for (intit3 = tree->topLevel.begin(); intit3 != tree->topLevel.end(); ++intit3) {
						z = *intit3;
						if ((z!=x)&&(z!=c)) {
							jscore += deltascore(d,c,x,z);
						}
					}
					cscore += centerscore(d,c,x);
				}
				if (sm.has_uv(c,x)) {
					std::cout << "bad stuff will happen\nc = "<<c<<", x = "<<x<<"\n";
				}
				assert(sm.AddPair(c,x,jscore,cscore));
				assert(sm.AddPair(x,c,jscore,cscore));
			}
		}
		for (intit = secondNeighbors[c].begin(); intit != secondNeighbors[c].end(); ++intit) {
			x = *intit;
			secondNeighbors[x].insert(c);
			if ((x!=a)&&(x!=b)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					for (intit3 = tree->topLevel.begin(); intit3 != tree->topLevel.end(); ++intit3) {
						z = *intit3;
						if ((z!=x)&&(z!=c)) {
							jscore += deltascore(d,c,x,z);
						}
					}
					cscore += centerscore(d,c,x);
				}
				assert(sm.AddPair(c,x,jscore,cscore));
				assert(sm.AddPair(x,c,jscore,cscore));
			}
		}

		//create NEW 2nd neighbor scores, if any
		for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
			x = *intit;
			if (x!=b) {
				for (intit2 = firstNeighbors[b].begin(); intit2 != firstNeighbors[b].end(); ++intit2) {
					y = *intit2;
					if (y!=a) {
						// x && y may have just become 2nd neighbors
						if ((x!=y)&&(secondNeighbors[x].find(y)==secondNeighbors[x].end())&&(firstNeighbors[x].find(y)==firstNeighbors[x].end())) {
							// make second neighbors and add the score
							secondNeighbors[x].insert(y);
							cscore = 0; jscore = 0;
							for (d=0;d<dim;d++) {
								for (intit3 = tree->topLevel.begin(); intit3 != tree->topLevel.end(); intit3++) {
									z = *intit3;
									if ((z!=x)&&(z!=y)) {
										jscore = jscore + deltascore(d,x,y,z);
									}
								}
								cscore = cscore + centerscore(d,x,y);
							}
							assert(sm.AddPair(x,y,jscore,cscore));
						}
						if ((x!=y)&&(secondNeighbors[y].find(x)==secondNeighbors[y].end())&&(firstNeighbors[y].find(x)==firstNeighbors[y].end())) {
							// make second neighbors and add the score
							secondNeighbors[y].insert(x);
							cscore = 0; jscore = 0;
							for (d=0;d<dim;d++) {
								for (intit3 = tree->topLevel.begin(); intit3 != tree->topLevel.end(); intit3++) {
									z = *intit3;
									if ((z!=y)&&(z!=x)) {
										jscore = jscore + deltascore(d,y,x,z);
									}
								}
								cscore = cscore + centerscore(d,y,x);
							}
							assert(sm.AddPair(y,x,jscore,cscore));
						}
					}
				}
			}
		}

		// delete a,b from weights (and degrees) in all dimensions
		for (d=0;d<dim;d++) {
			if (D[d].numV<(a+1)) {
				w[d].degrees.erase(a);
				w[d].selfMissing.erase(a);
				w[d].nV.erase(a);
			}
			if (D[d].numV<(b+1)) {
				w[d].degrees.erase(b);
				w[d].selfMissing.erase(b);
				w[d].nV.erase(b);
			}
			// now let us go through the neighbors of a
			for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(a);
				w[d].dat[a].erase(x);
			}
			for (intit = secondNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(a);
				w[d].dat[a].erase(x);
			}
			w[d].dat.erase(a);
			// now let us go through the neighbors of b
			for (intit = firstNeighbors[b].begin(); intit != firstNeighbors[b].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(b);
				w[d].dat[b].erase(x);
			}
			for (intit = secondNeighbors[b].begin(); intit != secondNeighbors[b].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(b);
				w[d].dat[b].erase(x);
			}
			w[d].dat.erase(b);
		}
		// delete a,b from neighbors
		for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
			x = *intit;
			firstNeighbors[x].erase(a);
		}
		firstNeighbors.erase(a);
		for (intit = firstNeighbors[b].begin(); intit != firstNeighbors[b].end(); ++intit) {
			x = *intit;
			firstNeighbors[x].erase(b);
		}
		firstNeighbors.erase(b);
		for (intit = secondNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
			x = *intit;
			secondNeighbors[x].erase(a);
		}
		secondNeighbors.erase(a);
		for (intit = secondNeighbors[b].begin(); intit != secondNeighbors[b].end(); ++intit) {
			x = *intit;
			secondNeighbors[x].erase(b);
		}
		secondNeighbors.erase(b);
		// should be done!
		

	
		numJoins++;
		if (DEBUGMODE) {
			std::cout << "joined "<<a<<" and "<<b<<" to form "<<c<<"\n";
			assert(! (1+sm.has_u(a)));
			assert(! (1+sm.has_u(b)));
		}
	}
	return numJoins;
};


bool Engine::initializeFirstLev() {
	std::set<int> emptySet;
	curLev = 0;
	int i,d,u,v,x,y,z;
	double jscore,cscore;
	Node* pn;
	std::map<int, graphData::destList*>::iterator it1;
	graphData::destList::iterator it2;
	for (i=0;i<D[0].numV;i++) { //TODO do the tree making inside the tree constructor
		pn = new Node(i,-1,1);
		pn->theta = new float[dim];
		tree->nodeMap[i] = pn;
		tree->numNodes = i+1;
		tree->topLevel.insert(i);
		pn->vertexSet.insert(i);
		pn->collapsed = 1;
		pn->vertsComputed = 1;
	}
	// initialize weights, degrees, selfMissing and first neighbors, and nV
	for (d=0;d<dim;d++) {
		for (it1 = D[d].edgeList.begin(); it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			w[d].nV[u]=1;
			if (firstNeighbors.find(u)==firstNeighbors.end()) firstNeighbors[u] = emptySet;
			if (w[d].degrees.find(u)==w[d].degrees.end()) w[d].degrees[u]=0;
			if (w[d].selfMissing.find(u)==w[d].selfMissing.end()) w[d].selfMissing[u]=0;
			for (it2 = (*it1).second->begin(); it2 != (*it1).second->end(); ++it2) {
				v = (*it2).first;
				w[d].nV[v]=1;
				w[d].AddPair(u,v,(*it2).second);
				w[d].AddPair(v,u,(*it2).second);
				if (w[d].degrees.find(v)==w[d].degrees.end()) w[d].degrees[v]=0;
				if (w[d].selfMissing.find(v)==w[d].selfMissing.end()) w[d].selfMissing[v]=0;
				w[d].degrees[u] += (*it2).second;
				w[d].degrees[v] += (*it2).second;
				if (firstNeighbors.find(v)==firstNeighbors.end()) firstNeighbors[v] = emptySet;
				firstNeighbors[u].insert(v);
				firstNeighbors[v].insert(u);
			}
		}
	}
	// initialize 2nd neighbors
	std::map<int, Node*>::iterator itnode;
	std::set<int>::iterator neighbit;
	std::set<int>::iterator intsetit;
	for (itnode = tree->nodeMap.begin(); itnode != tree->nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			if (secondNeighbors.find(x)==secondNeighbors.end()) secondNeighbors[u] = emptySet;
			set_difference_update(secondNeighbors[x],firstNeighbors[y],firstNeighbors[x]);
			secondNeighbors[x].erase(x);
		}
	}
	// initialize scores
	for (itnode = tree->nodeMap.begin(); itnode != tree->nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			// if score(x,y) doesnt exist, compute it
			if (! sm.has_uv(x,y)) {
				// add up score contributions
				jscore = 0; cscore = 0;
				// go through all dimensions
				for (d=0;d<dim;d++) {
					// go through all top level groups z != x,y
					for (intsetit = tree->topLevel.begin(); intsetit != tree->topLevel.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore(d,x,y,z);
						}
					}
					cscore = cscore + centerscore(d,x,y);
				}
				sm.AddPair(x,y,jscore,cscore);
			}
		}
		for (neighbit = secondNeighbors[x].begin(); neighbit != secondNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			// if score(x,y) doesnt exist, compute it
			if (! sm.has_uv(x,y)) {
				// add up score contributions
				jscore = 0; cscore = 0;
				// go through all dimensions
				for (d=0;d<dim;d++) {
					// go through all top level groups z != x,y
					for (intsetit = tree->topLevel.begin(); intsetit != tree->topLevel.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore(d,x,y,z);
						}
					}
					cscore = cscore + centerscore(d,x,y);
				}
				sm.AddPair(x,y,jscore,cscore);
			}
		}
	}
	//done initializations
	return 1;
};

/*
bool Engine::doStage() {
	curLev++;
	int a,b,i,d;
	std::set<int> aset,bset,xset;
	std::set<int> *p1set;
	std::set<int> *p2set;
	std::set<int>::iterator setiter,setiter2,setiter3;
	// dataMap iterators
	for (d=dim;d<dim;d++) {
		p1set = w[d].grkeys();
		for (setiter = p1set->begin(); setiter != p1set->end(); ++setiter) aset.insert(*setiter);
		delete p1set;
	}
	// aset is the set of all groups to use as 'a' (all active groups in case of undirected graph)
	for (setiter = aset.begin(); setiter != aset.end(); ++setiter) {
		a = (*setiter);
		// let us go through all the dimensions to get the 1st and 2nd neighbors in bset
		for (d=0;d<dim;d++) {
			p1set = w[d].neighbors(a);
			for (setiter2 = p1set->begin(); setiter2 != p1set->end(); ++setiter2) {
				bset.insert(*setiter2);
				p2set = w[d].neighbors(*setiter2);
				for (setiter3 = p2set->begin(); setiter3 != p2set->end(); ++setiter3) bset.insert(*setiter3);
				delete p2set;
			}
			delete p1set;
		}
		for (setiter2 = bset.begin(); setiter2 != bset.end(); ++setiter2) {
			b = *setiter2;
			if (! sm.has_uv(a,b)) {
				// calculate score for a,b
				for (d=0;d<dim;d++) {
					// run over x
					for (setiter3 = aset.begin();
					sm.Addto(u,v,deltascore(d,a,b,x));
			
		
		


}; */

#endif
