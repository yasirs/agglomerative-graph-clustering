#ifndef ENGINE_HPP
#define ENGINE_HPP
#include <map>
#if ISVC
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"
#include "scoremap.hpp"
#include "dataMap.hpp"
#include "mysetfuncs.hpp"
#include "dataMapOther.hpp"
#include "nodetreeOther.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
// #include <boost/math/special_functions/gamma.hpp>


struct FPairHash {
  std::size_t operator() (const std::pair<float,float>& p) const {
    return ((int) (1000*p.first) + (int) (100*p.second) );
  }
};

struct FPairEqual {
  bool operator() (const std::pair<float,float>& a, const std::pair<float,float>& b) const {
	if (a.first != b.first) return 0;
	else if (a.second != b.second) return 0;
	else return 1;
  }
};

float mySafeLog(float x) {
	if (x<1e-8) return 1e10;
	else return log(x);
};



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


float gammaFunction(float a) {
	if (a==1) return 1;
	else return gamma(a);
};

float gammaFunction(int a) {
	if (a==1) return 1;
	else return gamma(a);
};

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
		Engine(graphData* G, graphData* Goriginal, graphData* GsoFar, int d);
		bool initializeScores();
		int run();
		bool doStage();
		bool cleanUp();
		float deltascore(int d, int a, int b, int x);
		float centerscore(int d, int a, int b);
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
		for (v=u+1; v<D[d].numV; v++) {
			if ((edges)||(! D[d].has_uv(u,v))) {
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
		for (v=u+1; v<D[d].numV; v++) {
			if ((edges)||(! D[d].has_uv(u,v))) {
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
	float s,dummy;
	std::ofstream file;
	file.open(fn,std::ios::out);
	t = D[d].numV-2;
	for (u=0; u<D[d].numV; u++) {
		n = D[d].degree(u);
		for (v=u+1; v<D[d].numV; v++) {
			if ((edges)||(! D[d].has_uv(u,v))) {
				m = D[d].degree(v);
				dmin = std::min(m,n);
				c = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
				s = 0;
				for (x = c; x<= dmin; x++) {
					dummy = 1.0;
					dummy = dummy / gammaFunction(1+x);
					dummy = dummy / gammaFunction(1+m-x);
					dummy = dummy / gammaFunction(1+n-x);
					dummy = dummy / gammaFunction(1+t-m-n+x);
					s = s + dummy;
					
				}
				s = s * gammaFunction(1+m) * gammaFunction(1+n) * gammaFunction(1+t-m) * gammaFunction(1+t-n) / gammaFunction(1+t);
				s = -log10(s);
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
		for (v=u+1; v<D[d].numV; v++) {
			if ((edges)||(! D[d].has_uv(u,v))) {
				bd = D[d].degree(v);
				file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << ad*bd << '\n';
			}
		}
	}
	file.close();
};




Engine::Engine(graphData* G, graphData* Goriginal, graphData* GsoFar, int d) {
	D = G;
	dim = d;
	w = new dataMapOther[dim];
	tree = new TreeClassOther(G,d);
	int x,y;
	// initialize weights, degrees, selfMissing and first neighbors, and nV
	for (d=0; d<dim; d++) {
		( (dataMapOther*)  &w[d]   )->initialize(D[d], Goriginal[d], GsoFar[d], firstNeighbors);
	}
	std::set<int> emptySet;
	// initialize 2nd neighbors
	std::map<int, Node*>::iterator itnode;
	std::set<int>::iterator neighbit;
	std::set<int>::iterator intsetit;
	for (itnode = tree->nodeMap.begin(); itnode != tree->nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		if (secondNeighbors.find(x)==secondNeighbors.end()) secondNeighbors[x] = emptySet;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			set_difference_update(secondNeighbors[x],firstNeighbors[y],firstNeighbors[x]);
		}
		secondNeighbors[x].erase(x);
	}
};










Engine::Engine(graphData* G, int d) {
	D = G;
	dim = d;
	w = new dataMap[dim];
	tree = new TreeClass(G,d);
	int x,y;
	// initialize weights, degrees, selfMissing and first neighbors, and nV
	for (d=0; d<dim; d++) {
		w[d].initialize(D[d], firstNeighbors);
	}
	std::set<int> emptySet;
	/*for (d=0;d<dim;d++) {
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
	}*/
	// initialize 2nd neighbors
	std::map<int, Node*>::iterator itnode;
	std::set<int>::iterator neighbit;
	std::set<int>::iterator intsetit;
	for (itnode = tree->nodeMap.begin(); itnode != tree->nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		if (secondNeighbors.find(x)==secondNeighbors.end()) secondNeighbors[x] = emptySet;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			set_difference_update(secondNeighbors[x],firstNeighbors[y],firstNeighbors[x]);
		}
		secondNeighbors[x].erase(x);
	}
};


Engine::~Engine() {
	delete[] w;
	delete tree;
};


float Engine::centerscore(int d, int a, int b) {
	assert(a!=b);
	float ans;
	if (D[d].gtype=='b') {
		float Tab, Eab, Eaa, Ebb, Taa, Tbb, Tcc;
		Tab = w[d].nV[a] * w[d].nV[b];
		Taa = w[d].nV[a] * (w[d].nV[a] - 1)/2.0f;
		Tbb = w[d].nV[b] * (w[d].nV[b] - 1)/2.0f;
		Tcc = Taa + Tbb + Tab;
		Eab = w[d].get_uv(a,b);
		Eaa = w[d].get_uv(a,a);
		Ebb = w[d].get_uv(b,b);
		float Etot = Eaa+Ebb+Eab;
		float Ttot = Taa+Tbb+Tab;
		float thetaTot = Etot/Ttot;
		float thetaA = Eaa/Taa;
		float thetaB = Ebb/Tbb;
		float thetaAB = Eab/Tab;
		ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
			- Eaa * mySafeLog(thetaA) - (Taa - Eaa) * mySafeLog(1-thetaA)
			- Ebb * mySafeLog(thetaB) - (Tbb - Ebb) * mySafeLog(1-thetaB)
			- Eab * mySafeLog(thetaAB) - (Tab - Eab) * mySafeLog(1-thetaAB);
		
	/*	ans =   lnBetaFunction(Eaa + Eab + Ebb + 1,Haa + Hbb + Tab - Eab + 1)
			- lnBetaFunction(Eaa + 1, Haa +1)
			- lnBetaFunction(Ebb + 1, Hbb +1)
			- lnBetaFunction(Eab+1,Tab - Eab +1);
		//for the random reference model
#if (! NOREFERENCE)
		ans -=  lnBetaFunction(Tcc*D[d].aveP + 1, Tcc*(1 - D[d].aveP) + 1)
			- lnBetaFunction(Taa*D[d].aveP + 1, Taa*(1 - D[d].aveP) + 1)
			- lnBetaFunction(Tbb*D[d].aveP + 1, Tbb*(1 - D[d].aveP) + 1)
			- lnBetaFunction(Tab*D[d].aveP + 1, Tab*(1 - D[d].aveP) + 1);
#endif*/
	} else if (D[d].gtype=='w') {
		float Tab, Eab, Eaa, Ebb, Haa, Hbb, Taa, Tbb;
		Eab = w[d].get_uv(a,b);
		Eaa = w[d].get_uv(a,a);
		Ebb = w[d].get_uv(b,b);
		Haa = w[d].selfMissing[a];
		Hbb = w[d].selfMissing[a];
		Tab = w[d].degrees[a] * w[d].degrees[b];
		Taa = Eaa + Haa;
		Tbb = Ebb + Hbb;
		float Etot = Eaa+Ebb+Eab;
		float Ttot = Taa+Tbb+Tab;
		float thetaTot = Etot/Ttot;
		float thetaA = Eaa/Taa;
		float thetaB = Ebb/Tbb;
		float thetaAB = Eab/Tab;
		ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
			- Eaa * mySafeLog(thetaA) - (Taa - Eaa) * mySafeLog(1-thetaA)
			- Ebb * mySafeLog(thetaB) - (Tbb - Ebb) * mySafeLog(1-thetaB)
			- Eab * mySafeLog(thetaAB) - (Tab - Eab) * mySafeLog(1-thetaAB);
		/*
		ans =   lnBetaFunction(Eaa + Eab + Ebb + 1,Haa + Hbb + Tab - Eab + 1)
			- lnBetaFunction(Eaa + 1, Haa +1)
			- lnBetaFunction(Ebb + 1, Hbb +1)
			- lnBetaFunction(Eab+1,Tab - Eab +1);
		//for the random reference model
#if (! NOREFERENCE)
		Tcc = Taa + Tbb + Tab;
		float Ebaraa, Ebarbb, Ebarab, Ebarcc;
		Ebarab = w[d].degrees[a] * w[d].degrees[b];
		Ebaraa = Taa/(2.0f*D[d].Etot);
		Ebarbb = Tbb/(2.0f*D[d].Etot);
		Ebarcc = Ebaraa + Ebarbb + Ebarbb;
		ans -=  lnBetaFunction(Ebarcc + 1, Tcc - Ebarcc + 1)
			- lnBetaFunction(Ebaraa + 1, Taa - Ebaraa + 1)
			- lnBetaFunction(Ebarbb + 1, Tbb - Ebarbb + 1)
			- lnBetaFunction(Ebarab + 1, Tab - Ebarab + 1);	
#endif*/
	}
	return ans;
};


float Engine::deltascore(int d, int a, int b, int x) {
	// maximum likelihood version
	assert((a!=b)&&(a!=x)&&(b!=x));
	float ans;
	float da, db, dx;
	if (D[d].gtype=='b') {
		float Eax, Ebx, Tax, Tbx;
		Eax = w[d].get_uv(a,x);
		Ebx = w[d].get_uv(b,x);
		da = w[d].nV[a];
		db = w[d].nV[b];
		dx = w[d].nV[x];
		Tax = da * dx;
		Tbx = db * dx;
		//Hax = Tax - Eax;
		//Hbx = Tbx - Ebx;
		float Etot = Eax+Ebx;
		float Ttot = Tax+Tbx;
		float thetaTot = Etot/Ttot;
		float thetaA = Eax/Tax;
		float thetaB = Ebx/Tbx;
		ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
			- Eax * mySafeLog(thetaA) - (Tax - Eax) * mySafeLog(1-thetaA)
			- Ebx * mySafeLog(thetaB) - (Tbx - Ebx) * mySafeLog(1-thetaB);
		//TODO::
		if (std::isnan(ans)) {
			std::cout << "isnan answer!\n";
		}
	} else if (D[d].gtype=='w') {
		std::cout << "max likelihood doesn't work for weighted graphs, I think\n";
		float Eax, Ebx, Tax, Tbx;
		Eax = w[d].get_uv(a,x);
		Ebx = w[d].get_uv(b,x);
		da = w[d].degrees[a];
		db = w[d].degrees[b];
		dx = w[d].degrees[x];
		Tax = da * dx;
		Tbx = db * dx;
		/*Ebarax = da * dx/(2.0f*D[d].Etot);
		Ebarbx = db * dx/(2.0f*D[d].Etot);
		Hbarax = da * dx - Ebarax;
		Hbarbx = db * dx - Ebarbx;*/
		float Etot = Eax+Ebx;
		float Ttot = Tax+Tbx;
		float thetaTot = Etot/Ttot;
		float thetaA = Eax/Tax;
		float thetaB = Ebx/Tbx;
		ans = Etot * mySafeLog(thetaTot) + (Ttot - Etot) * mySafeLog(1-thetaTot)
			- Eax * mySafeLog(thetaA) - (Tax - Eax) * mySafeLog(1-thetaA)
			- Ebx * mySafeLog(thetaB) - (Tbx - Ebx) * mySafeLog(1-thetaB);
/*
#if (! NOREFERENCE)
		ans -=	 (	lnBetaFunction(Ebarax+Ebarbx+1.0f,Hbarax+Hbarbx+1.0f)
				-lnBetaFunction(Ebarax+1.0f,Hbarax+1.0f)
				-lnBetaFunction(Ebarbx+1.0f,Hbarbx+1.0f)
			 );
#endif*/
	}
	else {
		std::cout << "graph type "<<D[d].gtype<<" not understood.\n";
		throw 1;
	}
	return ans;
};


int Engine::run() {
	std::set<int> emptySet, tempSet;
	std::set<int> possSet1, possSet2;
	std::set<int>::iterator intit, intit2, intit3, intsetit;
	int numJoins = 0;
	int a,b,c,x,y,z,d, tint;
	float cscore, jscore;
	scoremap::pairScore pscore;
	std::map<int, scoremap::smap>::iterator smOut;
	std::map<int, scoremap::twoScores>::iterator smIn;
	std::map<int, std::map<int,float> >::iterator datOut;
	std::map<int, float>::iterator datIn;
	//while (sm.hasPos()) {
	while (! sm.isEmpty()) {
		// join and do stuff
		pscore = sm.popBestScore();
		a = pscore.u; b = pscore.v;
		sm.erase(b,a);
		// let us create new group c and create heirarchical relations
                c = tree->makeMergeNode(a,b);
		/*c = tree->numNodes;
		c++;
		tree->numNodes = c;
		pnode = new Node(c,-1,0);
		pnode->theta = new float[dim];
		pnode->thNum = new float[dim];
		pnode->thDen = new float[dim];
		tree->nodeMap[c] = pnode;
		tree->nodeMap[a]->parent = c;
		tree->nodeMap[b]->parent = c;
		tree->nodeMap[c]->childSet.insert(a);
		tree->nodeMap[c]->childSet.insert(b);
		tree->topLevel.erase(a);
		tree->topLevel.erase(b);
		tree->topLevel.insert(c);
		*/
		if ((tree->nodeMap[a]->collapsed)&&(tree->nodeMap[b]->collapsed)&&(pscore.s.centerMscore>0)) {
                        assert( tree->nodeMap[c]->collapseNode(tree->nodeMap) );
			/*pnode->collapsed = 1;
			for(intit= tree->nodeMap[a]->vertexSet.begin(); intit!= tree->nodeMap[a]->vertexSet.end(); ++intit) {
				pnode->vertexSet.insert(*intit);
			}
			for(intit= tree->nodeMap[b]->vertexSet.begin(); intit!= tree->nodeMap[b]->vertexSet.end(); ++intit) {
				pnode->vertexSet.insert(*intit);
			}
			pnode->vertsComputed = 1;*/
		} else {
			tree->nodeMap[c]->collapsed = 0;
			tree->nodeMap[c]->vertsComputed = 0;

		}
		// compute d,w,n,m for c
                tree->nodeMap[c]->makeDataforMerged(a,b,w,tree,D);
		/*for (d=0;d<dim;d++) {
			wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
			assert(w[d].AddPair(c,c,wc));
			w[d].degrees[c] = w[d].degrees[a] + w[d].degrees[b];
			w[d].selfMissing[c] = w[d].selfMissing[a] + w[d].selfMissing[b] + (w[d].degrees[a] * w[d].degrees[b]) - w[d].get_uv(a,b);
			w[d].nV[c] = w[d].nV[a] + w[d].nV[b];
			if (D[d].gtype=='w') {
				wab = w[d].get_uv(a,b);
				tree->nodeMap[c]->thNum[d] = wab;
				tree->nodeMap[c]->thDen[d] = w[d].degrees[a] * w[d].degrees[b];
			} else if (D[d].gtype=='b') {
				wab = w[d].get_uv(a,b);
				tree->nodeMap[c]->thNum[d] = wab;
				tree->nodeMap[c]->thDen[d] = w[d].nV[a] * w[d].nV[b];
			} else {
				std::cerr << "graph type "<<D[d].gtype<<" not yet supported (during theta calculation).\n";
				throw 1;
			}
			if (pnode->collapsed) {
				pnode->thNum[d] += tree->nodeMap[a]->thNum[d] + tree->nodeMap[b]->thNum[d];
				pnode->thDen[d] += tree->nodeMap[a]->thDen[d] + tree->nodeMap[b]->thDen[d];
			}
			theta = pnode->thNum[d] / pnode->thDen[d];
			if (std::isnan(theta)) theta = 0;
			tree->nodeMap[c]->theta[d] = theta;
			//TODO:: delete the following, only for debugging
			if (theta<0) {
				std::cerr << "bad theta being written!\n";
			}

		}*/
		firstNeighbors[c] = emptySet;
		secondNeighbors[c] = emptySet;
		set_union_update(firstNeighbors[c],firstNeighbors[a],firstNeighbors[b]);
		firstNeighbors[c].erase(a); firstNeighbors[c].erase(b); 
		tempSet = emptySet;
		set_union_update(tempSet,secondNeighbors[a],secondNeighbors[b]);
		set_difference_update(secondNeighbors[c],tempSet,firstNeighbors[c]);
		secondNeighbors[c].erase(a); secondNeighbors[c].erase(b);
		for (d=0;d<dim;d++) {
			w[d].addMergedData(a,b,c,firstNeighbors[c]);
			/*for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
				x = (*intit);
				w[d].AddPair(c,x,w[d].get_uv(a,x) + w[d].get_uv(b,x));
				w[d].AddPair(x,c,w[d].get_uv(x,a) + w[d].get_uv(b,x));
			}*/
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
		// add the c to the relevent groups neighbors!
		for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
			x = (*intit);
			firstNeighbors[x].insert(c);
		}
		for (intit = secondNeighbors[c].begin(); intit != secondNeighbors[c].end(); ++intit) {
			x = (*intit);
			secondNeighbors[x].insert(c);
		}

		// update all affected scores (at least one of nodes has to be a neighbor of c)
		for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
			x = (*intit);
			smOut = sm.scores.find(x);
			if (smOut != sm.scores.end()) {
				for (smIn = (*smOut).second.scoreDest.begin(); smIn != (*smOut).second.scoreDest.end(); ++smIn) {
					y = (*smIn).first;
					if (y!=c) {
						jscore =0;
						for (d=0;d<dim;d++) {
							jscore = jscore + deltascore(d,x,y,c)-deltascore(d,x,y,a)-deltascore(d,x,y,b);
						}
						if (firstNeighbors[c].count(y)) {
							assert(sm.AddTo(x,y,jscore/2));
						} else {
							assert(sm.AddTo(x,y,jscore));
						}
					}
				}
			}
		}
		
		//TODO:: debug check
		assert(sm.has_u(c)<0);
		// create new c,x scores
		for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
			x = *intit;
			if ((x!=a)&&(x!=b)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					// go through the union set of neighbors
					std::set<int> neighbUnion;
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[c]);
					for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(c==z)) ) {
							jscore = jscore + deltascore(d,c,x,z);
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
					// go through the union set of neighbors
					std::set<int> neighbUnion;
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[c]);
					for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(c==z)) ) {
							jscore = jscore + deltascore(d,x,c,z);
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
								// go through the union set of neighbors
								std::set<int> neighbUnion;
								set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
								for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
									z = *intsetit;
									if (! ((x==z)||(y==z)) ) {
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
								// go through the union set of neighbors
								std::set<int> neighbUnion;
								set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
								for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
									z = *intsetit;
									if (! ((x==z)||(y==z)) ) {
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
		for (d=0; d<dim; d++) {
			assert(w[d].allErase(a,b, D[d].numV));
		}
		/*for (d=0;d<dim;d++) {
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
		}*/
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
			tint = sm.has_u(a);
			if (tint!=-1) {
				std::cerr << "Error! deleted vertex "<<a<<" has score with "<<tint<<" !\n";
				throw 1;
			}
			tint = sm.has_u(b);
			if (tint!=-1) {
				std::cerr << "Error! deleted vertex "<<b<<" has score with "<<tint<<" !\n";
				throw 1;
			
			}
		}
	}
	return numJoins;
};


bool Engine::initializeScores() {
	std::set<int> emptySet;
	curLev = 0;
	int d,x,y,z;
	float jscore,cscore;
	std::map<int, graphData::destList*>::iterator it1;
	graphData::destList::iterator it2;
	std::map<int, Node*>::iterator itnode;
	std::set<int>::iterator neighbit;
	std::set<int>::iterator intsetit;
	/*for (i=0;i<D[0].numV;i++) { //TODO do the tree making inside the tree constructor
		pn = new Node(i,-1,1);
		pn->theta = new float[dim];
		pn->thNum = new float[dim];
		pn->thDen = new float[dim];
		for (d=0;d<dim;d++) {
			pn->theta[d] = 0;
			pn->thDen[d] = 0;
			pn->thNum[d] = 0;
		}
		tree->nodeMap[i] = pn;
		tree->numNodes = i+1;
		tree->topLevel.insert(i);
		pn->vertexSet.insert(i);
		pn->collapsed = 1;
		pn->vertsComputed = 1;
	}*/
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
					// go through the union set of neighbors
					std::set<int> neighbUnion;
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
					for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore(d,x,y,z);
						}
					}
					cscore = cscore + centerscore(d,x,y);
				}
				assert(sm.AddPair(x,y,jscore,cscore));
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
					// go through the union set of neighbors
					std::set<int> neighbUnion;
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
					for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						z = *intsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore(d,x,y,z);
						}
					}
					cscore = cscore + centerscore(d,x,y);
				}
				assert(sm.AddPair(x,y,jscore,cscore));
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
