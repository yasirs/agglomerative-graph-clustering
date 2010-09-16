#ifndef ENGINE_HPP
#define ENGINE_HPP
#include <map>
#if ISVC
#include <unordered_map>
#include <unordered_set>
#else
#include <tr1/unordered_map>
#include <tr1/unordered_set>
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
		float deltascore1(int d, int a, int b, int x) {
			return deltascore(d,a,b,x);
		}
		float debug_getscore(int d,int a,int b);
		float centerscore(int d, int a, int b);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
		void printJaccardFile(const char* filename, int d, bool edges);
		void printHyperGeomFile(const char* filename, int d, bool edges);
		void printDegreeProdFile(const char* filename, int d, bool edges);
		void printCommonNeighbFile(const char* filename, int d, bool edges);
};


float Engine::debug_getscore(int d,int a,int b) {
	float s = 0;
	for(std::set<int>::iterator nit (tree->topLevel.begin()); nit != tree->topLevel.end(); ++nit) {
		if (! ((*nit==a)|(*nit==b))  )
		s = s + deltascore(d,a,b,*nit);
	}
	return s;
}


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
	return this->w[d].MLcenterscore(a,b);
};


float Engine::deltascore(int d, int a, int b, int x) {
	// maximum likelihood version
	assert((a!=b)&&(a!=x)&&(b!=x));
	return this->w[d].MLdeltascore(a,b,x);
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
	scoremap::scDestType::iterator smIn;
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
		if ((tree->nodeMap[a]->collapsed)&&(tree->nodeMap[b]->collapsed)&&(pscore.s.centerMscore>=0)) {
                        assert( tree->nodeMap[c]->collapseNode(tree->nodeMap) );
		} else {
			tree->nodeMap[c]->collapsed = 0;
			tree->nodeMap[c]->vertsComputed = 0;

		}
		// compute neighbours
		firstNeighbors[c] = emptySet;
		secondNeighbors[c] = emptySet;
		set_union_update(firstNeighbors[c],firstNeighbors[a],firstNeighbors[b]);
		firstNeighbors[c].erase(a); firstNeighbors[c].erase(b); 
		tempSet = emptySet;
		set_union_update(tempSet,secondNeighbors[a],secondNeighbors[b]);
		set_difference_update(secondNeighbors[c],tempSet,firstNeighbors[c]);
		secondNeighbors[c].erase(a); secondNeighbors[c].erase(b);

		// compute d,w,n,m for c

		for (d=0;d<dim;d++) {
			w[d].addMergedData(a,b,c,firstNeighbors[c]);
		}

		tree->nodeMap[c]->writeThetaforMerged(a,b,w,tree,D);

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
							assert(sm.AddTo(y,x,jscore/2));
						} else {
							assert(sm.AddTo(x,y,jscore));
							assert(sm.AddTo(y,x,jscore));
						}
					}
				}
			}
		}
		
		//TODO:: debug check
		//assert(sm.has_u(c)<0);
		// create new c,x scores
		for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
			x = *intit;
			if ((x!=a)&&(x!=b)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					// go through the union set of neighbors
					//std::set<int> neighbUnion;
					//set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[c]);
					set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[c]);
					for (int z; neighbEnum.next(z);) {
					//for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						//z = *intsetit;
						if (! ((x==z)||(c==z)||(a==z)||(b==z)) ) {
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
					//std::set<int> neighbUnion;
					//set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[c]);
					set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[c]);
					for (int z; neighbEnum.next(z);) {
					//for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
						//z = *intsetit;
						if (! ((x==z)||(c==z)||(a==z)||(b==z)) ) {
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
								//std::set<int> neighbUnion;
								//set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
								set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[y]);
								for (int z; neighbEnum.next(z);) {
								//for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
									//z = *intsetit;
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
								//std::set<int> neighbUnion;
								//set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
								set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[y]);
								for (int z; neighbEnum.next(z);) {
								//for (intsetit = neighbUnion.begin(); intsetit != neighbUnion.end(); ++intsetit) {
									//z = *intsetit;
									if (! ((x==z)||(y==z)) ) {
										jscore = jscore + deltascore1(d,y,x,z);
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
		//for (d=0; d<dim; d++) {
		//	assert(w[d].allErase(a,b, D[d].numV));
		//}
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
	//std::tr1::unordered_set<int> neighbUnion;
	//std::tr1::unordered_set<int>::iterator unsetit;
	std::set<int> neighbUnion;
	std::set<int>::iterator unsetit;
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
					neighbUnion.clear();
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
					for (unsetit = neighbUnion.begin(); unsetit != neighbUnion.end(); ++unsetit) {
						z = *unsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore1(d,x,y,z);
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
					neighbUnion.clear();
					set_union_update(neighbUnion, firstNeighbors[x], firstNeighbors[y]);
					for (unsetit = neighbUnion.begin(); unsetit != neighbUnion.end(); ++unsetit) {
						z = *unsetit;
						if (! ((x==z)||(y==z)) ) {
							jscore = jscore + deltascore1(d,x,y,z);
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
