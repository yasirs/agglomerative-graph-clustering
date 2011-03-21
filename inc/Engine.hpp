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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

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
		dataMap** w;
		scoremap sm;
		int* NodeMembership;
		TreeClass *tree;
		std::map<int,std::set<int> > firstNeighbors;
		std::map<int,std::set<int> > secondNeighbors;
		//std::map<int, float> groupDegrees;
		struct mergeRecord {
			int child1; int child2; int merged;
			mergeRecord(int i, int j, int k) {
				child1 = i; child2 = j; merged = k;
			}
		};
		std::list<mergeRecord> mergeList;
		~Engine();
		Engine(graphData* G, int d);
		Engine(graphData* G, graphData* Goriginal, graphData* GsoFar, int d);
		bool initializeScoresML();
		int runML();
		int runFB();
		void passFB();
		void passML();
		bool doStage();
		bool cleanUp();
		float deltascoreML(int d, int a, int b, int x);
		float deltascoreFB(int d, int a, int b, int x);
		//float debug_getscore(int d,int a,int b);
		float centerscoreML(int d, int a, int b);
		float centerscoreFB(int d, int a, int b);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
		void printJaccardFile(const char* filename, int d, bool skipEdges);
		void printHyperGeomFile(const char* filename, int d, bool skipEdges);
		void printDegreeProdFile(const char* filename, int d, bool skipEdges);
		void printCommonNeighbFile(const char* filename, int d, bool skipEdges);
};


void Engine::printCommonNeighbFile(const char* fn, int d, bool skipEdges) {
	int u,v, c;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		for (v=u+1; v<D[d].numV; v++) {
			if (skipEdges) if (D[d].has_uv(u,v)) continue;
			c = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
			file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << c << '\n';
		}
	}
	file.close();
};


void Engine::printJaccardFile(const char* fn, int d, bool skipEdges) {
	int u,v, aub, aib, ad;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		ad = D[d].degree(u);
		for (v=u+1; v<D[d].numV; v++) {
			if (skipEdges) if (D[d].has_uv(u,v)) continue;
			aib = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
			aub = (ad * D[d].degree(v)) - aib;
			file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << (aib+0.0f)/aub << '\n';
		}
	}
	file.close();
};

void Engine::printHyperGeomFile(const char* fn, int d, bool skipEdges) {
	unsigned int u, v, m, n, c, t;
	float s;
	std::ofstream file;
	file.open(fn,std::ios::out);
	t = D[d].numV-2;
	for (u=0; u<D[d].numV; u++) {
		n = D[d].degree(u);
		for (v=u+1; v<D[d].numV; v++) {
			if (skipEdges) if (D[d].has_uv(u,v)) continue;
			m = D[d].degree(v);
			//dmin = std::min(m,n);
			c = num_common_keys( *(D[d].edgeList[u]), *(D[d].edgeList[v]) );
			s = gsl_cdf_hypergeometric_Q(c, n, t-n, m) + gsl_ran_hypergeometric_pdf(c, n, t-n, m);
			s = -log10(s);
			file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << s << '\n';
		}
	}
	file.close();
};

void Engine::printDegreeProdFile(const char* fn, int d, bool skipEdges) {
	int u,v, ad, bd;
	std::ofstream file;
	file.open(fn,std::ios::out);
	for (u=0; u<D[d].numV; u++) {
		ad = D[d].degree(u);
		for (v=u+1; v<D[d].numV; v++) {
			if (skipEdges) if (D[d].has_uv(u,v)) continue;
			bd = D[d].degree(v);
			file << D[d].int2Name[u] << '\t' << D[d].int2Name[v] << '\t' << ad*bd << '\n';
		}
	}
	file.close();
};




Engine::Engine(graphData* G, graphData* Goriginal, graphData* GsoFar, int d) {
	D = G;
	dim = d;
	w = new dataMap*[dim];
	tree = new TreeClassOther(G,dim);
	int x,y;
	// initialize weights, degrees, selfMissing and first neighbors, and nV
	for (int dd=0; dd<dim; dd++) {
		w[dd] = new dataMapOther;
		( (dataMapOther*)  w[dd]   )->initialize(D[dd], Goriginal[dd], GsoFar[dd], firstNeighbors);
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
	w = new dataMap*[dim];
	tree = new TreeClass(G,d);
	int x,y;
	// initialize weights, degrees, selfMissing and first neighbors, and nV
	for (d=0; d<dim; d++) {
		w[d] = new dataMap;
		w[d]->initialize(D[d], firstNeighbors);
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


Engine::~Engine() {
	for (int i=0; i<dim; i++) {
		delete w[i];
	}
	delete[] w;
	delete tree;
};


float Engine::centerscoreFB(int d, int a, int b) {
	assert(a!=b);
	return this->w[d]->FBcenterscore(a,b);
};


float Engine::deltascoreFB(int d, int a, int b, int x) {
	assert((a!=b)&&(a!=x)&&(b!=x));
	return this->w[d]->FBdeltascore(a,b,x);
};



float Engine::centerscoreML(int d, int a, int b) {
	assert(a!=b);
	return this->w[d]->MLcenterscore(a,b);
};


float Engine::deltascoreML(int d, int a, int b, int x) {
	assert((a!=b)&&(a!=x)&&(b!=x));
	return this->w[d]->MLdeltascore(a,b,x);
};


int Engine::runML() {
	std::set<int> emptySet, tempSet;
	std::set<int> possSet1, possSet2;
	std::set<int>::iterator intit, intit2, intit3, intsetit;
	int numJoins = 0;
	int a,b,c,x,y,d, tint;
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
			w[d]->addMergedData(a,b,c,firstNeighbors[c]);
		}

		tree->nodeMap[c]->writeThetaforMerged(a,b,w,tree,D);

		// delete a,b scores
		for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
			x = *intit;
			if ((b!=x)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[a]->party)) {
				sm.erase(a,x);
				sm.erase(x,a);
			}
		}
		for (intit = secondNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
			x = *intit;
			if ((b!=x)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[a]->party)) {
				sm.erase(a,x);
				sm.erase(x,a);
			}
		}
		for (intit = firstNeighbors[b].begin(); intit != firstNeighbors[b].end(); ++intit) {
			x = *intit;
			if ((a!=x)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[b]->party)) {
				sm.erase(b,x);
				sm.erase(x,b);
			}
		}
		for (intit = secondNeighbors[b].begin(); intit != secondNeighbors[b].end(); ++intit) {
			x = *intit;
			if ((a!=x)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[b]->party)) {
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
							if ((not D[d].multiPartite)or(this->tree->nodeMap[x]->party != this->tree->nodeMap[c]->party))
								jscore = jscore + deltascoreML(d,x,y,c)-deltascoreML(d,x,y,a)-deltascoreML(d,x,y,b);
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
			if ((x!=a)&&(x!=b)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[c]->party)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					// go through the union set of neighbors
					set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[c]);
					for (int z; neighbEnum.next(z);) {
						if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[c]->party)) {
							if (! ((x==z)||(c==z)||(a==z)||(b==z)) ) {
								jscore = jscore + deltascoreML(d,c,x,z);
							}
						}
					}
					cscore += centerscoreML(d,c,x);
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
			if ((x!=a)&&(x!=b)&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[c]->party)) {
				cscore = 0;
				jscore = 0;
				for (d=0;d<dim;d++) {
					// go through the union set of neighbors
					set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[c]);
					for (int z; neighbEnum.next(z);) {
						if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[c]->party)) {
							if (! ((x==z)||(c==z)||(a==z)||(b==z)) ) {
								jscore = jscore + deltascoreML(d,x,c,z);
							}
						}
					}
					cscore += centerscoreML(d,c,x);
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
						if ((x!=y)&&(secondNeighbors[x].find(y)==secondNeighbors[x].end())&&(firstNeighbors[x].find(y)==firstNeighbors[x].end())&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[y]->party)) {
							// make second neighbors and add the score
							secondNeighbors[x].insert(y);
							cscore = 0; jscore = 0;
							for (d=0;d<dim;d++) {
								// go through the union set of neighbors
								set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[y]);
								for (int z; neighbEnum.next(z);) {
									if (! ((x==z)||(y==z)) ) {
										if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[x]->party)) {
											jscore = jscore + deltascoreML(d,x,y,z);
										}
									}
								}
								cscore = cscore + centerscoreML(d,x,y);
							}
							assert(sm.AddPair(x,y,jscore,cscore));
						}
						if ((x!=y)&&(secondNeighbors[y].find(x)==secondNeighbors[y].end())&&(firstNeighbors[y].find(x)==firstNeighbors[y].end())&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[y]->party)) {
							// make second neighbors and add the score
							secondNeighbors[y].insert(x);
							cscore = 0; jscore = 0;
							for (d=0;d<dim;d++) {
								// go through the union set of neighbors
								set_union_Enumerator<int> neighbEnum(firstNeighbors[x],firstNeighbors[y]);
								for (int z; neighbEnum.next(z);) {
									if (! ((x==z)||(y==z)) ) {
										if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[x]->party)) {
											jscore = jscore + deltascoreML(d,y,x,z);
										}
									}
								}
								cscore = cscore + centerscoreML(d,y,x);
							}
							assert(sm.AddPair(y,x,jscore,cscore));
						}
					}
				}
			}
		}

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
		this->mergeList.push_back(mergeRecord(a,b,c));
		if (DEBUGMODE) {
			std::cout << "joined ";
			if (D->int2Name.find(a) == D->int2Name.end()) std::cout << a; else std::cout << D->int2Name[a];
			std::cout  <<" and ";
			if (D->int2Name.find(b) == D->int2Name.end()) std::cout << b; else std::cout << D->int2Name[b];
			std::cout <<" to form "<< c <<", for score = "<< pscore.s.joinScore <<"\n";
			tint = sm.has_u(a);
			if (tint!=-1) {
				std::cerr << "Error! deleted vertex number "<< a <<" has score with "<<tint<<" !\n";
				throw 1;
			}
			tint = sm.has_u(b);
			if (tint!=-1) {
				std::cerr << "Error! deleted vertex number "<< b <<" has score with "<<tint<<" !\n";
				throw 1;
			
			}
		}
	}
	return mergeList.size();
};


bool Engine::initializeScoresML() {
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
	// initialize scores
	for (itnode = tree->nodeMap.begin(); itnode != tree->nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			// if score(x,y) doesnt exist, compute it
			if ((! sm.has_uv(x,y))&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[y]->party)) {
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
							if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[x]->party)) {
								jscore = jscore + deltascoreML(d,x,y,z);
							}
						}
					}
					cscore = cscore + centerscoreML(d,x,y);
				}
				assert(sm.AddPair(x,y,jscore,cscore));
			}
		}
		for (neighbit = secondNeighbors[x].begin(); neighbit != secondNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			// if score(x,y) doesnt exist, compute it
			if ((! sm.has_uv(x,y))&&(this->tree->nodeMap[x]->party == this->tree->nodeMap[y]->party)) {
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
							if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[x]->party)) {
								jscore = jscore + deltascoreML(d,x,y,z);
							}
						}
					}
					cscore = cscore + centerscoreML(d,x,y);
				}
				assert(sm.AddPair(x,y,jscore,cscore));
			}
		}
	}
	//done initializations
	return 1;
};

void Engine::passFB() {
	int a,b,c,d,z;
	float csc, dsc;
	// initialize toplevel to the vertices
	this->tree->topLevel.clear();
	for (int i=0; i<D[0].numV; i++) {
		tree->topLevel.insert(i);
	}
	std::list<mergeRecord>::iterator mergeit(mergeList.begin());
	std::set<int>::iterator intit;
	for (; mergeit != mergeList.end(); mergeit++) {
		a = mergeit->child1;
		b = mergeit->child2;
		c = mergeit->merged;
		// decide whether to do this merge
		dsc = 0;
		for (d=0;d<dim;d++) {
			for (intit = tree->topLevel.begin(); intit != tree->topLevel.end(); intit++) {
				z = *intit;
				if ((a!=z) & (b!=z) & (c!=z)) {
					if ((not D[d].multiPartite)or(this->tree->nodeMap[z]->party != this->tree->nodeMap[a]->party)) {
						dsc += deltascoreFB(d,a,b,z);
					}
				}
			}
		}
		if (dsc>=0) {
			// let us do the merge
			tree->topLevel.erase(a);
			tree->topLevel.erase(b);
			tree->topLevel.insert(c);
			
			// do we need to merge them?
			if ((tree->nodeMap[a]->collapsed) & (tree->nodeMap[b]->collapsed)) {
				// let us compute the center score
				csc = 0;
				for (d=0;d<dim;d++) {
					csc += centerscoreFB(d,a,b);
				}
				if (csc>=0) {
					// let us merge them
					assert(tree->nodeMap[c]->collapseNode(tree->nodeMap));
				} else {
					tree->nodeMap[c]->collapsed = 0;
					tree->nodeMap[c]->vertsComputed = 0;
				}
			}
			tree->nodeMap[c]->writeThetaforMerged(a,b,w,tree,D);
		} else {
			// break out of the loop
			break;
		}
		if (DEBUGMODE) {
			if (dsc>=0) {
				std::cout << "passFB: joined ";
				if (D->int2Name.find(a) == D->int2Name.end()) std::cout << a; else std::cout << D->int2Name[a];
				std::cout  <<" and ";
				if (D->int2Name.find(b) == D->int2Name.end()) std::cout << b; else std::cout << D->int2Name[b];
				std::cout <<" to form "<< c <<", for score = "<< dsc <<"\n";
				if (csc>=0) {
					std::cout << "passFB: collapsed "<< c << ", for centerscore = " << csc <<"\n";
				}
			}
			else {
				std::cout << "passFB: did not join ";
				if (D->int2Name.find(a) == D->int2Name.end()) std::cout << a; else std::cout << D->int2Name[a];
				std::cout  <<" and ";
				if (D->int2Name.find(b) == D->int2Name.end()) std::cout << b; else std::cout << D->int2Name[b];
				std::cout <<" to form "<< c <<", for score = "<< dsc <<"\n";
			
			}
		}
	}
					
}



#endif
