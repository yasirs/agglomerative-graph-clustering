#ifndef ENGINEDATA_HPP
#define ENGINEDATA_HPP
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


#define DEBUGMODE 1

class Engine{
	public:
		graphData *D;
		int dim;
		int curLev;
		dataMap* w;
		scoremap sm;
		int* NodeMembership;
		TreeClass tree;
		std::map<int,std::set<int> > firstNeighbors;
		std::map<int,std::set<int> > secondNeighbors;
		//std::map<int, float> groupDegrees;
		
		Engine(graphData* G, int d);
		double deltascore(int d, int a, int b, int x);
		bool initializeFirstLev();
		int run();
		bool doStage();
		bool cleanUp();
		double deltaScore(int d, int a, int b, int x);
		double centerscore(int d, int a, int b);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
		
};

int Engine::run() {
	std::set<int> emptySet, tempSet;
	std::set<int> possSet1, possSet2;
	std::set<int>::iterator intit, intit2, intit3;
	int numJoins = 0;
	int a,b,c,x,y,z,d;
	float wc,mc,dc;
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
		// let us create new group c and create heirarchical relations
		c = tree.numNodes;
		c++;
		tree.numNodes = c;
		pnode = new Node(-1,c);
		tree.nodeMap[c] = pnode;
		tree.nodeMap[a]->parent = c;
		tree.nodeMap[b]->parent = c;
		tree.nodeMap[c]->childSet.insert(a);
		tree.nodeMap[c]->childSet.insert(b);
		tree.topLevel.erase(a);
		tree.topLevel.erase(b);
		tree.topLevel.insert(c);
		//TODO: collapsed
		// compute d,w,n,m for c
		for (d=0;d<dim;d++) {
			wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
			assert(w[d].AddPair(c,c,wc));
		}
		firstNeighbors[c] = emptySet;
		secondNeighbors[c] = emptySet;
		set_union_update(firstNeighbors[c],firstNeighbors[a],firstNeighbors[b]);
		tempSet = emptySet;
		set_union_update(tempSet,secondNeighbors[a],secondNeighbors[b]);
		set_difference_update(secondNeighbors[c],tempSet,firstNeighbors[c]);
		for (d=0;d<dim;d++) {
			for (intit = firstNeighbors[c].begin(); intit != firstNeighbors[c].end(); ++intit) {
				w[d].AddPair(c,x,w[d].get_uv(a,x) + w[d].get_uv(b,x));
				w[d].AddPair(x,c,w[d].get_uv(x,a) + w[d].get_uv(b,x));
			}
		}
		// delete a,b scores
		for (intit = firstNeighbors[a].begin(); intit != secondNeighbors[a].end(); ++intit) {
			if (intit==firstNeighbors[a].end()) intit = secondNeighbors[a].begin(); //NOTE: this is a hack
			x = *intit;
			sm.erase(a,x);
		}
		for (intit = firstNeighbors[b].begin(); intit != secondNeighbors[b].end(); ++intit) {
			if (intit==firstNeighbors[b].end()) intit = secondNeighbors[b].begin(); //NOTE: this is a hack
			x = *intit;
			sm.erase(b,x);
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
		

		// create new c,x scores
		for (intit = firstNeighbors[c].begin(); intit != secondNeighbors[c].end(); ++intit) {
			if (intit==firstNeighbors[c].end()) intit = secondNeighbors[c].begin(); //NOTE: this is a hack
			x = *intit;
			cscore = 0;
			jscore = 0;
			for (d=0;d<dim;d++) {
				for (intit3 = tree.topLevel.begin(); intit3 != tree.topLevel.end(); ++intit3) {
					z = *intit3;
					if ((z!=x)and(z!=c)) {
						jscore += deltascore(d,c,x,z);
					}
				}
				cscore += centerscore(d,c,x);
			}
			assert(sm.AddPair(c,x,jscore,cscore));
		}

		//create NEW 2nd neighbor scores, if any
		for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
			for (intit2 = firstNeighbors[b].begin(); intit2 != firstNeighbors[b].end(); ++intit2) {
				x = *intit; y = *intit2;
				// x and y may have just become 2nd neighbors
				if (secondNeighbors[x].find(y)==secondNeighbors[x].end()) {
					// make second neighbors and add the score
					secondNeighbors[x].insert(y);
					cscore = 0; jscore = 0;
					for (d=0;d<dim;d++) {
						for (intit3 = tree.topLevel.begin(); intit3 != tree.topLevel.end(); intit3++) {
							z = *intit3;
							if ((z!=x)and(z!=y)) {
								jscore = jscore + deltascore(d,x,y,z);
							}
						}
						cscore = cscore + centerscore(d,x,y);
					}
					assert(sm.AddPair(x,y,jscore,cscore));
				}
				if (secondNeighbors[y].find(x)==secondNeighbors[y].end()) {
					// make second neighbors and add the score
					secondNeighbors[y].insert(x);
					cscore = 0; jscore = 0;
					for (d=0;d<dim;d++) {
						for (intit3 = tree.topLevel.begin(); intit3 != tree.topLevel.end(); intit3++) {
							z = *intit3;
							if ((z!=y)and(z!=x)) {
								jscore = jscore + deltascore(d,y,x,z);
							}
						}
						cscore = cscore + centerscore(d,y,x);
					}
					assert(sm.AddPair(y,x,jscore,cscore));
				}
			}
		}

		// delete a,b from weights (and degrees) in all dimensions
		for (d=0;d<dim;d++) {
			w[d].degrees.erase(a);
			w[d].degrees.erase(b);
			w[d].selfMissing.erase(a);
			w[d].selfMissing.erase(b);
			// now let us go through the neighbors of a
			for (intit = firstNeighbors[a].begin(); intit != firstNeighbors[a].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(a);
			}
			w[d].dat.erase(a);
			// now let us go through the neighbors of b
			for (intit = firstNeighbors[b].begin(); intit != firstNeighbors[b].end(); ++intit) {
				x = *intit;
				w[d].dat[x].erase(b);
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
			std::cout << "joined "<<a<<", and "<<b<<"to form "<<c<<"\n";
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
	std::map<int,graphData::destList>::iterator it1;
	graphData::destList::iterator it2;
	for (i=0;i<D[0].numV;i++) {
		pn = new Node(i,-1);
		tree.nodeMap[0] = pn;
		tree.topLevel.insert(tree.nodeMap.size());
	}
	// initialize weights, degrees, selfMissing and first neighbors
	for (d=0;d<dim;d++) {
		for (it1 = D[d].edgeList.begin(); it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			if (firstNeighbors.find(u)==firstNeighbors.end()) firstNeighbors[u] = emptySet;
			if (w[d].degrees.find(u)==w[d].degrees.end()) w[d].degrees[u]=0;
			if (w[d].selfMissing.find(u)==w[d].selfMissing.end()) w[d].selfMissing[u]=0;
			for (it2 = (*it1).second.begin(); it2 != (*it1).second.end(); ++it2) {
				v = (*it2).first;
	
				w[d].AddPair(u,v,(*it2).second);
				w[d].AddPair(v,u,(*it2).second);
				if (w[d].degrees.find(u)==w[d].degrees.end()) w[d].degrees[u]=0;
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
	for (itnode = tree.nodeMap.begin(); itnode != tree.nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighbit) {
			y = *neighbit;
			if (secondNeighbors.find(x)==secondNeighbors.end()) secondNeighbors[u] = emptySet;
			set_difference_update(secondNeighbors[x],firstNeighbors[y],firstNeighbors[x]);
		}
	}
	// initialize scores
	for (itnode = tree.nodeMap.begin(); itnode != tree.nodeMap.end(); ++itnode) {
		x = (*itnode).second->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != secondNeighbors[x].end(); ++neighbit) {
			// NOTE: this is a hack to go through two sets, 1st and 2nd neighbrs
			if (neighbit == firstNeighbors[x].end()) neighbit = secondNeighbors[x].begin();
			y = *neighbit;
			// if score(x,y) doesnt exist, compute it
			if (not sm.has_uv(x,y)) {
				// add up score contributions
				jscore = 0; cscore = 0;
				// go through all dimensions
				for (d=0;d<dim;d++) {
					// go through all top level groups z != x,y
					for (intsetit = tree.topLevel.begin(); intsetit != tree.topLevel.end(); ++intsetit) {
						z = *intsetit;
						if (not ((x==z)or(y==z)) ) {
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
			if (not sm.has_uv(a,b)) {
				// calculate score for a,b
				for (d=0;d<dim;d++) {
					// run over x
					for (setiter3 = aset.begin();
					sm.Addto(u,v,deltascore(d,a,b,x));
			
		
		


}; */

#endif
