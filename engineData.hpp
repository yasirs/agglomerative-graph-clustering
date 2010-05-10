#ifndef ENGINEDATA_HPP
#define ENGINEDATA_HPP
#include <map>
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"
#include "scoremap.hpp"
#include "dataMap.hpp"

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
		std::map<int, float> groupDegrees;
		std::map<int, float> selfMissing;
		
		Engine(graphData* G, int d);
		double deltascore(int d, int a, int b, int x);
		bool initializeFirstLev();
		int run();
		bool doStage();
		bool cleanUp();
		double deltaScore(int d, int a, int b, int x);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
		
};

int Engine::run() {
	int numJoins = 0;
	int a,b,c;
	Node* pnode;
	scoremap::pairScore pscore;
	while (sm.hasPos()) {
		// join and do stuff
		pscore = sm.popBestScore();
		a = pscore.u; b = pscore.v;
		// let us create new group c
		c = tree.numNodes;
		c++;
		tree.numNodes = c;
		pnode = new Node(-1,c);
		tree.nodeVec[a]->parent = c;
		tree.nodeVec[b]->parent = c;
		c.childSet.insert(a);
		c.childSet.insert(b);
		tree.topLevel.erase(a);
		tree.topLevel.erase(b);
		tree.topLevel.insert(c);



		numJoins++;
	}
	return numJoins;
};


bool Engine::initializeFirstLev() {
	std::set<int> emptySet;
	curLev = 0;
	int i,d,u,v,x,y,z;
	double jscore,cscore;
	Node* pn;
	for (i=0;i<G[0].numV;i++) {
		pn = new Node(i,-1);
		tree.nodeVec.push_back(pn);
		tree.topLevel.insert(tree.nodeVec.size());
	}
	std::map<int,graphData::destList>::iterator it1;
	destList::iterator it2;
	// initialize weights, degrees, selfMissing and first neighbors
	for (d=0;d<dim;d++) {
		for (it1 = D[d].edgeList.begin(); it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			if (firstNeighbors.find(u)==firstNeighbors.end()) firstNeighbors[u] = emptySet;
			if (groupDegrees.find(u)==groupDegrees.end()) groupDegrees[u]=0;
			if (selfMissing.find(u)==selfMissing.end()) selfMissing[u]=0;
			for (it2 = (*it1).second.begin(); it2 != (*it1).second.end(); ++it2) {
				v = (*it2).first;
	
				w[d].AddPair(u,v,(*it2).second);
				w[d].AddPair(v,u,(*it2).second);
				if (groupDegrees.find(v)==groupDegrees.end()) groupDegrees[v]=0;
				if (selfMissing.find(v)==selfMissing.end()) selfMissing[v]=0;
				groupDegrees[u] += (*it2).second;
				groupDegrees[v] += (*it2).second;
				if (firstNeighbors.find(v)==firstNeighbors.end()) firstNeighbors[v] = emptySet;
				firstNeighbors[u].insert(v);
				firstNeighbors[v].insert(u);
			}
		}
	}
	// initialize 2nd neighbors
	std::vector<Node*>::iterator itnode;
	std::set<int>::iterator neighbit;
	std::set<int>::iterator intsetit;
	for (itnode = tree.nodeVec.begin(); itnode != tree.nodeVec.end(); ++itnode) {
		x = itnode->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != firstNeighbors[x].end(); ++neighit) {
			y = *neighbit;
			if (secondNeighbors.find(x)==secondNeighbors.end()) secondNeighbors[u] = emptySet;
			set_difference_update(secondNeighbors[x],firstNeighbors[y],firstNeighbors[x]);
		}
	}
	// initialize scores
	for (itnode = tree.nodeVec.begin(); itnode != tree.nodeVec.end(); ++itnode) {
		x = itnode->nid;
		for (neighbit = firstNeighbors[x].begin(); neighbit != secondNeighbors[x].end(); ++neighit) {
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
