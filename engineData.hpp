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
		Engine(graphData* G, int d);
		scoremap::twoScores deltascore(int d, int a, int b, int x);
		bool initializeFirstLev();
		bool doStage();
		bool cleanUp();
		double deltaScore(int d, int a, int b, int x);
		std::set<int> getNeighborsVertex(int i);
		std::set<int> getNeighborsNode(int i);
};


bool Engine::initializeFirstLev() {
	curLev = 0;
	int i,d,u,v;
	Node* pn;
	for (i=0;i<G[0].numV;i++) {
		pn = new Node(-1);
		tree.nodeVec.push_back(pn);
		tree.topLevel.insert(tree.nodeVec.size());
	}
	it1 = std::map<int,graphData::destList>::iterator;
	it2 = destList::Iterator;
	for (d=0;d<dim;d++) {
		for (it1 = D[d].edgeList.begin(); it1 != D[d].edgeList.end(); ++it1) {
			u = (*it1).first;
			for (it2 = (*it1).second.begin(); it2 != (*it1).second.end(); ++it2) {
				v = (*it2).first;
				w[d].AddPair(u,v,(*it2).second);
				w[d].AddPair(v,u,(*it2).second);
			}
		}
	}
	return 1;
};


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
			
		
		


};

#endif
