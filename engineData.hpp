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
		int* NodeMembership;
		TreeClass tree;
		Engine(graphData* G, int d);
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
			}
		}
	}
	return 1;
};


bool Engine::doStage() {
	curLev++;
	std::map<int, std::map<int,float> >::iterator it1;
	std::map<int,float>::iterator it2;
	

};

#endif
