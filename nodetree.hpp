#ifndef NODETREE_HPP
#define NODETREE_HPP
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <utility>
#include <fstream>
#include <iostream>
#include "graphData.hpp"

class Node{
	public:
		int parent;
		int nid;
		std::set<int> childSet;
		std::set<int> vertexSet;
		bool isTerm;
		bool collapsed;
		bool vertsComputed;
		float *theta;
		std::set<int>* getAllVerts(std::map<int, Node*>& mp);
		Node(int i, int j, bool ist) {
			nid = i;
			parent = j;
			isTerm = ist;
		}
};



	
class TreeClass{
	public:
		graphData* D;
		std::map<int, Node*> nodeMap;
		std::set<int> topLevel;
		int numNodes;
		bool writeCompHierEdges(const char* fn);
		bool writeCollapsedHierEdges(const char* fn);
		bool writeHRG(const char* fn);
		bool writeNodeTypes(const char* fn);
		std::pair<int,int> getLCA(const int i, const int j);
		TreeClass(graphData* G) {
			D = G;
			numNodes=0;
		};
		~TreeClass() {
			std::map<int, Node*>::iterator nit;
			for (nit=nodeMap.begin(); nit != nodeMap.end(); ++nit) {
				delete (*nit).second;
			}
		};
};


std::pair<int,int> TreeClass::getLCA(int i, int j) {
	int u,v;
	u = i; v = j;
	int udone = 0; int vdone = 0;
	while (1) {
		if (u>v) {
			if (! vdone) {
				if (topLevel.find(v)==topLevel.end()) {
					vdone = 1;
				} else {
					v = nodeMap[v]->parent;
				}
			}
		} else {
			if (! udone) {
				if (topLevel.find(u)==topLevel.end()) {
					udone = 1;
				} else {
					u = nodeMap[u]->parent;
				}
			}
		}
		if ((u==v)||((udone)&&(vdone))) return std::pair<int,int>(u,v);
	}
};


bool TreeClass::writeNodeTypes(const char* fn) {
	std::ofstream file;
	std::map<int, Node*>::iterator mapit;
	file.open(fn,std::ios::out);
	file << "ID\tType\n";
	if (! file.is_open()) return 0;
	for (mapit = nodeMap.begin(); mapit != nodeMap.end(); mapit++) {
		// process the node
		if ((*mapit).second->isTerm) {
			file << D[0].int2Name[(*mapit).first] << "\tVertex\n";
			if (topLevel.find((*mapit).first)!=topLevel.end()) {
				file << "I" << (*mapit).first << "\tInternal\n";
			}
		} else {
			file << "I" << (*mapit).first << "\tInternal\n";
		}
	}
	return 1;
};



bool TreeClass::writeCollapsedHierEdges(const char* fn) {
	std::stack<int> st;
	int n;
	std::ofstream file;
	std::set<int>::iterator intit,intit2;
	file.open(fn,std::ios::out);
	if (! file.is_open()) return 0;
	for (intit = topLevel.begin(); intit != topLevel.end(); intit++) {
		if (nodeMap[*intit]->isTerm) {
			file << "I"<<(*intit) << "\t" << D[0].int2Name[(*intit)] << "\n";
		} else {
			st.push(*intit);
			while ((st.size()>0)) {
				n = st.top(); st.pop();
				if (nodeMap[n]->collapsed) {
					if (nodeMap[n]->vertexSet.size()>1) {
						for (intit2 = nodeMap[n]->vertexSet.begin(); intit2 != nodeMap[n]->vertexSet.end(); intit2++) {
							file << "I"<<n << "\t" << D[0].int2Name[*intit2] << "\n";
						}
					}
				} else {
					for (intit2 = nodeMap[n]->childSet.begin(); intit2 != nodeMap[n]->childSet.end(); ++intit2) {
						st.push(*intit2);
						// process the element
						if (nodeMap[*intit2]->isTerm) {
							file << "I"<<n << "\t" << D[0].int2Name[*intit2] << "\n";
						} else {
							file << "I"<<n << "\t" <<"I" << (*intit2) << "\n";
						}
							
					}
				}
			}
		}
	}
	file.close();
	return 1;
};



bool TreeClass::writeCompHierEdges(const char* fn) {
	std::stack<int> st;
	int n;
	std::ofstream file;
	std::set<int>::iterator intit,intit2;
	file.open(fn,std::ios::out);
	if (! file.is_open()) return 0;
	for (intit = topLevel.begin(); intit != topLevel.end(); intit++) {
		st.push(*intit);
		while ((st.size()>0)) {
			n = st.top(); st.pop();
			for (intit2 = nodeMap[n]->childSet.begin(); intit2 != nodeMap[n]->childSet.end(); ++intit2) {
				st.push(*intit2);
				// process the element
				file << n << "\t" << (*intit2) << "\n";
			}
		}
	}
	file.close();
	return 1;
};
			

#endif
