#ifndef NODETREE_HPP
#define NODETREE_HPP
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <fstream>


class Node{
	public:
		int parent;
		int nid;
		std::set<int> childSet;
		std::set<int> vertexSet;
		bool isTerm;
		bool collapsed;
		Node(int i, int j, bool ist) {
			nid = i;
			parent = j;
			isTerm = ist;
		}
};
		
class TreeClass{
	public:
		std::map<int, Node*> nodeMap;
		std::set<int> topLevel;
		int numNodes;
		bool writeCompHierEdges(const char* fn);
		bool writeHRG(const char* fn);
		bool writeNodeTypes(const char* fn);
		TreeClass() {
			numNodes=0;
		};
		~TreeClass() {
			std::map<int, Node*>::iterator nit;
			for (nit=nodeMap.begin(); nit != nodeMap.end(); ++nit) {
				delete (*nit).second;
			}
		};
};		


bool TreeClass::writeNodeTypes(const char* fn) {
	std::ofstream file;
	std::map<int, Node*>::iterator mapit;
	file.open(fn,std::ios::out);
	file << "ID\tType\n";
	if (not file.is_open()) return 0;
	for (mapit = nodeMap.begin(); mapit != nodeMap.end(); mapit++) {
		// process the node
		file << (*mapit).first << "\t";
		if ((*mapit).second->isTerm) file << "Vertex\n";
			else file << "Internal\n";
	}
	return 1;
};



bool TreeClass::writeCompHierEdges(const char* fn) {
	std::stack<int> st;
	int n;
	std::ofstream file;
	std::set<int>::iterator intit,intit2;
	file.open(fn,std::ios::out);
	if (not file.is_open()) return 0;
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
