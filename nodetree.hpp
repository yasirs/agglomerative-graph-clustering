#ifndef NODETREE_HPP
#define NODETREE_HPP
#include <vector>
#include <set>


class Node{
	public:
		int parent;
		int nid;
		std::set<int> childSet;
		std::set<int> vertexSet;
		Node(int i, int j) {
			nid = i;
			parent = j;
		}
};
		
class TreeClass{
	public:
		std::vector<Node*> nodeVec;
		std::set<int> topLevel;
		int numNodes;
		TreeClass() {
			numNodes=0;
		}
};		


#endif
