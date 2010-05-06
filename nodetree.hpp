#ifndef NODETREE_HPP
#define NODETREE_HPP
#include <vector>
#include <set>


class Node{
	public:
		int parent;
		std::set<int> childSet;
		std::set<int> vertexSet;
		Node(int i) {
			parent = i;
		}
};
		
class TreeClass{
	public:
		std::vector<Node*> nodeVec;
		std::set<int> topLevel;
};		


#endif
