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
#include "dataMap.hpp"


class TreeClass;


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
		float *thDen;
		float *thNum;
		std::set<int>* getAllVerts(std::map<int, Node*>& mp);
		Node(int i, int j, bool ist) {
			nid = i;
			parent = j;
			isTerm = ist;
		}
		Node(int nodeID, int parentID, bool isTerminal, int vertID, int dimension);
		Node(int nodeID, int parentID, bool isTerminal, int dimension);
		bool collapseNode(std::map<int,Node*> &nmap);
		virtual bool makeThetaforMerged(int a, int b, dataMap* w, TreeClass* tree, graphData* D);
		virtual ~Node();
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
		int dim;
		/*TreeClass(graphData* G) {
			D = G;
			numNodes=0;
		};*/
		TreeClass(graphData* G, int dimension);
		TreeClass() {}
		virtual ~TreeClass();
		Node* getNode(int n);
		virtual int makeMergeNode(int a, int b);
		virtual Node* returnNode(int i) {
			return nodeMap[i];
		}
};


Node::~Node() {
	delete[] theta;
	delete[] thNum;
	delete[] thDen;
}



bool Node::makeThetaforMerged(int a, int b, dataMap* w, TreeClass* tree, graphData* D) {
	float wc, wab;
	float x;
	for (int d=0;d<(tree->dim);d++) {
		//wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
		//assert(w[d].AddPair(this->nid,this->nid,wc));
		w[d].degrees[this->nid] = w[d].degrees[a] + w[d].degrees[b];
		w[d].selfMissing[this->nid] = w[d].selfMissing[a] + w[d].selfMissing[b] + (w[d].degrees[a] * w[d].degrees[b]) - w[d].get_uv(a,b);
		//w[d].nV[this->nid] = w[d].nV[a] + w[d].nV[b];
		w[d].nV.push_back(w[d].nV[a] + w[d].nV[b]);
		if (D[d].gtype=='w') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].degrees[a] * w[d].degrees[b];
		} else if (D[d].gtype=='b') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].nV[a] * w[d].nV[b];
		} else {
			std::cerr << "graph type "<<D[d].gtype<<" not yet supported (during theta calculation).\n";
			throw 1;
		}
		if (this->collapsed) {
			this->thNum[d] += tree->nodeMap[a]->thNum[d] + tree->nodeMap[b]->thNum[d];
			this->thDen[d] += tree->nodeMap[a]->thDen[d] + tree->nodeMap[b]->thDen[d];
		}
		x = this->thNum[d] / this->thDen[d];
		if (std::isnan(x)) x = 0;
		this->theta[d] = x;
		//TODO:: delete the following, only for debugging
		if (theta<0) {
			std::cerr << "bad theta being written!\n";
		}
	}
	return 1;
}

bool Node::collapseNode(std::map<int,Node*> &nmap) {
	this->collapsed = 1;
	// will return 1 if all children have their verts computed and we can compute for this node as well, 0 if children don't have their verts computed
	std::set<int>::iterator nit;
	std::set<int>::iterator setit;
	int allChildVerts = 1;
	for (nit=this->childSet.begin();
	    (nit != this->childSet.end()) and (allChildVerts);
	    nit++) {
		if (! nmap[ *nit ]->vertsComputed)
			allChildVerts = 0;
	}
	if (allChildVerts) {
		Node* cnode;
		for (nit=this->childSet.begin(); nit != this->childSet.end(); nit++) {
			cnode = nmap[ *nit ];
			for (setit = cnode->vertexSet.begin(); setit != cnode->vertexSet.end(); ++setit) {
				this->vertexSet.insert( *setit );
			}
		}
		this->vertsComputed = 1;
		return 1;
	}
	return 0;
};

Node::Node(int nodeID, int parentID, bool isTerminal, int dimension) {
	this->nid = nodeID;
	this->parent = parentID;
	this->isTerm = isTerminal;
	this->theta = new float[dimension];
	this->thDen = new float[dimension];
	this->thNum = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->theta[d] = 0;
		this->thDen[d] = 0;
		this->thNum[d] = 0;
	}
	if (isTerminal) {
		this->collapsed = 1;
	}
}

Node::Node(int nodeID, int parentID, bool isTerminal, int vertID, int dimension) {
	if (! isTerminal) {
		std::cerr << "bad call to Node constructor!, given vertex ID for non-terminal node!\n";
		throw(1);
	}
	this->nid = nodeID;
	this->parent = parentID;
	this->isTerm = isTerminal;
	this->theta = new float[dimension];
	this->thDen = new float[dimension];
	this->thNum = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->theta[d] = 0;
		this->thDen[d] = 0;
		this->thNum[d] = 0;
	}
	if (isTerminal) {
		this->collapsed = 1;
	}
	this->vertexSet.insert(vertID);
	this->vertsComputed = 1;
}

	

int TreeClass::makeMergeNode(int a, int b) {
	int c = this->numNodes;
	this->numNodes = c + 1;
	Node* pnode = new Node(c, -1, 0, dim);
	this->nodeMap[c] = pnode;
	this->nodeMap[a]->parent = c;
	this->nodeMap[b]->parent = c;
	pnode->childSet.insert(a);
	pnode->childSet.insert(b);
	this->topLevel.erase(a);
	this->topLevel.erase(b);
	this->topLevel.insert(c);
	return c;
}


TreeClass::TreeClass(graphData* G, int dimension) {
	this->dim = dimension;
	D = G;
	for (int i=0; i<D[0].numV; i++) {
		this->nodeMap[i] = new Node(i, -1, 1, i, this->dim);
		this->topLevel.insert(i);
	}
	this->numNodes = D[0].numV;
}
TreeClass::~TreeClass() {
	std::map<int, Node*>::iterator nit;
	for (nit=nodeMap.begin(); nit != nodeMap.end(); ++nit) {
		delete (*nit).second;
	}
	nodeMap.clear();
	topLevel.clear();
};



Node* TreeClass::getNode(int i) {
	Node* n;
	n = nodeMap[i];
	return n;
}


std::pair<int,int> TreeClass::getLCA(int i, int j) {
	int u,v;
	u = i; v = j;
	int udone = 0; int vdone = 0;
	while (1) {
		if (u>v) {
			if (! vdone) {
				if (topLevel.find(v)!=topLevel.end()) {
					vdone = 1;
					// v is done, so u has to go to topLevel also
					while (! udone) {
						if (topLevel.find(u) == topLevel.end()) {
							u = nodeMap[u]->parent;
						} else {
							udone = 1;
						}
					}
				} else {
					v = nodeMap[v]->parent;
				}
			}
		} else {
			if (! udone) {
				if (topLevel.find(u)!=topLevel.end()) {
					udone = 1;
					// u is done, so v has to go to topLevel also
					while (! vdone) {
						if (topLevel.find(v) == topLevel.end()) {
							v = nodeMap[v]->parent;
						} else {
							vdone = 1;
						}
					}
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
