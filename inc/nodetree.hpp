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
#include "modelParam.hpp"


class TreeClass;


class Node{
	public:
		int parent;
		int nid;
		int party;
		std::set<int> childSet;
		std::set<int> vertexSet;
		bool isTerm;
		bool collapsed;
		bool vertsComputed;
		ModelParamBase** params;
		std::set<int>* getAllVerts(std::map<int, Node*>& mp);
		Node(int i, int j, int partytype, bool ist) {
			nid = i;
			parent = j;
			isTerm = ist;
			party = partytype;
		}
		Node(int nodeID, int parentID, int partytype, bool isTerminal, int vertID, int dimension, graphData* D);
		Node(int nodeID, int parentID, int partytype, bool isTerminal, int dimension, graphData* D);
		bool collapseNode(std::map<int,Node*> &nmap);
		virtual bool writeThetaforMerged(int a, int b, dataMap** w, TreeClass* tree, graphData* D);
		virtual void destroy(int di);
		virtual ~Node();
		virtual bool isOther() {return 0;}
};

class TreeClass{
	public:
		virtual bool isOther() {return false; }
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

		TreeClass(graphData* G, int dimension);
		TreeClass() {}
		virtual ~TreeClass();
		Node* getNode(int n);
		virtual int makeMergeNode(int a, int b);
		virtual Node* returnNode(int i) {
			return nodeMap[i];
		}
};

void Node::destroy(int di) {
	int d;
	for (d=0;d<di;d++) {
		delete params[d];
	}
}



Node::~Node() {
	delete[] params;
}



bool Node::writeThetaforMerged(int a, int b, dataMap** w, TreeClass* tree, graphData* D) {
	for (int d=0;d<(tree->dim);d++) {
		this->params[d]->calculate(w[d]->get_uv(a,b),w[d]->datvert[a],w[d]->datvert[b]);
		if (this->collapsed) {
			this->params[d]->collapse(tree->nodeMap[a]->params[d],tree->nodeMap[b]->params[d]);
		}
		this->params[d]->cleanup();
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

Node::Node(int nodeID, int parentID, int partytype, bool isTerminal, int dimension, graphData* D) {
	this->party = partytype;
	this->nid = nodeID;
	this->parent = parentID;
	this->isTerm = isTerminal;
	this->params = new ModelParamBase*[dimension];
	for (int d=0; d<dimension; d++) {
		if (D[d].gtype == 'b') {
			this->params[d] = new BinomialParam;
		} else if (D[d].gtype == 'p') {
			this->params[d] = new PoissonParam;
		} else if (D[d].gtype == 'w') {
			this->params[d] = new WParam;
		} else if (D[d].gtype == 'd') {
			this->params[d] = new DcorrParam;
		} else if (D[d].gtype == 'g') {
			this->params[d] = new GaussianParam;
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
		if (isTerminal)
			this->params[d]->init();
	}
	if (isTerminal) {
		this->collapsed = 1;
	}
}

Node::Node(int nodeID, int parentID, int partytype, bool isTerminal, int vertID, int dimension, graphData* D) {
	this->party = partytype;
	if (! isTerminal) {
		std::cerr << "bad call to Node constructor!, given vertex ID for non-terminal node!\n";
		throw(1);
	}
	this->nid = nodeID;
	this->parent = parentID;
	this->isTerm = isTerminal;
	this->params = new ModelParamBase* [dimension];
	for (int d=0; d<dimension; d++) {
		if (D[d].gtype == 'b') {
			this->params[d] = new BinomialParam;
		} else if (D[d].gtype == 'p') {
			this->params[d] = new PoissonParam;
		} else if (D[d].gtype =='w') {
			this->params[d] = new WParam;
		} else if (D[d].gtype =='d') {
			this->params[d] = new DcorrParam;
		} else if (D[d].gtype =='g') {
			this->params[d] = new GaussianParam;
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
		this->params[d]->init();
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
	int ptype = this->nodeMap[a]->party;
	assert(ptype == this->nodeMap[b]->party);
	Node* pnode = new Node(c, -1, ptype, 0, dim,D);
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
	for (unsigned int i=0; i<D[0].numV; i++) {
		this->nodeMap[i] = new Node(i, -1, D[0].typeList[i], 1, i, this->dim, D);
		this->topLevel.insert(i);
	}
	this->numNodes = D[0].numV;
}
TreeClass::~TreeClass() {
	std::map<int, Node*>::iterator nit;
	for (nit=nodeMap.begin(); nit != nodeMap.end(); ++nit) {
		nit->second->destroy(dim);
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
