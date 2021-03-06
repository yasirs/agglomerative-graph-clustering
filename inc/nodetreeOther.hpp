#ifndef NODETREEOTHER_HPP
#define NODETREEOTHER_HPP
#include "nodetree.hpp"
#include "dataMapOther.hpp"
#include "modelParam.hpp"

class TreeClassOther;

class NodeOther: public Node {
	public:
		ModelParamBase** paramsOriginal;
		ModelParamBase** paramsSoFar;
		virtual bool writeThetaforMerged(int a, int b, dataMap** w, TreeClass* tree, graphData* D);
		NodeOther(int nodeID, int parentID, int partytype, bool isTerminal, int vertID, int dimension, graphData* D);
		NodeOther(int nodeID, int parentID, int partytype, bool isTerminal, int dimension, graphData* D);
		virtual void destroy(int di);
		virtual ~NodeOther();
		virtual bool isOther() {return 1;}
};


class TreeClassOther: public TreeClass {
	public:
		virtual bool isOther() { return true; }
		virtual int makeMergeNode(int a, int b);
		TreeClassOther(graphData* G, int dimension);
		virtual ~TreeClassOther();
		NodeOther* returnNode(int i) {
			return (NodeOther*) nodeMap[i];
		}
};

TreeClassOther::~TreeClassOther() {
	// nothing to do
	// the base class destructor should be called after this
}

int TreeClassOther::makeMergeNode(int a, int b) {
	int c = this->numNodes;
	this->numNodes = c + 1;
	int ptype = this->nodeMap[a]->party;
	assert(ptype == this->nodeMap[b]->party);
	NodeOther* pnode = new NodeOther(c, -1, ptype, 0, dim, D);
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



TreeClassOther::TreeClassOther(graphData* G, int dimension) {
	this->dim = dimension;
	D = G;
	for (unsigned int i=0; i<D[0].numV; i++) {
		this->nodeMap[i] = new NodeOther(i, -1, D[0].typeList[i], 1, i, this->dim, D);
		this->topLevel.insert(i);
	}
	this->numNodes = D[0].numV;
}

NodeOther::NodeOther(int nodeID, int parentID, int partytype, bool isTerminal, int dimension, graphData* D) :Node(nodeID, parentID, partytype, isTerminal, dimension, D) {
	this->paramsOriginal = new ModelParamBase*[dimension];
	this->paramsSoFar = new ModelParamBase*[dimension];
	for (int d=0; d<dimension; d++) {
		if (D[d].gtype == 'b') {
			this->paramsOriginal[d] = new BinomialParam;
			this->paramsSoFar[d] = new BinomialParam;
		} else if (D[d].gtype == 'p') {
			this->paramsOriginal[d] = new PoissonParam;
			this->paramsSoFar[d] = new PoissonParam;
		} else if (D[d].gtype == 'w') {
			this->paramsOriginal[d] = new WParam;
			this->paramsSoFar[d] = new WParam;
		} else if (D[d].gtype == 'd') {
			this->paramsOriginal[d] = new DcorrParam;
			this->paramsSoFar[d] = new DcorrParam;
		} else if (D[d].gtype == 'g') {
			this->paramsOriginal[d] = new GaussianParam;
			this->paramsSoFar[d] = new GaussianParam;
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
		if (isTerminal) {
			this->paramsOriginal[d]->init();
			this->paramsSoFar[d]->init();
		}
	}
}

NodeOther::NodeOther(int nodeID, int parentID, int partytype, bool isTerminal, int vertID, int dimension, graphData* D): Node(nodeID, parentID, partytype, isTerminal, vertID, dimension, D) {
	if (! isTerminal) {
		std::cerr << "bad call to NodeOther constructor!, given vertex ID for non-terminal node!\n";
		throw(1);
	}
	this->paramsOriginal = new ModelParamBase*[dimension];
	this->paramsSoFar = new ModelParamBase*[dimension];
	for (int d=0; d<dimension; d++) {
		if (D[d].gtype == 'b') {
			this->paramsOriginal[d] = new BinomialParam;
			this->paramsSoFar[d] = new BinomialParam;
		} else if (D[d].gtype == 'p') {
			this->paramsOriginal[d] = new PoissonParam;
			this->paramsSoFar[d] = new PoissonParam;
		} else if (D[d].gtype == 'w') {
			this->paramsOriginal[d] = new WParam;
			this->paramsSoFar[d] = new WParam;
		} else if (D[d].gtype == 'd') {
			this->paramsOriginal[d] = new DcorrParam;
			this->paramsSoFar[d] = new DcorrParam;
		} else if (D[d].gtype == 'g') {
			this->paramsOriginal[d] = new GaussianParam;
			this->paramsSoFar[d] = new GaussianParam;
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
		paramsOriginal[d]->init();
		paramsSoFar[d]->init();
	}
}


void NodeOther::destroy(int di) {
	Node::destroy(di);
	int d;
	for (d=0;d<di;d++) {
		delete paramsSoFar[d];
		delete paramsOriginal[d];
	}
}



NodeOther::~NodeOther() {
	delete[] paramsSoFar;
	delete[] paramsOriginal;
}

bool NodeOther::writeThetaforMerged(int a, int b, dataMap** ww, TreeClass* tree, graphData* D) {
	dataMapOther **w = (dataMapOther**) ww;

	for (int d=0;d<(tree->dim);d++) {
		this->params[d]->calculate(w[d]->get_uv(a,b),w[d]->datvert[a],w[d]->datvert[b]);
		this->paramsOriginal[d]->calculate(w[d]->get_uvOriginal(a,b),w[d]->oDatvert[a],w[d]->oDatvert[b]);
		this->paramsSoFar[d]->calculate(w[d]->get_uvSoFar(a,b),w[d]->sDatvert[a],w[d]->sDatvert[b]);
		if (this->collapsed) {
			this->params[d]->collapse(((NodeOther*) tree->nodeMap[a])->params[d], ((NodeOther*) tree->nodeMap[b])->params[d]);
			this->paramsOriginal[d]->collapse(((NodeOther*) tree->nodeMap[a])->paramsOriginal[d], ((NodeOther*) tree->nodeMap[b])->paramsOriginal[d]);
			this->paramsSoFar[d]->collapse(((NodeOther*) tree->nodeMap[a])->paramsSoFar[d], ((NodeOther*) tree->nodeMap[b])->paramsSoFar[d]);
		}
		this->paramsOriginal[d]->cleanup();
		this->paramsSoFar[d]->cleanup();
		this->params[d]->cleanup();
		//this->params[d]->bestfromSoFar(this->paramsOriginal[d],this->paramsSoFar[d]);
		
	}
	return 1;
}







#endif
