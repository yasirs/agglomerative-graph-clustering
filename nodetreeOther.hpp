#ifndef NODETREEOTHER_HPP
#define NODETREEOTHER_HPP
#include "nodetree.hpp"
#include "dataMapOther.hpp"


class TreeClassOther;

class NodeOther: public Node {
	public:
		float *thetaOther;
		float *thDenOther;
		float *thNumOther;
		virtual bool makeDataforMerged(int a, int b, dataMapOther* w, TreeClassOther* tree, graphData* D);
		NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension);
		NodeOther(int nodeID, int parentID, bool isTerminal, int dimension);
		virtual ~NodeOther();
		
};

class TreeClassOther: public TreeClass {
	public:
		std::map<int, NodeOther*> nodeMap;	
		virtual int makeMergeNode(int a, int b);
		TreeClassOther(graphData* G, int dimension);
		~TreeClassOther();
};

TreeClassOther::~TreeClassOther() {
	// nothing to do
	// the base class destructor should be called after this
}

int TreeClassOther::makeMergeNode(int a, int b) {
	int c = this->numNodes;
	this->numNodes = c + 1;
	NodeOther* pnode = new NodeOther(c, -1, 0, dim);
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
	for (int i=0; i<D[0].numV; i++) {
		this->nodeMap[i] = new NodeOther(i, -1, 1, i, this->dim);
		this->topLevel.insert(i);
	}
	this->numNodes = D[0].numV;
}

NodeOther::NodeOther(int nodeID, int parentID, bool isTerminal, int dimension) :Node(nodeID, parentID, isTerminal, dimension) {
	this->thetaOther = new float[dimension];
	this->thDenOther = new float[dimension];
	this->thNumOther = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->thetaOther[d] = 0;
		this->thDenOther[d] = 0;
		this->thNumOther[d] = 0;
	}
}

NodeOther::NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension): Node(nodeID, parentID, isTerminal, vertID, dimension) {
	if (! isTerminal) {
		std::cerr << "bad call to NodeOther constructor!, given vertex ID for non-terminal node!\n";
		throw(1);
	}
	this->theta = new float[dimension];
	this->thDen = new float[dimension];
	this->thNum = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->theta[d] = 0;
		this->thDen[d] = 0;
		this->thNum[d] = 0;
	}
}






NodeOther::~NodeOther() {
	delete[] thetaOther;
	delete[] thNumOther;
	delete[] thDenOther;
}

bool NodeOther::makeDataforMerged(int a, int b, dataMapOther* w, TreeClassOther* tree, graphData* D) {
	float wc, wab;
	float x;
	for (int d=0;d<(tree->dim);d++) {
		wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
		assert(w[d].AddPair(this->nid,this->nid,wc));
		wc = w[d].get_uvOther(a,a) + w[d].get_uvOther(b,b) + w[d].get_uvOther(a,b);
		assert(w[d].AddPairOther(this->nid,this->nid,wc));
		w[d].degrees[this->nid] = w[d].degrees[a] + w[d].degrees[b];
		w[d].selfMissing[this->nid] = w[d].selfMissing[a] + w[d].selfMissing[b] + (w[d].degrees[a] * w[d].degrees[b]) - w[d].get_uv(a,b);
		w[d].nV[this->nid] = w[d].nV[a] + w[d].nV[b];
		w[d].oDegrees[this->nid] = w[d].oDegrees[a] + w[d].oDegrees[b];
		w[d].oSelfMissing[this->nid] = w[d].oSelfMissing[a] + w[d].oSelfMissing[b] + (w[d].oDegrees[a] * w[d].oDegrees[b]) - w[d].get_uv(a,b);
		w[d].oNV[this->nid] = w[d].oNV[a] + w[d].oNV[b];
		if (D[d].gtype=='w') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].degrees[a] * w[d].degrees[b];
			wab = w[d].get_uvOther(a,b);
			this->thNumOther[d] = wab;
			this->thDenOther[d] = w[d].oDegrees[a] * w[d].oDegrees[b];
		} else if (D[d].gtype=='b') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].nV[a] * w[d].nV[b];
			wab = w[d].get_uvOther(a,b);
			this->thNumOther[d] = wab;
			this->thDenOther[d] = w[d].oNV[a] * w[d].oNV[b];
		} else {
			std::cerr << "graph type "<<D[d].gtype<<" not yet supported (during theta calculation).\n";
			throw 1;
		}
		if (this->collapsed) {
			this->thNum[d] += tree->nodeMap[a]->thNum[d] + tree->nodeMap[b]->thNum[d];
			this->thDen[d] += tree->nodeMap[a]->thDen[d] + tree->nodeMap[b]->thDen[d];
			this->thNumOther[d] += tree->nodeMap[a]->thNumOther[d] + tree->nodeMap[b]->thNumOther[d];
			this->thDenOther[d] += tree->nodeMap[a]->thDenOther[d] + tree->nodeMap[b]->thDenOther[d];
		}
		x = this->thNum[d] / this->thDen[d];
		if (std::isnan(x)) x = 0;
		this->theta[d] = x;
		x = this->thNumOther[d] / this->thDenOther[d];
		if (std::isnan(x)) x = 0;
		this->thetaOther[d] = x;
		//TODO:: delete the following, only for debugging
		if (theta<0) {
			std::cerr << "bad theta being written!\n";
		}
		if (thetaOther<0) {
			std::cerr << "bad theta being written!\n";
		}
	}
	return 1;
}







#endif
