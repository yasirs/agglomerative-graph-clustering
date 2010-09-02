#ifndef NODETREEOTHER_HPP
#define NODETREEOTHER_HPP
#include "nodetree.hpp"
#include "dataMapOther.hpp"


class TreeClassOther;

class NodeOther: public Node {
	public:
		float *thetaOriginal;
		float *thDenOriginal;
		float *thNumOriginal;
		float *thetaSoFar;
		float *thDenSoFar;
		float *thNumSoFar;
		virtual bool makeThetaforMerged(int a, int b, dataMap* w, TreeClass* tree, graphData* D);
		NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension);
		NodeOther(int nodeID, int parentID, bool isTerminal, int dimension);
		virtual ~NodeOther();
		
};

class TreeClassOther: public TreeClass {
	public:
		//std::map<int, NodeOther*> nodeMap;	
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
	this->thetaOriginal = new float[dimension];
	this->thDenOriginal = new float[dimension];
	this->thNumOriginal = new float[dimension];
	this->thetaSoFar = new float[dimension];
	this->thDenSoFar = new float[dimension];
	this->thNumSoFar = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->thetaOriginal[d] = 0;
		this->thDenOriginal[d] = 0;
		this->thNumOriginal[d] = 0;
		this->thetaSoFar[d] = 0;
		this->thDenSoFar[d] = 0;
		this->thNumSoFar[d] = 0;
	}
}

NodeOther::NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension): Node(nodeID, parentID, isTerminal, vertID, dimension) {
	if (! isTerminal) {
		std::cerr << "bad call to NodeOther constructor!, given vertex ID for non-terminal node!\n";
		throw(1);
	}
	this->thetaOriginal = new float[dimension];
	this->thDenOriginal = new float[dimension];
	this->thNumOriginal = new float[dimension];
	this->thetaSoFar = new float[dimension];
	this->thDenSoFar = new float[dimension];
	this->thNumSoFar = new float[dimension];
	for (int d=0;d<dimension; d++) {
		this->thetaOriginal[d] = 0;
		this->thDenOriginal[d] = 0;
		this->thNumOriginal[d] = 0;
		this->thetaSoFar[d] = 0;
		this->thDenSoFar[d] = 0;
		this->thNumSoFar[d] = 0;
	}
}






NodeOther::~NodeOther() {
	delete[] thetaSoFar;
	delete[] thNumSoFar;
	delete[] thDenSoFar;
	delete[] thetaOriginal;
	delete[] thNumOriginal;
	delete[] thDenOriginal;
}

bool NodeOther::makeThetaforMerged(int a, int b, dataMap* ww, TreeClass* tree, graphData* D) {
	float wc, wab;
	float x;
	dataMapOther *w = (dataMapOther*) ww;

	for (int d=0;d<(tree->dim);d++) {
		//wc = w[d].get_uv(a,a) + w[d].get_uv(b,b) + w[d].get_uv(a,b);
		//assert(w[d].AddPair(this->nid,this->nid,wc));
		//wc = w[d].get_uvOriginal(a,a) + w[d].get_uvOriginal(b,b) + w[d].get_uvOriginal(a,b);
		//assert(w[d].AddPairOriginal(this->nid,this->nid,wc));
		//wc = w[d].get_uvSoFar(a,a) + w[d].get_uvSoFar(b,b) + w[d].get_uvSoFar(a,b);
		//assert(w[d].AddPairSoFar(this->nid,this->nid,wc));
		//w[d].degrees[this->nid] = w[d].degrees[a] + w[d].degrees[b];
		//w[d].selfMissing[this->nid] = w[d].selfMissing[a] + w[d].selfMissing[b] + (w[d].degrees[a] * w[d].degrees[b]) - w[d].get_uv(a,b);
		//w[d].nV[this->nid] = w[d].nV[a] + w[d].nV[b];
		//w[d].nV.push_back(w[d].nV[a]+w[d].nV[b]);
		/*w[d].oDegrees[this->nid] = w[d].oDegrees[a] + w[d].oDegrees[b];
		w[d].oSelfMissing[this->nid] = w[d].oSelfMissing[a] + w[d].oSelfMissing[b] + (w[d].oDegrees[a] * w[d].oDegrees[b]) - w[d].get_uvOriginal(a,b);
		w[d].oNV[this->nid] = w[d].oNV[a] + w[d].oNV[b];
		w[d].sDegrees[this->nid] = w[d].sDegrees[a] + w[d].sDegrees[b];
		w[d].sSelfMissing[this->nid] = w[d].sSelfMissing[a] + w[d].sSelfMissing[b] + (w[d].sDegrees[a] * w[d].sDegrees[b]) - w[d].get_uvSoFar(a,b);
		w[d].sNV[this->nid] = w[d].sNV[a] + w[d].sNV[b];*/
		if (D[d].gtype=='w') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].degrees[a] * w[d].degrees[b];
			wab = w[d].get_uvSoFar(a,b);
			this->thNumSoFar[d] = wab;
			this->thDenSoFar[d] = w[d].sDegrees[a] * w[d].sDegrees[b];
			wab = w[d].get_uvOriginal(a,b);
			this->thNumOriginal[d] = wab;
			this->thDenOriginal[d] = w[d].oDegrees[a] * w[d].oDegrees[b];
		} else if (D[d].gtype=='b') {
			wab = w[d].get_uv(a,b);
			this->thNum[d] = wab;
			this->thDen[d] = w[d].nV[a] * w[d].nV[b];
			wab = w[d].get_uvOriginal(a,b);
			this->thNumOriginal[d] = wab;
			this->thDenOriginal[d] = w[d].oNV[a] * w[d].oNV[b];
			wab = w[d].get_uvSoFar(a,b);
			this->thNumSoFar[d] = wab;
			this->thDenSoFar[d] = w[d].sNV[a] * w[d].sNV[b];
		} else {
			std::cerr << "graph type "<<D[d].gtype<<" not yet supported (during theta calculation).\n";
			throw 1;
		}
		if (this->collapsed) {
			this->thNum[d] += ( (NodeOther*) tree->nodeMap[a])->thNum[d] + ( (NodeOther*) tree->nodeMap[b])->thNum[d];
			this->thDen[d] += ( (NodeOther*) tree->nodeMap[a])->thDen[d] + ( (NodeOther*) tree->nodeMap[b])->thDen[d];
			this->thNumOriginal[d] += ( (NodeOther*) tree->nodeMap[a])->thNumOriginal[d] + ( (NodeOther*) tree->nodeMap[b])->thNumOriginal[d];
			this->thDenOriginal[d] += ( (NodeOther*) tree->nodeMap[a])->thDenOriginal[d] + ( (NodeOther*) tree->nodeMap[b])->thDenOriginal[d];
			this->thNumSoFar[d] += ( (NodeOther*) tree->nodeMap[a])->thNumSoFar[d] +( (NodeOther*) tree->nodeMap[b])->thNumSoFar[d];
			this->thDenSoFar[d] +=( (NodeOther*) tree->nodeMap[a])->thDenSoFar[d] + ( (NodeOther*) tree->nodeMap[b])->thDenSoFar[d];
		}
		/*
		x = this->thNum[d] / this->thDen[d];
		if (std::isnan(x)) x = 0;
		this->theta[d] = x;
		*/
		x = this->thNumOriginal[d] / this->thDenOriginal[d];
		if (std::isnan(x)) x = 0;
		this->thetaOriginal[d] = x;
		x = this->thNumSoFar[d] / this->thDenSoFar[d];
		if (std::isnan(x)) x = 0;
		this->thetaSoFar[d] = x;
		//TODO:: delete the following, only for debugging
		if (this->thetaOriginal[d]<0) {
			std::cerr << "bad theta_Original being written!\n";
		}
		if (this->thetaSoFar[d]<0) {
			std::cerr << "bad theta_SoFar being written!\n";
		}
		this->theta[d] = std::max(0.0f, (this->thetaOriginal[d] - this->thetaSoFar[d])/(1 - this->thetaSoFar[d]));
	}
	return 1;
}







#endif
