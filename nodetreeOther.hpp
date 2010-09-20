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
		/*float *thetaOriginal;
		float *thDenOriginal;
		float *thNumOriginal;
		float *thetaSoFar;
		float *thDenSoFar;
		float *thNumSoFar;*/
		virtual bool writeThetaforMerged(int a, int b, dataMap* w, TreeClass* tree, graphData* D);
		NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension, graphData* D);
		NodeOther(int nodeID, int parentID, bool isTerminal, int dimension, graphData* D);
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
	NodeOther* pnode = new NodeOther(c, -1, 0, dim, D);
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
		this->nodeMap[i] = new NodeOther(i, -1, 1, i, this->dim, D);
		this->topLevel.insert(i);
	}
	this->numNodes = D[0].numV;
}

NodeOther::NodeOther(int nodeID, int parentID, bool isTerminal, int dimension, graphData* D) :Node(nodeID, parentID, isTerminal, dimension, D) {
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
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
	}
	/*this->thetaOriginal = new float[dimension];
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
	*/
}

NodeOther::NodeOther(int nodeID, int parentID, bool isTerminal, int vertID, int dimension, graphData* D): Node(nodeID, parentID, isTerminal, vertID, dimension, D) {
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
		} else {
			std::cerr << "dont know what to do with graph type "<<D[d].gtype<<" while making node\n";
			throw(1);
		}
	}
	/*
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
	}*/
}






NodeOther::~NodeOther() {
	//delete[] params;
	delete[] paramsSoFar;
	delete[] paramsOriginal;
	/*delete[] thetaSoFar;
	delete[] thNumSoFar;
	delete[] thDenSoFar;
	delete[] thetaOriginal;
	delete[] thNumOriginal;
	delete[] thDenOriginal;
	*/
}

bool NodeOther::writeThetaforMerged(int a, int b, dataMap* ww, TreeClass* tree, graphData* D) {
	dataMapOther *w = (dataMapOther*) ww;

	for (int d=0;d<(tree->dim);d++) {
		this->params[d]->calculate(w[d].get_uv(a,b),w[d].datvert[a],w[d].datvert[b]);
		this->paramsOriginal[d]->calculate(w[d].get_uvOriginal(a,b),w[d].oDatvert[a],w[d].oDatvert[b]);
		this->paramsSoFar[d]->calculate(w[d].get_uvSoFar(a,b),w[d].sDatvert[a],w[d].sDatvert[b]);
		/*if (D[d].gtype=='w') {
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
		}*/
		if (this->collapsed) {
			this->params[d]->collapse(((NodeOther*) tree->nodeMap[a])->params[d], ((NodeOther*) tree->nodeMap[b])->params[d]);
			this->paramsOriginal[d]->collapse(((NodeOther*) tree->nodeMap[a])->paramsOriginal[d], ((NodeOther*) tree->nodeMap[b])->paramsOriginal[d]);
			this->paramsSoFar[d]->collapse(((NodeOther*) tree->nodeMap[a])->paramsSoFar[d], ((NodeOther*) tree->nodeMap[b])->paramsSoFar[d]);
			/*this->thNum[d] += ( (NodeOther*) tree->nodeMap[a])->thNum[d] + ( (NodeOther*) tree->nodeMap[b])->thNum[d];
			this->thDen[d] += ( (NodeOther*) tree->nodeMap[a])->thDen[d] + ( (NodeOther*) tree->nodeMap[b])->thDen[d];
			this->thNumOriginal[d] += ( (NodeOther*) tree->nodeMap[a])->thNumOriginal[d] + ( (NodeOther*) tree->nodeMap[b])->thNumOriginal[d];
			this->thDenOriginal[d] += ( (NodeOther*) tree->nodeMap[a])->thDenOriginal[d] + ( (NodeOther*) tree->nodeMap[b])->thDenOriginal[d];
			this->thNumSoFar[d] += ( (NodeOther*) tree->nodeMap[a])->thNumSoFar[d] +( (NodeOther*) tree->nodeMap[b])->thNumSoFar[d];
			this->thDenSoFar[d] +=( (NodeOther*) tree->nodeMap[a])->thDenSoFar[d] + ( (NodeOther*) tree->nodeMap[b])->thDenSoFar[d];
			*/
		}
		/*
		x = this->thNum[d] / this->thDen[d];
		if (std::isnan(x)) x = 0;
		this->theta[d] = x;
		*/
		this->params[d]->cleanup();
		this->paramsOriginal[d]->cleanup();
		this->paramsSoFar[d]->cleanup();
		/*x = this->thNumOriginal[d] / this->thDenOriginal[d];
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
		}*/
		this->params[d]->bestfromSoFar(this->paramsOriginal[d],this->paramsSoFar[d]);
		//this->theta[d] = std::max(0.0f, (this->thetaOriginal[d] - this->thetaSoFar[d])/(1 - this->thetaSoFar[d]));
	}
	return 1;
}







#endif
