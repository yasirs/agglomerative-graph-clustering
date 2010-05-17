#include "nodetree.hpp"
#include "graphData.hpp"


graphData* getResidual(graphData* D, TreeClass* tree, const int dim) {
	int d;
	int t,p,n,i;
	std::stack<int> st;
	graphData Dnew = new graphData[dim];
	std::set<int>::iterator intit, intit2, intit3;
	for (d=0;d<dim;d++) {
		Dnew[d].int2Name = D[d].int2Name;
		Dnew[d].name2Int = D[d].name2Int;
	}
	// let us go through the top level first
	for (intit = tree->topLevel.begin(); intit != tree->topLevel.end(); intit++) {
		// top level
		t = (*intit);
		// let us compute the vertex sets via tree search
		st.push(t);
		while (st.size()>0) {
			if (! tree->nodeMap[t].collapsed) {
				// we need to compute its vertex set
				for (intit2 = tree->nodeMap[t]->childSet.begin(); intit2 != tree->nodeMap[t]->childSet.end(); ++intit2) {
					

