#include <stack>

class NCNZGenerator{
	public:
		NCNZGenerator(TreeClass* t,int mydim);
		bool isDone();
		int goNext();
	private:
		std::stack<int> toVisit;
		TreeClass* tree;
		int d;
};

NCNZGenerator::NCNZGenerator(TreeClass* t, int mydim) {
	this->tree = t;
	for (std::set<int>::iterator tit(tree->topLevel.begin()); tit != tree->topLevel.end(); ++tit) {
		this->toVisit.push(*tit);
	}
	this->d = mydim;
}

bool NCNZGenerator::isDone() {
	return (this->toVisit.empty());
}

int NCNZGenerator::goNext() {
	if (this->toVisit.empty()) return 0;
	int curNode = this->toVisit.pop();
	Node* pnode = this->tree->nodeMap[curNode];
	while (1) {
		if (! pnode->collapsed) {
			for (std::set<int>::iterator cit(pnode->childSet.begin()); cit != pnode->childSet.end(); ++cit) {
				toVisit.push(*cit);
			}
		}
		if (pnode->theta[d] != 0)
			return curNode;
		else {
			if (this->toVisit.empty()) {
				return 0;
			} else {
				curNode = this->toVisit.pop();
				pnode = this->tree->nodeMap[curNode];
			}
		}
	}
}

		



class AllChildVertGenerator{
	public:
		AllChildVertGenerator(TreeClass* t, int nodenum);
		bool isDone();
		int goNext();
	private:
		std::stack<int> toVisit;
		TreeClass* tree;
};

bool AllChildVertGenerator::isDone() {
	return toVisit.empty();
}


AllChildVertGenerator::AllChildVertGenerator(TreeClass* t, int nodenum) {
	this->tree = t;
	toVisit.push(nodenum);
}

int AllChildVertGenerator::goNext() {
	int curNode;
	Node* pnode;
	if (toVisit.empty()) {
		return 0;
	}
	curNode = toVisit.pop();
	pnode = tree->nodeMap[curNode];
	while (1) {
		if (pnode->isTerm) {
			// nothing to do but return this
			return curNode;
		} else {
			if (pnode->vertsComputed) {
				std::set<int>::iterator vit (pnode->vertexSet.begin());
				for (; vit != pnode->vertexSet.end(); vit++) {
					toVisit.push(*vit);
				}
				curNode = toVisit.pop();
				assert( tree->nodeMap[curNode]->isTerm);
				return curNode;
			} else {
				// curNode is not a terminal, neither are its vertices computed
				std::set<int>::iterator nit (pnode->childSet.begin());
				for (; nit != pnode->vertexSet.end(); nit++) {
					toVisit.push(*nit);
				}
				// pop and continue
				curNode = toVisit.pop();
				pnode = tree->nodeMap[curNode];
			}
		}
	}



			
