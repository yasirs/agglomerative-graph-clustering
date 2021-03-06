#ifndef SCOREMAP_HPP
#define SCOREMAP_HPP
#include <map>
#include <list>
#include <set>
#include <cmath>
#include <cassert>
#include <iostream>
#include <tr1/unordered_map>
#include "graphData.hpp"
#include "nodetree.hpp"



#define BIGNEG -9e300
#define EPS 1e-6

class scoremap{
	public:
		typedef struct twoScorestruct{
			float joinScore;
			float centerMscore;
			twoScorestruct (float j, float c) {
				joinScore = j;
				centerMscore = c;
			}
			twoScorestruct() {}
		} twoScores;
		typedef std::tr1::unordered_map<int, twoScores> scDestType;
		typedef struct {
			twoScores bestP;
			int bestK;
			scDestType scoreDest;
		} smap;
		typedef struct {
			int u;
			int v;
			twoScores s;
		} pairScore;
		scoremap() {
			bestP = twoScorestruct(BIGNEG,BIGNEG);
			bestK = -1;
		}
		int bestK;
		twoScores bestP;
		std::map<int, smap> scores;
		bool has_uv(int u, int v); //done
		pairScore getBestScore(); //done
		pairScore popBestScore(); //done
		twoScores get_uv(int u, int v);
		smap get_u(int u) { return scores[u]; }
		bool hasPos(); //done
		bool isEmpty(); //done
		int AddPair(int u, int v, float jscore, float cscore); //done
		int AddPair(int u, int v, twoScores score); 
		bool AddTo(int u, int v, float scorej); // done
		bool eraseAll(); //done
		bool erase(int u,int v); //done
		int has_u(int u);
		int size();
		std::set<int>* allPartners_u(int u);
};

bool operator<(const scoremap::twoScores& s1, const scoremap::twoScores& s2);
bool operator>(const scoremap::twoScores& s1, const scoremap::twoScores& s2);



int scoremap::size() {
	return scores.size();
}

bool scoremap::AddTo(int u, int v, float scorej) {
	//TODO:: 
	//std::cout << "before:"<<scores[u].scoreDest[v].joinScore;
	std::map<int, smap>::iterator smOut = scores.find(u);
	scDestType::iterator smIn = (*smOut).second.scoreDest.find(v);
	if (smIn==(*smOut).second.scoreDest.end()) {
		return 0;
	} else {
		(*smIn).second.joinScore += scorej;
		twoScores ts1, ts2;
		if (((*smIn).second) > ((*smOut).second.bestP)) { // new bestP is this one!
			(*smOut).second.bestP = (*smIn).second;
			(*smOut).second.bestK = (*smIn).first;
			if ((*smIn).second > this->bestP) {
				this->bestP = (*smIn).second;
				this->bestK = (*smOut).first;
			}
		} else {
			if ( (*smOut).second.bestK == (*smIn).first ) { // this was the old bestP, is it still so?
				(*smOut).second.bestP = (*smIn).second;
				for (scDestType::iterator inIt=(*smOut).second.scoreDest.begin();
					inIt != (*smOut).second.scoreDest.end(); ++inIt ) {
					if ((*inIt).second > (*smOut).second.bestP) {
						(*smOut).second.bestP = (*inIt).second;
						(*smOut).second.bestK = (*inIt).first;
					}
				}
				if (this->bestK == (*smOut).first) {
					this->bestP = (*smOut).second.bestP;
					for (std::map<int, smap>::iterator outIt=scores.begin(); outIt != scores.end(); outIt++) {
						if ((*outIt).second.bestP > this->bestP) {
							this->bestP = (*outIt).second.bestP;
							this->bestK = (*outIt).first;
						}
					}
				}
			}
		}				
		return 1;
	}
};


scoremap::twoScores scoremap::get_uv(int u, int v) {
	if (has_uv(u,v)) {
		return scores[u].scoreDest[v];
	} else {
		std::cerr << "Error!! doesnt have "<<u<<", "<<v<<"\n";
		throw 1;
	}
};

bool scoremap::has_uv(int u, int v) {
	if (scores.find(u)==scores.end()) {
		return 0;
	} else if (scores[u].scoreDest.find(v)==scores[u].scoreDest.end()) {
		return 0;
	} else return 1;
};

std::set<int>* scoremap::allPartners_u(int u) {
	std::set<int>* ans;
	ans = new std::set<int>;
	scDestType::iterator it;
	if (scores.find(u)==scores.end()) {std::cout << "doesn't have primary key "<<u<<"\n"; return ans;}
	for (it=scores[u].scoreDest.begin();it!=scores[u].scoreDest.end(); ++it) {
		(*ans).insert((*it).first);
	}
	return ans;
};



int scoremap::has_u(int u) {
	std::map<int, smap>::iterator OutIt;
	scDestType::iterator InIt;
	if (scores.find(u)!=scores.end()) return (*(scores[u].scoreDest.begin())).first;
	for (OutIt = scores.begin(); OutIt != scores.end(); ++OutIt) {
		for (InIt = (*OutIt).second.scoreDest.begin();InIt != (*OutIt).second.scoreDest.end(); ++InIt) {
			if ((*InIt).first==u) return (*OutIt).first;
		}
	}
	return -1;
};

bool operator<(const scoremap::twoScores& s1, const scoremap::twoScores& s2) {
	if (fabs(s1.joinScore-s2.joinScore)<EPS) {
		return (s1.centerMscore<s2.centerMscore);
	} else {
		return (s1.joinScore<s2.joinScore);
	}
};

bool operator>(const scoremap::twoScores& s1, const scoremap::twoScores& s2) {
	if (fabs(s1.joinScore-s2.joinScore)<EPS) {
		return (s1.centerMscore>s2.centerMscore);
	} else {
		return (s1.joinScore>s2.joinScore);
	}
};


bool scoremap::erase(int u, int v) {
	std::map<int, scoremap::smap>::iterator it1;
	scDestType::iterator it2;
	std::map<int, scoremap::smap>::iterator outIt(scores.find(u));
	//TODO
	if (outIt==scores.end())
		throw 1;
	if (bestK==u) {
		if ((*outIt).second.bestK==v) {
			(*outIt).second.scoreDest.erase(v);
			(*outIt).second.bestP = twoScorestruct(BIGNEG,BIGNEG);
			(*outIt).second.bestK = -1;
			for (it2=(*outIt).second.scoreDest.begin(); it2 != (*outIt).second.scoreDest.end(); ++it2) {
				if ((*it2).second>(*outIt).second.bestP) {
					(*outIt).second.bestP = (*it2).second;
					(*outIt).second.bestK = (*it2).first;
				}
			}
			bestP = twoScorestruct(BIGNEG,BIGNEG);
			bestK = -1;
			for (it1=scores.begin(); it1 != scores.end(); ++it1) {
				if ((*it1).second.bestP > bestP) {
					bestP = (*it1).second.bestP;
					bestK = (*it1).first;
				}
			}
		} else {
			(*outIt).second.scoreDest.erase(v);
		}
	} else {
		if ((*outIt).second.bestK==v) {
			(*outIt).second.scoreDest.erase(v);
			(*outIt).second.bestP =  twoScorestruct(BIGNEG,BIGNEG);
			(*outIt).second.bestK = -1;
			for (it2=(*outIt).second.scoreDest.begin(); it2 != (*outIt).second.scoreDest.end(); ++it2) {
				if ((*it2).second>(*outIt).second.bestP) {
					(*outIt).second.bestP = (*it2).second;
					(*outIt).second.bestK = (*it2).first;
				}
			}
		} else {
			(*outIt).second.scoreDest.erase(v);
		}
	}
	if ((*outIt).second.scoreDest.empty()) {
		scores.erase(u);
	}
	if (scores.empty()) {
		bestP =  twoScorestruct(BIGNEG,BIGNEG);
		bestK=-1;
	}
	return 1;
};

bool scoremap::eraseAll() {
	if (scores.empty())
		return 0;
	else {
		scores.erase(scores.begin(), scores.end());
		bestP =  twoScorestruct(BIGNEG,BIGNEG);
		bestK = -1;
		return 1;
	}
};

int scoremap::AddPair(int u, int v, twoScores score) {
	// to return 0 if the pair already existed, 1 if it didnt exist and is not the best, and 2 if it didnt exist and is the best
	if (scores.find(u)!=scores.end()) {
		if (scores[u].scoreDest.find(v)!=scores[u].scoreDest.end()) {
			scores[u].scoreDest[v] =  score;
			if (score>scores[u].bestP) {
				scores[u].bestP = score; scores[u].bestK = v;
			}
			return 0;
		} else {
			scores[u].scoreDest[v] = score;
		}
	} else {
		smap stnew;
		stnew.bestP = twoScorestruct(BIGNEG,BIGNEG); stnew.bestK = -1; stnew.scoreDest[v]=score;
		scores[u] = stnew;
	}
	if (score>scores[u].bestP) {
		scores[u].bestP = score; scores[u].bestK = v;
		if (score>bestP) {
			bestP = score; bestK = u;
			return 2;
		} else
			return 1;
	} else {
		return 1;
	}
};






int scoremap::AddPair(int u, int v, float scorej, float scorec) {
	if (u==v) throw 1; //TODO this was for debugging, possibly leave them in
	twoScores score = twoScorestruct(scorej, scorec);
	return (this->AddPair(u,v,score));
};

bool scoremap::isEmpty() {
	if (scores.empty()) 
		return 1;
	else
		return 0;
};

bool scoremap::hasPos() {
	if ((! scores.empty())&&(bestP>twoScorestruct(0,0)))
		return 1;
	else
		return 0;
};


scoremap::pairScore scoremap::getBestScore() {
	pairScore ans;
	ans.u = bestK;
	ans.v = scores[bestK].bestK;
	ans.s = scores[bestK].scoreDest[ans.v];
	return ans;
}

scoremap::pairScore scoremap::popBestScore() {
	twoScores nb; int nk;
	pairScore ans;
	ans.u = bestK;
	ans.v = scores[bestK].bestK;
	ans.s = scores[bestK].scoreDest[ans.v];
	std::map<int, scoremap::smap>::iterator it1;
	scDestType::iterator it2;

	//now let us erase the last best
	scores[bestK].scoreDest.erase(ans.v);
	nb =  twoScorestruct(BIGNEG,BIGNEG);
	nk = -1;
	for (it2 = scores[bestK].scoreDest.begin(); it2 != scores[bestK].scoreDest.end(); ++it2) {
		if (nb<(*it2).second) {
			nb = (*it2).second;
			nk = (*it2).first;
		}
	}
	if (scores[bestK].scoreDest.empty()) {
		scores.erase(bestK);
	} else {
		scores[bestK].bestK = nk; scores[bestK].bestP = nb;
	}
	nb =  twoScorestruct(BIGNEG,BIGNEG);
	nk = -1;
	for (it1 = scores.begin(); it1 != scores.end(); ++it1) {
		if (nb<(*it1).second.bestP) {
			nb = (*it1).second.bestP;
			nk = (*it1).first;
		}
	}
	bestK = nk; bestP = nb;
	return ans;
};






#endif
