#ifndef SCOREMAP_HPP
#define SCOREMAP_HPP
#include <map>
#include <list>
#include <set>
#include <cmath>
#include "graphData.hpp"
#include "nodetree.hpp"



#define BIGNEG -9e300
#define EPS 1e-6

class scoremap{
	public:
		typedef struct twoScorestruct{
			double joinScore;
			double centerMscore;
			twoScorestruct (double j, double c) {
				joinScore = j;
				centerMscore = c;
			}
			twoScorestruct() {}
		} twoScores;
		typedef struct {
			twoScores bestP;
			int bestK;
			std::map<int, twoScores> scoreDest;
		} smap;
		typedef struct {
			int u;
			int v;
			twoScores s;
		} pairScore;
		int bestK;
		twoScores bestP;
		std::map<int, smap> scores;
		bool has_uv(int u, int v); //done
		pairScore getBestScore(); //done
		pairScore popBestScore(); //done
		bool hasPos(); //done
		bool isempty(); //done
		int AddPair(int u, int v, double jscore, double cscore); //done
		int AddPair(int u, int v, twoScores score); 
		int AddTo(int u, int v, double scorej); // done
		bool eraseAll(); //done
		bool erase(int u,int v); //done
		int has_u(int u);
};

int scoremap::has_u(int u) {
	bool ans = 0;
	std::map<int, smap>::iterator OutIt;
	std::map<int, twoScores>::iterator InIt;
	if (scores.find(u)!=scores.end()) return (*(scores[u].scoreDest.begin())).first;
	for (OutIt = scores.begin(); OutIt != scores.end(); ++OutIt) {
		for (InIt = (*OutIt).second.scoreDest.begin();InIt != (*OutIt).second.scoreDest.end(); ++InIt) {
			if ((*InIt).first==u) return (*OutIt).first;
		}
	}
	return -1;
}

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
	std::map<int, twoScores>::iterator it2;
	if (bestK==u) {
		if (scores[u].bestK==v) {
			scores[u].scoreDest.erase(v);
			scores[u].bestP = twoScorestruct(BIGNEG,BIGNEG);
			scores[u].bestK = -1;
			for (it2=scores[u].scoreDest.begin(); it2 != scores[u].scoreDest.end(); ++it2) {
				if ((*it2).second>scores[u].bestP) {
					scores[u].bestP = (*it2).second;
					scores[u].bestK = (*it2).first;
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
			scores[u].scoreDest.erase(v);
		}
	} else {
		if (scores[u].bestK==v) {
			scores[u].scoreDest.erase(v);
			scores[u].bestP =  twoScorestruct(BIGNEG,BIGNEG);
			scores[u].bestK = -1;
			for (it2=scores[u].scoreDest.begin(); it2 != scores[u].scoreDest.end(); ++it2) {
				if ((*it2).second>scores[u].bestP) {
					scores[u].bestP = (*it2).second;
					scores[u].bestK = (*it2).first;
				}
			}
		} else {
			scores[u].scoreDest.erase(v);
		}
	}
	if (scores[u].scoreDest.empty()) {
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
	}
};


int scoremap::AddTo(int u, int v, double scorej) {
	double nsj;
	double nsc;
	nsj = scores[u].scoreDest[v].joinScore + scorej;
	nsc = scores[u].scoreDest[v].joinScore;
	assert( not(this->AddPair(u,v,nsj,nsc)));
}



int scoremap::AddPair(int u, int v, double scorej, double scorec) {
	if (u==v) throw 1; //TODO this was for debugging, possibly leave them in
	twoScores score = twoScorestruct(scorej, scorec);
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
	}
};

bool scoremap::isempty() {
	if (scores.empty()) 
		return 1;
	else
		return 0;
};

bool scoremap::hasPos() {
	if ((not scores.empty())and(bestP>twoScorestruct(0,0)))
		return 1;
	else
		return 0;
};

bool scoremap::has_uv(int u, int v) {
	std::map<int, scoremap::smap>::iterator it1;
	std::map<int, twoScores>::iterator it2;
	it1 = scores.find(u);
	if (it1 != scores.end()) {
		it2 = (*it1).second.scoreDest.find(v);
		if (it2 != (*it1).second.scoreDest.end())
			return 1;
		else
			return 0;
	}
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
	twoScores nb; int nk,i;
	pairScore ans;
	ans.u = bestK;
	ans.v = scores[bestK].bestK;
	ans.s = scores[bestK].scoreDest[ans.v];
	std::map<int, scoremap::smap>::iterator it1;
	std::map<int, twoScores>::iterator it2;

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
