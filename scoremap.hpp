#ifndef SCOREMAP_HPP
#define SCOREMAP_HPP
#include <map>
#include <list>
#include <set>
#include "graphData.hpp"
#include "nodetree.hpp"



#define BIGNEG -9e300

class scoremap{
	public:
		typedef struct {
			double bestP;
			int bestK;
			std::map<int, double> scoreDest;
		} smap;
		typedef struct {
			int u;
			int v;
			float s;
		} pairScore;
		int bestK;
		double bestP;
		std::map<int, smap> scores;
		bool has_uv(int u, int v); //done
		pairScore getBestScore(); //done
		pairScore popBestScore(); //done
		bool hasPos(); //done
		bool isempty(); //done
		int AddPair(int u, int v, double score); //done
		bool eraseAll(); //done
		bool erase(int u,int v); //done
};

bool scoremap::erase(int u, int v) {
	std::map<int, scoremap::smap>::iterator it1;
	std::map<int, double>::iterator it2;
	if (bestK==u) {
		if (scores[u].bestK==v) {
			scores[u].scoreDest.erase(v);
			scores[u].bestP = BIGNEG;
			scores[u].bestK = -1;
			for (it2=scores[u].scoreDest.begin(); it2 != scores[u].scoreDest.end(); ++it2) {
				if ((*it2).second>scores[u].bestP) {
					scores[u].bestP = (*it2).second;
					scores[u].bestK = (*it2).first;
				}
			}
			bestP = BIGNEG;
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
			scores[u].bestP = BIGNEG;
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
		bestP = BIGNEG;
		bestK=-1;
	}	
};

bool scoremap::eraseAll() {
	if (scores.empty())
		return 0;
	else {
		scores.erase(scores.begin(), scores.end());
		bestP = BIGNEG;
		bestK = -1;
		return 1;
	}
};

int scoremap::AddPair(int u, int v, double score) {
	// to return 0 if the pair already existed, 1 if it didnt exist and is not the best, and 2 if it didnt exist and is the best
	if (scores.find(u)!=scores.end()) {
		if (scores[u].scoreDest.find(v)!=scores[u].scoreDest.end()) {
			scores[u].scoreDest[v] = score;
			if (score>scores[u].bestP) {
				scores[u].bestP = score; scores[u].bestK = v;
			}
			return 0;
		} else {
			scores[u].scoreDest[v] = score;
		}
	} else {
		smap stnew;
		stnew.bestP = score; stnew.bestK = v; stnew.scoreDest[v]=score;
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
	if ((not scores.empty())and(bestP>0))
		return 1;
	else
		return 0;
};

bool scoremap::has_uv(int u, int v) {
	std::map<int, scoremap::smap>::iterator it1;
	std::map<int, double>::iterator it2;
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
	double nb; int nk,i;
	pairScore ans;
	ans.u = bestK;
	ans.v = scores[bestK].bestK;
	ans.s = scores[bestK].scoreDest[ans.v];
	std::map<int, scoremap::smap>::iterator it1;
	std::map<int, double>::iterator it2;

	//now let us erase the last best
	scores[bestK].scoreDest.erase(ans.v);
	nb = BIGNEG;
	nk = -1;
	for (it2 = scores[bestK].scoreDest.begin(); it2 != scores[bestK].scoreDest.end(); ++it2) {
		if (nb>(*it2).second) {
			nb = (*it2).second;
			nk = (*it2).first;
		}
	}
	if (scores[bestK].scoreDest.empty()) {
		scores.erase(bestK);
	} else {
		scores[bestK].bestK = nk; scores[bestK].bestP = nb;
	}
	nb = BIGNEG;
	nk = -1;
	for (it1 = scores.begin(); it1 != scores.end(); ++it1) {
		if (nb>(*it1).second.bestP) {
			nb = (*it1).second.bestP;
			nk = (*it1).first;
		}
	}
	bestK = nk; bestP = nb;
	return ans;
};






#endif
