#ifndef GRAPHDATA_HPP
#define GRAPHDATA_HPP
#include <map>


class graphData{
	public:
		int numV;
		char gtype;
		typedef std::map<int, float> destList;
		std::map<int, destList> edgeList;
		bool readWeighted(const char* filename);
		bool readBinary(const char* filename);
		int degree(int i);
		std::set<int> neighbors(int i);
};


#endif

