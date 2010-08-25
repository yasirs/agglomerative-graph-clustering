#ifndef LABELDATA_HPP
#define LABELDATA_HPP
#include "graphData.hpp"
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <tr1/unordered_map>

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiter);



class labelData{
	public:
		typedef std::tr1::unordered_map<int, float*> destList;
		std::map<int, std::string> int2Name;
		std::map<int, destList*> edgeList;
		int levels;
		//float get_uv(int u, int v);
		//void set_uv(int u, int v, float w);
		//bool has_uv(int u, int v);
		//bool delete_uv(int u, int v);
		//int Add_uv(int u, int v, float w); 
		labelData(graphData* Gorginal, const char* filename, int reg);
		void putSoFar(graphData* GsoFar, int l);
		void write(const char* filename);
		~labelData() {
			std::map<int, destList*>::iterator it;
			destList::iterator Init;
			for (it=edgeList.begin(); it!=edgeList.end(); ++it) {
				for (Init = it->second->begin(); Init != it->second->end(); ++Init) {
					delete[] Init->second;
				}
				(*it).second->clear();
				delete (*it).second;
			}
			edgeList.clear();
		}
};

void labelData::putSoFar(graphData* GsoFar, int levthis) {
	int u,v;
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = outIt->second->begin(); inIt != outIt->second->end(); ++inIt) {
			v = inIt->first;
			inIt->second[levthis] = GsoFar->get_uv(u,v);
		}
	}
};

void labelData::write(const char* fn) {
	int u,v;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	std::string uname, vname;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		uname = int2Name[u];
		for (inIt = outIt->second->begin(); inIt != outIt->second->end(); ++inIt) {
			v = inIt->first;
			vname = int2Name[v];
			file << uname << '\t' << vname;
			for (int i=0; i <= levels; i++) {
				file << '\t' << inIt->second[i];
			}
			file << '\n';
		}
	}
	file.close();
};



labelData::labelData(graphData* Goriginal, const char* filename, int reg) {
	this->levels = reg;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	int linesread=0;
	std::string St1, St2;
	file.open(filename,std::ios::in);
	if (! file.is_open()) throw( 0);
	this->int2Name = Goriginal->int2Name;
	while (! file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			linesread++;
			if (Goriginal->name2Int.find(tok[0])==Goriginal->name2Int.end()) {
				// do nothing
				continue;
				throw( 0);
			} else {
				u = Goriginal->name2Int[tok[0]];
			}
			if (Goriginal->name2Int.find(tok[1])==Goriginal->name2Int.end()) {
				// do nothing
				continue;
				throw( 0);
			} else {
				v = Goriginal->name2Int[tok[1]];
			}
			std::stringstream temp;
			temp << tok[2];
			temp >> weight;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			float* fp = new float[levels+1];
			fp[0] = weight;
			(*edgeList[u])[v] = fp;
		}
	}
};








#endif

