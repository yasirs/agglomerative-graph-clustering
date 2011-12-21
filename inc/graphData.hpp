#ifndef GRAPHDATA_HPP
#define GRAPHDATA_HPP
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <exception>
#include <tr1/unordered_map>
#include <boost/lexical_cast.hpp>

#ifndef DEBUGMODE
#define DEBUGMODE 0
#endif


void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiter);

std::string getstvec(std::vector<std::string> a, int i) {
	return a[i];
}

void putstvec(std::vector<std::string> a, int i, const char* st) {
	a[i] = st;
}

class graphData{
	public:
		unsigned int numV,multiPartite;
		char gtype;
		std::string getType() { return std::string(1,gtype); }
		bool setType(std::string ss) {
			if (ss.length()!=1) {
				throw(std::runtime_error("type should be of length 1, not" + boost::lexical_cast<std::string>(ss.length())));
				if (ss.length()>1) gtype = ss[0];
				return false;
			}
			gtype = ss[0];
			return true;
		}
		float aveP;
		float Etot;
		typedef std::tr1::unordered_map<int, float> destList;
		std::map<int, std::string> int2Name;
		std::map<std::string, int> name2Int;
		std::map<int, destList*> edgeList;
		std::map<int, int> typeList;
		float get_uv(int u, int v);
		void set_uv(int u, int v, float w);
		bool has_uv(int u, int v);
		bool delete_uv(int u, int v);
		int Add_uv(int u, int v, float w); 
		bool readBinaryBasedOnOld(graphData* Gorginal, const char* filename);
		bool readGeneralBasedOnOld(graphData* Gorginal, const char* filename, int nsoFar, const char* sep ="\t");
		bool readWeighted(const char* filename);
		bool readBinary(const char* filename);
		bool readGeneral(const char* filename, const char* sep="\t"); //something wrong here, doesnt work??
		void writeBoth(const char* filename);
		void writeSingle(const char* filename);
		void writeSingle_noname(const char* filename);
		int degree(int i);
		int get_int_from_name(const char * vname);
		void print_name_from_int(int i);
		//std::set<int> neighbors(int i);
		~graphData() {
			std::map<int, destList*>::iterator it;
			for (it=edgeList.begin(); it!=edgeList.end(); ++it) {
				(*it).second->clear();
				delete (*it).second;
			}
			edgeList.clear();
		}
		void copyNoEdges(graphData& Dnew);
		bool hasName(const std::string& name);
		bool hasName(const char* name);
};

bool graphData::hasName(const std::string& name) {
	return (name2Int.find(name)!=name2Int.end());
}

bool graphData::hasName(const char* name) {
	return (name2Int.find(std::string(name))!=name2Int.end());
}

void graphData::print_name_from_int(int i) {
	std::cout << int2Name[i] << "\n";
}

int graphData::get_int_from_name(const char* vname) {
	if (name2Int.find(std::string(vname))==name2Int.end()) {
		std::cout << "Error:: not found vertex named " << vname << "in the graph!\n";
		return 0;
	} else {
		return name2Int[std::string(vname)];
	}
}

void graphData::copyNoEdges(graphData& Dnew) {
	Dnew.int2Name = this->int2Name;
	Dnew.name2Int = this->name2Int;
	Dnew.gtype = this->gtype;
	Dnew.numV = this->numV;
}


int graphData::degree(int i) {
	return edgeList[i]->size();
};




void graphData::writeSingle_noname(const char* fn) {
	int u,v;
	float weight;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = edgeList[u]->begin(); inIt != edgeList[u]->end(); ++inIt) {
			v = inIt->first;
			weight = inIt->second;
			if (u<v) {
				if (gtype=='w') {
					file << u <<'\t' << v  << '\t' << weight << '\n';
				}
				else if (gtype=='b') {
					file << u <<'\t' <<  v << '\n';
				}
				else if (gtype=='p') {
					file << u <<'\t' << v  << '\t' << weight << '\n';
				}
				else if (gtype=='d') {
					file << u <<'\t' << v  << '\t' << weight << '\n';
				}
				else if (gtype=='g') {
					file << u <<'\t' << v  << '\t' << weight << '\n';
				}
				else {
					std::cerr << gtype << " type of graph not recognized in writing graphs\n";
					throw 1;
				}
			}
		}
	}
	file.close();
};


void graphData::writeSingle(const char* fn) {
	int u,v;
	float weight;
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
			weight = inIt->second;
			if (u<v) {
				if (gtype=='w') {
					file << uname <<'\t' << vname << '\t' << weight << '\n';
				}
				else if (gtype=='b') {
					file << uname <<'\t' << vname << '\t' << weight <<'\n';
				}
				else if (gtype=='p') {
					file << uname <<'\t' << vname  << '\t' << weight << '\n';
				}
				else if (gtype=='d') {
					file << uname <<'\t' << vname  << '\t' << weight << '\n';
				}
				else if (gtype=='g') {
					file << uname <<'\t' << vname  << '\t' << weight << '\n';
				}
				else {
					std::cerr << gtype << " type of graph not recognized in writing graphs\n";
					throw 1;
				}
			}
		}
	}
	file.close();
};



void graphData::writeBoth(const char* fn) {
	int u,v;
	float weight;
	std::ofstream file;
	file.open(fn,std::ios::out);
	std::map<int, destList*>::iterator outIt;
	destList::iterator inIt;
	for (outIt = edgeList.begin(); outIt != edgeList.end(); ++outIt) {
		u = outIt->first;
		for (inIt = edgeList[u]->begin(); inIt != edgeList[u]->end(); ++inIt) {
			v = inIt->first;
			weight = inIt->second;
			if (gtype=='w') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
			}
			else if (gtype=='b') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\n';
			}
			else if (gtype=='p') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
			}
			else if (gtype=='d') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
			}
			else if (gtype=='g') {
				file << int2Name[u] <<'\t' << int2Name[v] << '\t' << weight << '\n';
			}
			else {
				std::cerr << gtype << " type of graph not recognized in writing graphs\n";
				throw 1;
			}
		}
	}
	file.close();	
};


bool graphData::has_uv(int u, int v) {
	if (edgeList.find(u)!=edgeList.end()) {
		if (edgeList[u]->find(v) != edgeList[u]->end()) {
			return 1;
		}
	}
	return 0;
};


float graphData::get_uv(int u, int v) {
	std::map<int,destList*>::iterator outIt;
	outIt = edgeList.find(u);
	if (outIt != edgeList.end()) {
		destList::iterator inIt;
		inIt = (*outIt).second->find(v);
		if (inIt != (*outIt).second->end()) {
			return (*inIt).second;
		}
	}
	return 0;
};

bool graphData::readGeneralBasedOnOld(graphData* Goriginal, const char* filename, int nsoFar, const char* sep) {
	// NOTE: remember to set the graph type
	this->gtype = 'u'; // u for unknown, it is a placeholder
	this->Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v,ii;
	float sum = 0;
	float weight;
	std::string Vt1, Vt2;
	int Vg1, Vg2;
	this->int2Name = Goriginal[0].int2Name;
	this->name2Int = Goriginal[0].name2Int;
	this->typeList = Goriginal[0].typeList;
	this->numV = Goriginal->numV;
	this->multiPartite = 0;
	file.open(filename,std::ios::in);
	if (! file.is_open()) {throw(std::runtime_error("file " + std::string(filename)+"not open.\n")); return 0;}
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok,sep);
		if (tok.size()==2) {
			Vt1 = tok[0];
			Vt2 = tok[1];
			weight = 1; Vg1=0; Vg2=0;
		} else if (tok.size()==3) {
			Vt1 = tok[0];
			Vt2 = tok[1];
			weight = boost::lexical_cast<float>(tok[2]);
			Vg1 = 0; Vg2 = 0;
		} else if (tok.size()==4) {
			Vt1 = tok[0];
			Vt2 = tok[2];
			weight = 1;
			Vg1 = boost::lexical_cast<int>(tok[1]);
			Vg2 = boost::lexical_cast<int>(tok[3]);
			this->multiPartite = 1;
		} else if (tok.size()==5) {
			Vt1 = tok[0];
			Vt2 = tok[2];
			weight = boost::lexical_cast<float>(tok[4]);
			Vg1 = boost::lexical_cast<int>(tok[1]);
			Vg2 = boost::lexical_cast<int>(tok[3]);
			this->multiPartite = 1;
		} else {
			// nothing to do for this line
			if (DEBUGMODE) throw(std::runtime_error("Error: don't know what to do with a line of " + boost::lexical_cast<std::string>(tok.size()) + "tokens at line = " + boost::lexical_cast<std::string>(strline)+"\n"));
			continue;
		}
		if (Vt1==Vt2) {
			throw(std::runtime_error("self loop in " + std::string(filename) +", line = " + boost::lexical_cast<std::string>(strline)+"\n"));
		}
		
		if (name2Int.find(Vt1)==name2Int.end()) {
			u = name2Int.size();
			name2Int[Vt1] = u;
			int2Name[u] = Vt1;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].int2Name[u] = Vt1;
				Goriginal[ii].name2Int[Vt1] = u;
			}
		} else {
			u = name2Int[Vt1];
		}
		if (name2Int.find(Vt2)==name2Int.end()) {
			v = name2Int.size();
			name2Int[Vt2] = v;
			int2Name[v] = Vt2;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].int2Name[v] = Vt2;
				Goriginal[ii].name2Int[Vt2] = v;
			}
		} else {
			v = name2Int[Vt2];
		}
		if (edgeList.find(u)==edgeList.end()) {
			pdl = new destList;
			edgeList[u] = pdl;
		}
		if (edgeList[u]->find(v)==edgeList[u]->end()) {
			(*edgeList[u])[v] = weight;
		} else {
			(*edgeList[u])[v] += weight;
		}


		if (edgeList.find(v)==edgeList.end()) {
			pdl = new destList;
			edgeList[v] = pdl;
		}
		if (edgeList[v]->find(u)==edgeList[v]->end()) {
			(*edgeList[v])[u] = weight;
		} else {
			(*edgeList[v])[u] += weight;
		}
		if (typeList.find(u)==typeList.end()) {
			typeList[u] = Vg1;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].typeList[u] = Vg1;
			}
		} else if ((typeList[u]==0)and(Vg1 != 0)) {
			typeList[u] = Vg1;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].typeList[u] = Vg1;
			}
		}
		if (typeList.find(v)==typeList.end()) {
			typeList[v] = Vg2;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].typeList[v] = Vg2;
			}
		} else if ((typeList[v]==0)and(Vg2 != 0)) {
			typeList[v] = Vg2;
			for (ii=0;ii<nsoFar;ii++) {
				Goriginal[ii].typeList[v] = Vg2;
			}
		}

		sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
				// weighted graphs
		Etot += weight;
	}
	numV = name2Int.size();
	for (ii=0;ii<nsoFar;ii++) {
		Goriginal[ii].numV = numV;
	}
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
#if DEBUGMODE
	std::cout << "read graph, V = " << numV << ", E = " << sum/2<<"\n";
#endif
	return 1;
};





/*bool graphData::readGeneralBasedOnOld(graphData* Goriginal, const char* filename, int nsoFar) {
	//gtype = 'b';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u, v, ii;
	float sum = 0;
	float weight;
	int linesread=0;
	std::string St1, St2;
	//std::istringstream temp;
	file.open(filename,std::ios::in);
	if (! file.is_open()) return 0;
	this->int2Name = Goriginal[0].int2Name;
	this->name2Int = Goriginal[0].name2Int;
	this->numV = Goriginal->numV;
	while (! file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			linesread++;
			if (name2Int.find(tok[0])==name2Int.end()) {
				u = name2Int.size();
				name2Int[tok[0]] = u;
				int2Name[u] = tok[0];
				for (ii=0;ii<nsoFar;ii++) {
					Goriginal[ii].int2Name[u] = tok[0];
					Goriginal[ii].name2Int[tok[0]] = u;
				}
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				v = name2Int.size();
				name2Int[tok[1]] = v;
				int2Name[v] = tok[1];
				for (ii=0;ii<nsoFar;ii++) {
					Goriginal[ii].int2Name[v] = tok[1];
					Goriginal[ii].name2Int[tok[1]] = v;
				}
			} else {
				v = name2Int[tok[1]];
			}
			if (tok.size()>2) {
				std::istringstream(tok[2]) >> weight;
			} else weight = 1;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = weight;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = weight;
			sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
					// weighted graphs
			Etot += weight;
		}
	}
	numV = name2Int.size();
	for (ii=0;ii<nsoFar;ii++) {
		Goriginal[ii].numV = Goriginal[ii].name2Int.size();
	}
	
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
	return 1;
};*/



bool graphData::readBinaryBasedOnOld(graphData* Goriginal, const char* filename) {
	gtype = 'b';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	int linesread=0;
	std::string St1, St2;
	//std::istringstream temp;
	file.open(filename,std::ios::in);
	if (! file.is_open()) return 0;
	this->int2Name = Goriginal->int2Name;
	this->name2Int = Goriginal->name2Int;
	this->numV = Goriginal->numV;
	while (! file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			linesread++;
			if (name2Int.find(tok[0])==name2Int.end()) {
				// do nothing
				continue;
				return 0;
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				// do nothing
				continue;
				return 0;
			} else {
				v = name2Int[tok[1]];
			}
			weight = 1;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = weight;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = weight;
			sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
					// weighted graphs
			Etot += weight;
		}
	}
	numV = name2Int.size();
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
	return 1;
};


bool graphData::readGeneral(const char* fn, const char* sep) {
	// NOTE: remember to set the graph type
	this->gtype = 'u'; // u for unknown, it is a placeholder
	this->Etot = 0.0f;
	this->multiPartite = 0;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	std::string Vt1, Vt2;
	int Vg1, Vg2;
	//std::istringstream temp;
	file.open(fn,std::ios::in);
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok,sep);
		if (tok.size()==2) {
			Vt1 = tok[0];
			Vt2 = tok[1];
			weight = 1; Vg1=0; Vg2=0;
		} else if (tok.size()==3) {
			Vt1 = tok[0];
			Vt2 = tok[1];
			weight = boost::lexical_cast<float>(tok[2]);
			Vg1 = 0; Vg2 = 0;
		} else if (tok.size()==4) {
			Vt1 = tok[0];
			Vt2 = tok[2];
			weight = 1;
			Vg1 = boost::lexical_cast<int>(tok[1]);
			Vg2 = boost::lexical_cast<int>(tok[3]);
			this->multiPartite = 1;
		} else if (tok.size()==5) {
			Vt1 = tok[0];
			Vt2 = tok[2];
			weight = boost::lexical_cast<float>(tok[4]);
			Vg1 = boost::lexical_cast<int>(tok[1]);
			Vg2 = boost::lexical_cast<int>(tok[3]);
			this->multiPartite = 1;
		} else {
			// nothing to do for this line
			if (DEBUGMODE) std::cout << "Error: don't know what to do with a line of "<<tok.size()<<"tokens!\n";
			continue;
		}

		if (Vt1==Vt2) {
			std::cout << "self loop in "<<fn<<", line = "<<strline<<"\n";
		}
		
		if (name2Int.find(Vt1)==name2Int.end()) {
			u = name2Int.size();
			name2Int[Vt1] = u;
			int2Name[u] = Vt1;
		} else {
			u = name2Int[Vt1];
		}
		if (name2Int.find(Vt2)==name2Int.end()) {
			v = name2Int.size();
			name2Int[Vt2] = v;
			int2Name[v] = Vt2;
		} else {
			v = name2Int[Vt2];
		}
		if (edgeList.find(u)==edgeList.end()) {
			pdl = new destList;
			edgeList[u] = pdl;
		}
		if (edgeList[u]->find(v)==edgeList[u]->end()) {
			(*edgeList[u])[v] = weight;
		} else {
			(*edgeList[u])[v] += weight;
		}


		if (edgeList.find(v)==edgeList.end()) {
			pdl = new destList;
			edgeList[v] = pdl;
		}
		if (edgeList[v]->find(u)==edgeList[v]->end()) {
			(*edgeList[v])[u] = weight;
		} else {
			(*edgeList[v])[u] += weight;
		}
		if (typeList.find(u)==typeList.end()) {
			typeList[u] = Vg1;
		} else if ((typeList[u]==0)and(Vg1 != 0)) {
			typeList[u] = Vg1;
		}
		if (typeList.find(v)==typeList.end()) {
			typeList[v] = Vg2;
		} else if ((typeList[v]==0)and(Vg2 != 0)) {
			typeList[v] = Vg2;
		}

		sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
				// weighted graphs
		Etot += 2*weight;
	}
	numV = name2Int.size();
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
#if DEBUGMODE
	std::cout << "read graph, V = " << numV << ", E = " << sum/2<<"\n";
#endif
	return 1;
};




bool graphData::readWeighted(const char* fn) {
	gtype = 'w';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	int u,v;
	float sum = 0;
	float weight;
	std::string St1, St2;
	//std::istringstream temp;
	file.open(fn,std::ios::in);
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			if (name2Int.find(tok[0])==name2Int.end()) {
				u = name2Int.size();
				name2Int[tok[0]] = u;
				int2Name[u] = tok[0];
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				v = name2Int.size();
				name2Int[tok[1]] = v;
				int2Name[v] = tok[1];
			} else {
				v = name2Int[tok[1]];
			}
			if (tok.size()>2) {
				std::istringstream(tok[2]) >> weight;
			} else weight = 1;
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = weight;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = weight;
			sum += 2; // NOTE : this shouldn't really matter as the sum isn't used for
					// weighted graphs
			Etot += weight;
		}
	}
	numV = name2Int.size();
	aveP = sum/(numV * numV); // shouldn't matter because aveP is to be used for binary
	return 1;
};

void graphData::set_uv(int u, int v, float w) {
	std::map<int, destList*>::iterator outIt;
	outIt = edgeList.find(u);
	if (outIt==edgeList.end()) {
		destList* pdest = new destList;
		(*pdest).insert((*pdest).begin(), std::pair<int,float>(v,w));
		edgeList.insert(outIt,std::pair<int, destList*>(u,pdest));
		//edgeList[u] = pdest;
	} else {
		(*(*outIt).second)[v] = w;
	}
};

bool graphData::delete_uv(int u, int v) {
	if (edgeList.find(u)==edgeList.end()) {
		return 0;
	} else {
		edgeList[u]->erase(v);
		return 1;
	}
}

int graphData::Add_uv(int u, int v, float w) {
	if (edgeList.find(u)==edgeList.end()) {
		edgeList[u] = new destList;
		(*edgeList[u])[v] = w;
		return 0;
	} else {
		if ((*edgeList[u]).find(v) != (*edgeList[u]).end()) {
			(*edgeList[u])[v] += w;
			return 1;
		} else {
			(*edgeList[u])[v] = w;
			return 0;
		}
	}
};

bool graphData::readBinary(const char* fn) {
	gtype = 'b';
	Etot = 0.0f;
	std::string strline;
	std::ifstream file;
	std::vector<std::string> tok;
	destList* pdl;
	//std::istringstream temp;
	int u,v;
	float sum = 0;
	file.open(fn,std::ios::in);
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,strline);
		tok.clear();
		my_Tokenize(strline,tok," \t");
		if (tok.size()>1) {
			if (name2Int.find(tok[0])==name2Int.end()) {
				u = name2Int.size();
				name2Int[tok[0]] = u;
				int2Name[u] = tok[0];
			} else {
				u = name2Int[tok[0]];
			}
			if (name2Int.find(tok[1])==name2Int.end()) {
				v = name2Int.size();
				name2Int[tok[1]] = v;
				int2Name[v] = tok[1];
			} else {
				v = name2Int[tok[1]];
			}
			if (edgeList.find(u)==edgeList.end()) {
				pdl = new destList;
				edgeList[u] = pdl;
			}
			(*edgeList[u])[v] = 1.0f;
			if (edgeList.find(v)==edgeList.end()) {
				pdl = new destList;
				edgeList[v] = pdl;
			}
			(*edgeList[v])[u] = 1.0f;
			sum += 2;
			Etot += 1.0f;
		}
	}
	numV = name2Int.size();
	aveP = sum/(numV * numV);
	return 1;
};

void my_Tokenize(const std::string& str,std::vector<std::string>& tokens,const std::string& delimiters = "")  // via http://www.geocities.com/eric6930/cplus.html unknown author presumed Eric.
{
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastPos);
    }
}





#endif

