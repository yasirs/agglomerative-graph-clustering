//#include "mysetfuncs.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <algorithm>


std::vector<std::string>* splitspaces(const std::string& str, const std::string& delimiters = " ")
{
    std::vector<std::string>* tokens = new std::vector<std::string>;
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens->push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
};



class myEdge{
	public:
		std::string v1, v2;
		myEdge(std::string line) {
			std::string myline = line;
			std::string endlines ("\n\r");
			size_t found;
			found = myline.find_last_not_of(endlines);
			if (found!=std::string::npos)
    				myline.erase(found+1);
			else
				myline.clear();
                        std::vector<std::string>* v;
                        v = splitspaces(myline," \t");
                        if ((*v)[0]< (*v)[1]) {
                                this->v1 = (*v)[0];
                                this->v2 = (*v)[1];
                        } else {
                                v2 = (*v)[0];
                                v1 = (*v)[1];
                        }
			assert(v1.compare(std::string())!=0);
			assert(v2.compare(std::string())!=0);
                }
		/*myEdge(const char* s) {
			std::string str(s);
			this->myEdge(s);
		}*/
                bool operator==(const myEdge& m) const {
                        if (v1.compare(m.v1)==0)
				return (v2.compare(m.v2)==0);
			else {
				if (v1.compare(m.v2)==0)
					return (v2.compare(m.v1)==0);
				else
					return 0;
			}
                }
		myEdge(int i, int j) {
			char dummy[20];
			sprintf(dummy,"%d",i);	this->v1 = std::string(dummy);
			sprintf(dummy,"%d",j);	this->v2 = std::string(dummy);
		}
		myEdge(const std::string& vert1, const std::string& vert2) {
			this->v1 = vert1; this->v2=vert2;
			assert(v1.compare(std::string())!=0);
			assert(v2.compare(std::string())!=0);
		}
};


struct EdgeHash{
	public:
		size_t operator() (const myEdge& E) const {
			return (std::tr1::hash<std::string>() (E.v1) + std::tr1::hash<std::string>() (E.v2));
			//return std::tr1::hash<std::string>() (E.v1);
		}
};





int main(int argc, char* argv[]) {
	std::tr1::unordered_set<myEdge, EdgeHash> edset;
	std::tr1::unordered_map<int, std::string> int2Name;
	std::tr1::unordered_map<std::string, int> name2Int;
	float holdout;
	unsigned long long int tstamp;
	std::string fnstem, fnin, fnout, line;
	unsigned int nlines, i, nhold, j, k, nverts;
	if (argc<2) {
		std::cerr << "Need name for input graph (.txt will be appended)\n";
		throw(1);
	}
	if (argc<3) {
		std::cerr << "need fraction of holdout\n";
		throw(1);
	}
	std::istringstream(argv[2]) >> holdout;
	if (holdout>1) {
		std::cerr << "fraction needs to be <1!\n";
		throw(1);
	}
	tstamp = clock(); //time(NULL);
	//std::cin >> tstamp;
	srand( tstamp );
	std::set<unsigned int> eset;
	fnstem = argv[1];
	fnin = fnstem + ".txt";
	std::ifstream file;
	std::ofstream filegold, fileedges;
	file.open(fnin.c_str(),std::ios::in);
	nlines = 0;
	if (! file.is_open()) return 0;
	while (!file.eof()) {
		getline(file,line);
		nlines++;
	}
	file.close();
	nhold = holdout * nlines;
	nhold = std::max(nhold,(unsigned ) 1);
	nhold = std::min(nhold,nlines-1);
	//std::cout << "nlines = "<< nlines << "\n";
	//std::cout << "nhold = "<< nhold << "\n";
	//std::cout << "holdout = "<< holdout << "\n";
	std::cout << tstamp;
	std::cout.flush();
	while (eset.size()<nhold) {
		i = rand() % nlines;
		eset.insert(i);
	}
	std::vector<unsigned int> evec (eset.begin(), eset.end());
	std::sort (evec.begin(), evec.end());
	file.open(fnin.c_str(),std::ios::in);
	fnout = fnstem + ".labels";
	filegold.open(fnout.c_str(), std::ios::out);
	fnout = fnstem + ".edges";
	fileedges.open(fnout.c_str(), std::ios::out);
	if (! file.is_open()) return 0;
	i = 0;
	std::vector<unsigned int>::iterator nit = evec.begin();
	while (!file.eof()) {
		getline(file,line);
		if (line.compare(std::string())!=0) {
			myEdge EE(line);
			if ((*nit)==i) {
				++nit;
				// this is a hold out edge
				filegold << line << "\t1\n"; 
			} else {
				fileedges << line << "\n";
			}
			edset.insert(EE);
			if (name2Int.find(EE.v1)==name2Int.end()) {
				int newint = int2Name.size();
				name2Int[EE.v1] = newint;
				int2Name[newint] = EE.v1;
			}
			if (name2Int.find(EE.v2)==name2Int.end()) {
				int newint = int2Name.size();
				name2Int[EE.v2] = newint;
				int2Name[newint] = EE.v2;
			}
			i++;
		}
	}
	file.close();
	fileedges.close();
	nverts = name2Int.size();
	std::tr1::unordered_set<myEdge, EdgeHash> negset;
	while (negset.size()<nhold) {
		j = rand() % nverts;
		k = rand() % nverts;
		if (j != k) {
			myEdge EE(int2Name[j], int2Name[k]);
			if (edset.find(EE)==edset.end()) {
				negset.insert(EE);
			} 
		}
	}
	std::tr1::unordered_set<myEdge, EdgeHash>::iterator edgeit(negset.begin());
	for (;edgeit != negset.end(); edgeit++) {
		filegold << (*edgeit).v1 << "\t" << (*edgeit).v2 << "\t-1\n";
		//std::cout << (*edgeit).v1 << "\t" << (*edgeit).v2 << "\t-1\n";
	}
	filegold.close();


	return 0;
};

