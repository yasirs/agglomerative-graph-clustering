#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <cstdlib>
#include <vector>
//#include "mysetfuncs.hpp"



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
		myEdge() {};
		myEdge(std::string& line) {
			std::vector<std::string>* v;
			v = splitspaces(line," \t");
			if ((*v)[0]< (*v)[1]) {
				v1 = (*v)[0];
				v2 = (*v)[1];
			} else {
				v2 = (*v)[0];
				v1 = (*v)[1];
			}
		}
		bool operator==(const myEdge& m) const {
			return ((v1.compare(m.v1)==0)and(v2.compare(m.v2)==0));
		}
		bool operator<(const myEdge& m) const {
			if (v1.compare(m.v1)==0) return (v2.compare(m.v2)<0);
			else return (v1.compare(m.v1)<0);
		}
			
};
	

int main(int argc, char* argv[]) {
	myEdge e;
	std::vector<std::string> *words;
	std::string fnstem, fnin, fnout, line;
	if (argc<3) {
		std::cerr << "Usage:\nwriteLabels scorefile goldfile outfile\n";
		throw(1);
	}
	std::ifstream scfile, gfile;
	std::ofstream ofile;
	std::set<myEdge> eset;
	gfile.open(argv[2],std::ios::in);
	if (! gfile.is_open()) return 0;
	while (!gfile.eof()) {
		getline(gfile,line);
		if (line.compare("")!=0) {
			eset.insert(myEdge(line));
		}
	}
	gfile.close();
	scfile.open(argv[1],std::ios::in);
	ofile.open(argv[3],std::ios::out);
	if (! scfile.is_open()) return 0;
	if (! ofile.is_open()) return 0;
	while (!scfile.eof()) {
		getline(scfile,line);
		if (line.compare("")!=0) {
			e = myEdge(line);
			words = splitspaces(line," \t");
			if (eset.find(e)==eset.end()) {
				//label 0
				ofile << (*words)[0] << '\t' << (*words)[1] << "\t0\n";
			} else {
				ofile << (*words)[0] << '\t' << (*words)[1] << "\t1\n";
			}
		}
	}
	scfile.close();
	ofile.close();
	return 0;
};

