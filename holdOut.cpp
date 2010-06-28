#include "mysetfuncs.hpp"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <cstdlib>

int main(int argc, char* argv[]) {
	float holdout;
	unsigned long long int tstamp;
	std::string fnstem, fnin, fnout, line;
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
	tstamp = time(NULL);
	srand( tstamp );
	std::set<int> eset;
	int nlines, i, n, nhold;
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
	nhold = std::max(nhold,1);
	nhold = std::min(nhold,nlines-1);
	//std::cout << "nlines = "<< nlines << "\n";
	//std::cout << "nhold = "<< nhold << "\n";
	//std::cout << "holdout = "<< holdout << "\n";
	std::cout << tstamp;
	std::cout.flush();
	n=0;
	while (eset.size()<nhold) {
		i = rand() % nlines;
		eset.insert(i);
	}
	std::vector<int> evec (eset.begin(), eset.end());
	std::sort (evec.begin(), evec.end());
	file.open(fnin.c_str(),std::ios::in);
	fnout = fnstem + ".gold";
	filegold.open(fnout.c_str(), std::ios::out);
	fnout = fnstem + ".edges";
	fileedges.open(fnout.c_str(), std::ios::out);
	if (! file.is_open()) return 0;
	i = 0;
	std::vector<int>::iterator nit = evec.begin();
	while (!file.eof()) {
		getline(file,line);
		if ((*nit)==i) {
			++nit;
			// this is a hold out edge
			filegold << line << "\n"; 
		} else {
			fileedges << line << "\n";
		}
		i++;
	}
	file.close();
	filegold.close();
	fileedges.close();
	return 0;
};

