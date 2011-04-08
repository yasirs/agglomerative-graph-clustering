#define DEBUGMODE 1
#define NOGSL 0
#define NOREFERENCE 0
#include "engineData.hpp"
#include "graphData.hpp"
#include <iostream>
#include <string>



int main(int argc, char* argv[]) {
	char dummy[100];
	if (argc==1) {
		std::cout << "Need filename for input graph\n";
		throw 1;
	}
	std::cin >> dummy;
	std::string fnstem,fnin,fnout;
	fnstem = argv[1];
	fnin = fnstem+".txt";
	// let us read in the graph
	graphData *G;
	G = new graphData[1];
	G[0].readBinary(fnin.c_str());
	std::cout << "read file!\n";
	Engine en(G,1);
	en.initializeFirstLev();
	en.run();
	std::cout << "done running\n";
	fnout = fnstem+".network";
	en.tree->writeCollapsedHierEdges(fnout.c_str());
	fnout = fnstem + ".nodes";
	en.tree->writeNodeTypes(fnout.c_str());
	delete[] G;
	//delete G;
	return 1;
}
