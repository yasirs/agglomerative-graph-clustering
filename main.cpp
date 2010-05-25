#define DEBUGMODE 1
#define NOGSL 1
#define NOREFERENCE 0
#include "Engine.hpp"
#include "graphData.hpp"
#include "residuals.hpp"
#include <iostream>
#include <string> 



int main(int argc, char* argv[]) {
	std::string fnstem,fnin,fnout;
	char dummy[100];
	if (argc==1) {
		std::cout << "Need filename for input graph\n";
		throw 1;
	}
	std::cin >> dummy;
	fnstem = argv[1];
	//fnstem = "grassland";
	fnin = fnstem+".txt";
	// let us read in the graph
	graphData *G,*Gnew;
	Engine *en;
	G = new graphData[1];
	G[0].readBinary(fnin.c_str());
	std::cout << "read file!\n";
	en = new Engine(G,1);
	en->initializeFirstLev();
	en->run();
	std::cout << "done running\n";
	fnout = fnstem + ".network";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	fnout = fnstem + ".nodes";
	en->tree->writeNodeTypes(fnout.c_str());
	std::cout << "getting the residual graph\n";
	Gnew = getResidual(G, en->tree, en->w, 1);
	std::cout << "Etot = "<< Gnew->Etot << "\n";
	int residint = 1;
	while(Gnew->Etot > .1f) {
		std::cout << "got the residual\n";
		delete en;
		en = new Engine(Gnew,1);
		std::cout << "running on the residual\n";
		en->initializeFirstLev();
		en->run();
		residint++;
		std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";
		fnout = fnstem + "2.network";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		delete[] G;
		G = Gnew;
		std::cout << "getting the residual graph\n";
		Gnew = getResidual(G, en->tree, en->w, 1);
		std::cout << "Etot = "<< Gnew->Etot << "\n";
	}
	std::cout << "done!\n";
	delete[] G;
	delete en;
	//delete G;
	return 1;
}
