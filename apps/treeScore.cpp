#define DEBUGMODE 1
#define NOGSL 1
#define NOREFERENCE 0
#define ISVC 0
#if ISVC
#include "C:\Users\Suhail\Documents\Visual Studio 2010\Projects\g\g\vchead.hpp"
#endif
#include "Engine.hpp"
#include "graphData.hpp"
#include "residuals.hpp"
#include "linkPredictor.hpp"
#include "linkPredictorOther.hpp"
#include "labelData.hpp"
#include <iostream>
#include <string> 


bool fexists(const char* filename) {
	std::ifstream ifile(filename);
	return ifile;
}


int main(int argc, char* argv[]) {
	std::string fnstem,fnin,fnout, fjoin;
	linkPredictorOther lp;
	int numResids;
	char gtype;
	bool doScores;
	labelData* Glabel;
	if (argc==1) {
		std::cout << "Need filename for input graph\n";
		throw 1;
	}
	if (argc==2) {
		std::cout << "Need filename containing tree joins\n";
		throw 1;
	}
	if (argc>3) {
		gtype = argv[3][0];
	} else {
		gtype = 'b';
		std::cout << "graph type = " << gtype << '\n';
	}
	fnstem = argv[1];
	fjoin = argv[2];
	fnin = fnstem+".edges";
	// let us read in the graph
	graphData *Goriginal, *Gnew, *GsoFar;
	Engine *en;
	Goriginal = new graphData[1];
	GsoFar = new graphData[1];
	Goriginal[0].readGeneral(fnin.c_str());
	Goriginal[0].gtype = gtype;
	std::cout << "read file!\n";
	fnin = fnstem + ".labels";
	doScores = fexists(fnin.c_str());
	if (doScores) {
		Glabel = new labelData(& Goriginal[0], fnin.c_str(), numResids+1+3);
	}
	Goriginal[0].copyNoEdges(GsoFar[0]);
	std::cout << "Initializing Engine, scores\n";
	en = new Engine(Goriginal,Goriginal, GsoFar, 1);
	en->initializeScoresML();

	// read the tree join
	std::cout << "starting to run\n";
	en->readJoins(fjoin.c_str());
	std::cout << "done running\n";
	lp.attach(en); // attached the engine
	std::cout << "Attached!\n";
	
	fnout = fnstem + ".clusters";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	lp.updateSoFarLazy(GsoFar);
	fnout = fnstem + ".soFar";
	GsoFar->writeSingle(fnout.c_str());
	if (doScores)
		Glabel->putSoFar(GsoFar, 1);

	fnin = fnstem + ".edges";
	if (doScores) {
		Glabel->populateLocal(fnin.c_str(), 1);
		fnout = fnstem+".labelScores";
		Glabel->write(fnout.c_str());
	}
	std::cout << "done!\n";
	delete[] Goriginal;
	delete[] GsoFar;
	delete en;
	//delete G;
	return 1;
}
