#define DEBUGMODE 0
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
#include <iostream>
#include <string> 



int main(int argc, char* argv[]) {
	std::string fnstem,fnin,fnout;
	char temp[80];
	linkPredictorOther lp;
	int numResids;
	if (argc==1) {
		std::cout << "Need filename for input graph\n";
		throw 1;
	}
	if (argc==2) {
		numResids = 6;
		std::cout << "Assuming default: "<< numResids << " residuals to compute.\n";
	} else {
		std::stringstream ssss;
		ssss << argv[2];
		ssss >> numResids;
	}
	fnstem = argv[1];
	//fnstem = "grassland";
	fnin = fnstem+".edges";
	// let us read in the graph
	graphData *Goriginal, *Gnew, *GPred, *GsoFar, *Glabel;
	Engine *en;
	Goriginal = new graphData[1];
	GsoFar = new graphData[1];
	Glabel = new graphData[1];
	Goriginal[0].readBinary(fnin.c_str());
	std::cout << "read file!\n";
	fnin = fnstem + ".label";
	Glabel[0].readBinaryBasedOnOld(& Goriginal[0],fnin.c_str());
	Goriginal[0].copyNoEdges(GsoFar[0]);
	en = new Engine(Goriginal,Goriginal, GsoFar, 1);
	en->initializeScores();
	/*
	// let us write the competing scores
	std::cout << "writing degree product\n";
	fnout = fnstem + ".dprod"; en->printDegreeProdFile(fnout.c_str(), 0, 0);
	//std::cout << "writing hyper geom\n";
	//fnout = fnstem + ".hyperg"; en->printHyperGeomFile(fnout.c_str(), 0, 0);
	std::cout << "writing jaccard\n";
	fnout = fnstem + ".jaccard"; en->printJaccardFile(fnout.c_str(), 0, 0);
	std::cout << "writing comm neigh\n";
	fnout = fnstem + ".cneighb"; en->printCommonNeighbFile(fnout.c_str(), 0, 0);
	*/

	// run the agglomerative algorithm
	std::cout << "starting to run\n";
	en->run();
	std::cout << "done running\n";
	lp.attach(en); // attached the engine
	std::cout << "Attached!\n";
	GPred = lp.makeEdgePred(Glabel);
	fnout = fnstem + ".scores0";
	GPred[0].writeSingle(fnout.c_str());
	delete[] GPred;

	fnout = fnstem + ".nodes";
	//en->tree->writeNodeTypes(fnout.c_str());
	



	fnout = fnstem + "0.clusters";
	en->tree->writeCollapsedHierEdges(fnout.c_str());

	int residint = 1;
	while(residint <= numResids) {
		std::stringstream dummy;
		dummy << residint;
		std::string sres;
		dummy >> sres;
		
		lp.updateSoFar(GsoFar);
		fnout = fnstem + sres + ".soFar";
		GsoFar->writeSingle(fnout.c_str());
		std::cout << "getting the residual graph\n";
		//Gnew = getResidual(G, lp);
		Gnew = residualDiff(Goriginal, GsoFar, 1);
		std::cout << "got the residual\n";
		std::cout << "Etot = "<< Gnew->Etot << "\n";
		delete en;
		en = new Engine(Gnew,Goriginal,GsoFar,1);
		en->initializeScores();
		std::cout << "running on the residual\n";
		en->run();
		lp.attach(en);

		GPred = lp.makeEdgePred(Goriginal);
		fnout = fnstem + ".scores" + sres;
		GPred[0].writeSingle(fnout.c_str());


		//fnout = fnout + "_noname";
		//GPred[0].writeSingle_noname(fnout.c_str());
		std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";

		fnout = fnstem + sres + ".clusters";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		residint++;
		delete[] GPred;
		delete[] Gnew;
	}
	std::cout << "done!\n";
	delete[] Goriginal;
	delete[] GsoFar;
	delete[] Glabel;
	delete en;
	//delete G;
	return 1;
}
