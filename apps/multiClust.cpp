#define DEBUGMODE 1
#define NOGSL 0
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
	std::string fnstem,fnin,fnout;
	linkPredictorOther lp;
	int numResids;
	char* gtype;
	bool doScores;
	labelData* Glabel;
	int nGraphs;
	if (argc<5) {
		std::cout << "Usage :"<<argv[0]<<" FileStem nGraphs Type_1 Type_2 ... Type_nGraphs nResiduals\n";
		std::cout << "Minimum nGraphs = 1\nMinimum nResiduals = 0\n";
		throw(1);
	}
	{ // block just to have temporary variable ssss
		std::stringstream ssss;
		ssss << argv[2];
		ssss >> nGraphs;
	}
	std::cout << "Expecting " << nGraphs << " input graphs.\n";
	if (nGraphs<1) {
		std::cout << "Need at least 1 input graph.\n";
		throw(1);
	}
	gtype = new char[nGraphs];
	if (argc < (3 + nGraphs)) {
		std::cout << "Usage :"<<argv[0]<<" FileStem nGraphs Type_1 Type_2 ... Type_nGraphs nResiduals\n";
		std::cout << "Minimum nGraphs = 1\nMinimum nResiduals = 0\n";
		throw(1);
	}
	for (int i=0;i<nGraphs;i++) {
		gtype[i] = argv[i+3][0];
	}
	{
		std::stringstream ssss;
		ssss << argv[nGraphs+3];
		ssss >> numResids;
	}

	fnstem = argv[1];
	graphData *Goriginal, *Gnew, *GsoFar;
	Engine *en;
	Goriginal = new graphData[nGraphs];
	GsoFar = new graphData[nGraphs];
	for (int ii=0; ii<nGraphs; ii++) {
		std::stringstream sss;
		sss << ii;
		std::string iis;
		sss >> iis;
		fnin = fnstem + iis + ".edges";
		if (ii==0)
			Goriginal[ii].readGeneral(fnin.c_str());
		else
			Goriginal[ii].readGeneralBasedOnOld(
				&(Goriginal[0]), 
				fnin.c_str(), 
				ii);
		Goriginal[ii].gtype = gtype[ii];
		Goriginal[ii].copyNoEdges(GsoFar[ii]);
	}
	std::cout << "read file!\n";
	fnin = fnstem + ".labels";
	doScores = fexists(fnin.c_str());
	if (doScores) {
		Glabel = new labelData(& Goriginal[0], fnin.c_str(), numResids+1+3);
	}
	en = new Engine(Goriginal, Goriginal, GsoFar, nGraphs);
	en->initializeScoresML();

	// run the agglomerative algorithm
	std::cout << "starting to run\n";
	en->runML();
	en->passFB();
	std::cout << "done running\n";
	lp.attach(en); // attached the engine
	std::cout << "Attached!\n";
	fnout = fnstem + ".nodes";
	en->tree->writeNodeTypes(fnout.c_str());
	
	fnout = fnstem + "0.clusters";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	lp.updateSoFarLazy(GsoFar);
	fnout = fnstem + "0.soFar";
	GsoFar[nGraphs-1].writeSingle(fnout.c_str());
	if (doScores)
		Glabel->putSoFar(GsoFar, nGraphs);

	int residint = 1;
	while(residint <= numResids) {
		std::stringstream dummy;
		dummy << residint;
		std::string sres;
		dummy >> sres;
		std::cout << "getting the residual graph\n";
		//Gnew = getResidual(G, lp);
		Gnew = residualDiff(Goriginal, GsoFar, nGraphs);
		fnout = fnstem + sres + ".residual"; Gnew[nGraphs-1].writeSingle(fnout.c_str());
		std::cout << "got the residual\n";
		for (int ii=0; ii<nGraphs; ii++) {
			std::cout << "Etot graph " << ii << "= " << Gnew[ii].Etot << "\n";
		}
		delete en;
		en = new Engine(Gnew,Goriginal,GsoFar, nGraphs);
		en->initializeScoresML();
		std::cout << "running on the residual\n";
		en->runML();
		en->passFB();
		lp.attach(en);
		//fnout = fnstem + ".scores" + sres;
		std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";
		fnout = fnstem + sres + ".clusters";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		lp.updateSoFarLazy(GsoFar);
		fnout = fnstem + sres + ".soFar";
		GsoFar[nGraphs-1].writeSingle(fnout.c_str());
		if (doScores) 
			Glabel->putSoFar(GsoFar, residint+1);
		if (residint>1) delete[] Gnew;
		residint++;
	}
	fnin = fnstem + ".edges";
	if (doScores) {
		Glabel->populateLocal(fnin.c_str(), residint+1);
		fnout = fnstem+".labelScores";
		Glabel->write(fnout.c_str());
	}
	std::cout << "done!\n";
	delete[] Goriginal;
	delete[] GsoFar;
	//delete[] Gnew;
	delete en;
	//delete G;
	return 1;
}
