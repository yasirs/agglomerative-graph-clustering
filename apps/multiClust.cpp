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
	char gtype;
	bool doScores;
	labelData* Glabel;
	int nGraphs;
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
	if (argc>3) {
		gtype = argv[3][0];
	} else {
		gtype = 'b';
		std::cout << "graph type = " << gtype << '\n';
	}
	if (argc>4) {
		std::stringstream ssss;
		ssss << argv[4];
		ssss >> nGraphs;
		std::cout << "Assuming " << nGraphs << " input graphs.\n";
	} else {
		nGraphs = 1;
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
			Goriginal[ii].readWeighted(fnin.c_str());
		else
			Goriginal[ii].readGeneralBasedOnOld(Goriginal, fnin.c_str(), ii);
		Goriginal[ii].gtype = gtype;
		Goriginal[ii].copyNoEdges(GsoFar[ii]);
	}
	std::cout << "read file!\n";
	fnin = fnstem + ".labels";
	doScores = fexists(fnin.c_str());
	if (doScores) {
		Glabel = new labelData(& Goriginal[0], fnin.c_str(), numResids+1+3);
	}
	en = new Engine(Goriginal,Goriginal, GsoFar, nGraphs);
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
		if (residint>1) delete[] Gnew;
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
	delete[] Gnew;
	delete en;
	//delete G;
	return 1;
}
