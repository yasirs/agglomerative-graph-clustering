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


#if !defined(str2int_INCLUDED)
#define str2int_INCLUDED
int str2int(const std::string& input) {
	unsigned int output, p10;
	int i;
	p10 = 1;
	output = 0;
	for (i = input.length()-1; i>=0; i--) {
		output += p10*(int(input[i])-48);
		p10 = p10 * 10;
	}
	return output;
}

std::string int2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	std::string str  = "";
	if (input == 0) { str = "0"; } else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}


#endif



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
	fnstem = argv[1];
	//fnstem = "grassland";
	fnin = fnstem+".edges";
	// let us read in the graph
	graphData *Goriginal, *Gnew, *GsoFar;
	Engine *en;
	Goriginal = new graphData[1];
	GsoFar = new graphData[1];
	Goriginal[0].readGeneral(fnin.c_str());
	Goriginal[0].gtype = gtype;
	std::cout << "read file!\n";
	std::cout << "Initial Etot = "<<Goriginal->Etot << "\n";
	fnin = fnstem + ".labels";
	doScores = fexists(fnin.c_str());
	if (doScores) {
		Glabel = new labelData(& Goriginal[0], fnin.c_str(), numResids+1+3);
	}
	Goriginal[0].copyNoEdges(GsoFar[0]);
	std::cout << "Initializing Engine, scores\n";
	en = new Engine(Goriginal,Goriginal, GsoFar, 1);
	en->initializeScoresML();
	std::cout << "writing hyper geom\n";
	fnout = fnstem + ".hyperg"; en->printHyperGeomFile(fnout.c_str(), 0, 0);
	if (doScores) {
		// let us write the competing scores
		//std::cout << "writing degree product\n";
		//fnout = fnstem + ".dprod"; en->printDegreeProdFile(fnout.c_str(), 0, 0);
		//std::cout << "writing hyper geom\n";
		//fnout = fnstem + ".hyperg"; en->printHyperGeomFile(fnout.c_str(), 0, 0);
		//std::cout << "writing jaccard\n";
		//fnout = fnstem + ".jaccard"; en->printJaccardFile(fnout.c_str(), 0, 0);
		//std::cout << "writing comm neigh\n";
		//fnout = fnstem + ".cneighb"; en->printCommonNeighbFile(fnout.c_str(), 0, 0);
	}

	// run the agglomerative algorithm
	std::cout << "starting to run\n";
	en->runML(true,false);
	fnout = fnstem +"_0.hrg"; en->printHRG(fnout.c_str(),0);
	en->passFB();
	std::cout << "done running\n";
	lp.attach(en); // attached the engine
	std::cout << "Attached!\n";
	//fnout = fnstem + ".nodes";
	//en->tree->writeNodeTypes(fnout.c_str());
	
	fnout = fnstem + "0.clusters";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	lp.updateSoFarLazy(GsoFar);
	fnout = fnstem + ".0soFar";
	GsoFar->writeSingle(fnout.c_str());
	if (doScores)
		Glabel->putSoFar(GsoFar, 1);

	int residint = 1;
	while(residint <= numResids) {
		std::stringstream dummy;
		dummy << residint;
		std::string sres;
		dummy >> sres;
		std::cout << "getting the residual graph\n";
		Gnew = residualDiff(Goriginal, GsoFar, 1);
		fnout = fnstem + sres + ".residual"; Gnew->writeSingle(fnout.c_str());
		std::cout << "got the residual\n";
		std::cout << "Etot = "<< Gnew->Etot << "\n";
		delete en;
		en = new Engine(Gnew,Goriginal,GsoFar,1);
		en->initializeScoresML();
		std::cout << "running on the residual\n";
		en->runML(true,false);
		fnout = fnstem +"_"+int2str(residint)+".hrg"; en->printHRG(fnout.c_str(),0);
		en->passFB();
		lp.attach(en);
		//fnout = fnstem + ".scores" + sres;
		//std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";
		fnout = fnstem + sres + ".clusters";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		lp.updateSoFarLazy(GsoFar);
		fnout = fnstem +  "." + sres + "soFar";
		if (doScores) 
			Glabel->putSoFar(GsoFar, residint+1);
		GsoFar->writeSingle(fnout.c_str());
		delete[] Gnew;
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
	delete en;
	//delete G;
	return 1;
}
