#define DEBUGMODE 1
#define NOREFERENCE 0
#include "Engine.hpp"
#include "graphData.hpp"
#include "residuals.hpp"
#include "linkPredictor.hpp"
#include "linkPredictorOther.hpp"
#include "labelData.hpp"
#include <iostream>
#include <string> 
#include <stdexcept>


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
		if (gtype!='b') {
			std::invalid_argument("Only binary graphs for now!\n");
		}
	} else {
		gtype = 'b';
		std::cout << "graph type = " << gtype << '\n';
	}
	fnstem = argv[1];
	fnin = fnstem+".edges";
	graphData *Goriginal, *Gnew, *GsoFar, *Gempty;
	Engine *en;
	Goriginal = new graphData[2];
	GsoFar = new graphData[2];
	Goriginal[0].readGeneral(fnin.c_str());
	Goriginal[0].gtype = gtype;
	Goriginal[0].makeComplimentary(Goriginal[1]); Goriginal[1].gtype = gtype;
	std::cout << "read file!\n";
	std::cout << "Initial Etot = "<<Goriginal->Etot << "\n";

	Goriginal[0].copyNoEdges(GsoFar[0]);
	Goriginal[0].copyNoEdges(GsoFar[1]);

	std::cout << "Initializing Engine, scores\n";
	en = new Engine(Goriginal, Goriginal, GsoFar, 2);
	en->initializeScoresML();

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
	lp.updateTrackEdgeHoles(GsoFar);
	fnout = fnstem + ".0soFarEdges"; GsoFar[0].writeSingle(fnout.c_str());
	fnout = fnstem + ".0soFarHoles"; GsoFar[1].writeSingle(fnout.c_str());


	// write the post-pass predictions for this round
	fnout = fnstem+"." + "0postpassEdges"; 
	Gempty = new graphData[2];
	Goriginal[0].copyNoEdges(Gempty[0]);
	Goriginal[1].copyNoEdges(Gempty[1]);	
	lp.updateSoFarLazy(Gempty); Gempty[0].writeSingle(fnout.c_str()); 
	fnout = fnstem+"."+"0postpassHoles"; Gempty[1].writeSingle(fnout.c_str());
	delete[] Gempty;
	


	int residint = 1;
	while(residint <= numResids) {
		std::stringstream dummy;
		dummy << residint;
		std::string sres;
		dummy >> sres;
		std::cout << "getting the residual graph\n";
		Gnew = residualDiff(Goriginal, GsoFar, 2);
		fnout = fnstem + sres + ".residualEdges"; Gnew[0].writeSingle(fnout.c_str());
		fnout = fnstem + sres + ".residualHoles"; Gnew[1].writeSingle(fnout.c_str());
		std::cout << "got the residual\n";
		std::cout << "Etot = "<< Gnew->Etot << "\n";
		delete en;
		en = new Engine(Gnew,Goriginal,GsoFar,2);
		en->initializeScoresML();
		std::cout << "running on the residual\n";
		en->runML(true,false);
		fnout = fnstem +"_"+int2str(residint)+".hrg"; en->printHRG(fnout.c_str(),0);
		
		en->passFB();
		lp.attach(en);

		
		// write the post-pass predictions for this round
		fnout = fnstem+"." + sres+"postpassEdges"; 
		Gempty = new graphData[2];
		Goriginal[0].copyNoEdges(Gempty[0]);
		Goriginal[1].copyNoEdges(Gempty[1]);	
		lp.updateSoFarLazy(Gempty); Gempty[0].writeSingle(fnout.c_str()); 
		fnout = fnstem+"."+sres+"postpassHoles"; Gempty[1].writeSingle(fnout.c_str());
		delete[] Gempty;
		

		//fnout = fnstem + ".scores" + sres;
		//std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";
		fnout = fnstem + sres + ".clusters";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		lp.updateTrackEdgeHoles(GsoFar);
		fnout = fnstem +  "." + sres + "soFarEdges"; GsoFar[0].writeSingle(fnout.c_str());
		fnout = fnstem +  "." + sres + "soFarHoles"; GsoFar[1].writeSingle(fnout.c_str());
		delete[] Gnew;
		residint++;
	}
	std::cout << "done!\n";
	delete[] Goriginal;
	delete[] GsoFar;
	delete en;
	//delete G;
	return 1;
}
