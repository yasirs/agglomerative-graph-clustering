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
	linkPredictor lp;
	char dummy[100];
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
	// std::cin >> dummy;
	fnstem = argv[1];
	//fnstem = "grassland";
	fnin = fnstem+".edges";
	// let us read in the graph
	graphData *Goriginal, *Gnew, *GPred, *GOld, *GsoFar;
	Engine *en;
	Goriginal = new graphData[1];
	Goriginal[0].readBinary(fnin.c_str());
	std::cout << "read file!\n";
	en = new Engine(Goriginal,1);
	en->initializeScores();
	// let us write the competing scores
	/*std::cout << "writing degree product\n";
	fnout = fnstem + ".dprod"; en->printDegreeProdFile(fnout.c_str(), 0, 0);
	std::cout << "writing hyper geom\n";
	fnout = fnstem + ".hyperg"; en->printHyperGeomFile(fnout.c_str(), 0, 0);
	std::cout << "writing jaccard\n";
	fnout = fnstem + ".jaccard"; en->printJaccardFile(fnout.c_str(), 0, 0);
	std::cout << "writing comm neigh\n";
	fnout = fnstem + ".cneighb"; en->printCommonNeighbFile(fnout.c_str(), 0, 0);
	*/

	// run the agglomerative algorithm
	std::cout << "starting to run\n";
	en->run();
	fnout = fnstem + ".network";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	std::cout << "done running\n";
	lp.attach(en); // attached the enginea

	/*GPred = lp.makeNonEdgePred(Goriginal);
	fnout = fnstem + ".scores0";
	GPred[0].writeSingle(fnout.c_str());
	delete[] GPred;*/

	GsoFar = lp.copyNoEdges(Goriginal);
	//fnout = fnstem + ".nodes";
	//en->tree->writeNodeTypes(fnout.c_str());
	std::cout << "getting the residual graph\n";
	//Gnew = getResidual(G, lp);
	Gnew = getORResidual(Goriginal, lp, GsoFar);
	
	fnout = fnstem + "0.residual";
	Gnew->writeSingle(fnout.c_str());
	fnout = fnstem + "0.clusters";
	en->tree->writeCollapsedHierEdges(fnout.c_str());
	updateSoFar(Goriginal, lp, GsoFar);
	std::cout << "Etot = "<< Gnew->Etot << "\n";
	int residint = 1;
	while(residint <= numResids) {
		std::cout << "got the residual\n";
		delete en;
		en = new Engine(Gnew,1);
		std::cout << "running on the residual\n";
		en->initializeScores();
		en->run();
		lp.attach(en);

		/*GPred = lp.makeNonEdgePred(Goriginal);
		std::ostringstream os;
		os << fnstem << ".scores" << residint;
		fnout = os.str();
		GPred[0].writeSingle(fnout.c_str());
		delete[] GPred;*/


		//fnout = fnout + "_noname";
		//GPred[0].writeSingle_noname(fnout.c_str());
		std::cout << "done!\nPrinting out the " << residint <<"th heirarchical network\n";
		std::stringstream dummy;
		dummy << residint;
		std::string sres;
		dummy >> sres;
		fnout = fnstem + sres + ".residual";
		GOld = Gnew;
		std::cout << "getting the residual graph\n";
		Gnew = getORResidual(Goriginal, lp, GsoFar);
		Gnew->writeSingle(fnout.c_str());
		fnout = fnstem + sres + ".clusters";
		en->tree->writeCollapsedHierEdges(fnout.c_str());
		updateSoFar(Goriginal, lp, GsoFar);
		delete[] GOld;
		std::cout << "Etot = "<< Gnew->Etot << "\n";
		residint++;
	}
	std::cout << "done!\n";
	delete[] Goriginal;
	delete en;
	//delete G;
	return 1;
}
