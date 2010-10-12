default: scoreSome.cpp ../inc/nodetree.hpp ../inc/nodetreeOther.hpp ../inc/dataMap.hpp ../inc/dataMapOther.hpp \
	../inc/mysetfuncs.hpp ../inc/mygenerators.hpp ../inc/linkPredictor.hpp ../inc/linkPredictorOther.hpp \
	../inc/residuals.hpp ../inc/scoremap.hpp ../inc/modelStats.hpp ../inc/modelParam.hpp
	g++ -I../inc/ /opt/local/lib/libgsl.a scoreSome.cpp -o ~/temp/scoreSome -O3

debug: scoreSome.cpp ../inc/nodetree.hpp ../inc/nodetreeOther.hpp ../inc/dataMap.hpp ../inc/dataMapOther.hpp \
	../inc/mysetfuncs.hpp ../inc/mygenerators.hpp ../inc/linkPredictor.hpp ../inc/linkPredictorOther.hpp \
	../inc/residuals.hpp ../inc/scoremap.hpp ../inc/modelStats.hpp ../inc/modelParam.hpp
	g++ -I../inc/ -L/opt/local/lib -lgsl scoreSome.cpp -o ~/temp/scoreSome -g -Wall

holdOutCV: holdOutCV.cpp ../inc/mysetfuncs.hpp 
	g++ -I../inc/ holdOutCV.cpp -o ~/temp/holdOutCV -O3 -Wall

writeLabels: writeLabels.cpp ../inc/mysetfuncs.hpp 
	g++ -I../inc/ writeLabels.cpp -o ~/temp/writeLabels -O3 -Wall



