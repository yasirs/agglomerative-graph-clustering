default: scoreSome.cpp ../inc/nodetree.hpp ../inc/nodetreeOther.hpp ../inc/dataMap.hpp ../inc/dataMapOther.hpp \
	../inc/mysetfuncs.hpp ../inc/mygenerators.hpp ../inc/linkPredictor.hpp ../inc/linkPredictorOther.hpp \
	../inc/residuals.hpp ../inc/scoremap.hpp ../inc/modelStats.hpp ../inc/modelParam.hpp
	g++ -I../inc/ /opt/local/lib/libgsl.a scoreSome.cpp -o ~/temp/scoreSome -O3

debug: scoreSome.cpp nodetree.hpp nodetreeOther.hpp dataMap.hpp dataMapOther.hpp \
	mysetfuncs.hpp mygenerators.hpp linkPredictor.hpp linkPredictorOther.hpp \
	residuals.hpp scoremap.hpp modelStats.hpp modelParam.hpp
	/usr/bin/g++-4.2 -lgsl -L/opt/local/lib scoreSome.cpp -o ~/temp/scoreSome -g -Wall


holdOutCV: holdOutCV.cpp ../inc/mysetfuncs.hpp 
	g++ -I../inc/ /opt/local/lib/libgsl.a holdOutCV.cpp -o ~/temp/holdOutCV -O3



writeLabels: writeLabels.cpp ../inc/mysetfuncs.hpp 
	g++ -I../inc/ /opt/local/lib/libgsl.a writeLabels.cpp -o ~/temp/writeLabels -O3



