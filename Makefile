default: scoreSome.cpp nodetree.hpp nodetreeOther.hpp dataMap.hpp dataMapOther.hpp \
	mysetfuncs.hpp mygenerators.hpp linkPredictor.hpp linkPredictorOther.hpp residuals.hpp scoremap.hpp
	g++ -lgsl -L/opt/local/lib scoreSome.cpp -o ~/temp/scoreSome -O3




debug: scoreSome.cpp nodetree.hpp nodetreeOther.hpp dataMap.hpp dataMapOther.hpp \
	mysetfuncs.hpp mygenerators.hpp linkPredictor.hpp linkPredictorOther.hpp residuals.hpp scoremap.hpp
	/usr/bin/g++-4.2 -lgsl -L/opt/local/lib scoreSome.cpp -o ~/temp/scoreSome -g -Wall
