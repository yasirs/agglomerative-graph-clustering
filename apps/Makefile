PATHTOLIBS = /opt/local/lib
MYGPP = /opt/local/bin/g++-mp-4.4
NATIVEGPP = /usr/bin/g++-4.2
OUTPATH = ~/temp

default: scoreSome.cpp ../inc/nodetree.hpp ../inc/nodetreeOther.hpp ../inc/dataMap.hpp ../inc/dataMapOther.hpp \
	../inc/mysetfuncs.hpp ../inc/mygenerators.hpp ../inc/linkPredictor.hpp ../inc/linkPredictorOther.hpp \
	../inc/residuals.hpp ../inc/scoremap.hpp ../inc/modelStats.hpp ../inc/modelParam.hpp
	$(MYGPP) -I../inc/ $(PATHTOLIBS)/libgsl.a scoreSome.cpp -o $(OUTPATH)/scoreSome -O3

debug: scoreSome.cpp ../inc/nodetree.hpp ../inc/nodetreeOther.hpp ../inc/dataMap.hpp ../inc/dataMapOther.hpp \
	../inc/mysetfuncs.hpp ../inc/mygenerators.hpp ../inc/linkPredictor.hpp ../inc/linkPredictorOther.hpp \
	../inc/residuals.hpp ../inc/scoremap.hpp ../inc/modelStats.hpp ../inc/modelParam.hpp
	$(NATIVEGPP) -I../inc/ $(PATHTOLIBS)/libgsl.a scoreSome.cpp -o $(OUTPATH)/scoreSome -g -Wall

holdOutCV: holdOutCV.cpp ../inc/mysetfuncs.hpp 
	$(MYGPP) -I../inc/ holdOutCV.cpp -o $(OUTPATH)/holdOutCV -O3 -Wall

writeLabels: writeLabels.cpp ../inc/mysetfuncs.hpp 
	$(MYGPP) -I../inc/ /opt/local/lib/gcc44/libstdc++.a writeLabels.cpp -o $(OUTPATH)/writeLabels -O3 -Wall


all: default holdOutCV writeLabels

