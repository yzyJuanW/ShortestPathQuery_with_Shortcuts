#
#	Makefile 
#	Author: zy based on Trasier
##
CC = gcc
CXX = g++
#yuxiang's method
LIBS = -D HELLINGER

#METISLIB = -I /homes/zgongae/metis-5.1.0/include -L /homes/zgongae/metis-5.1.0/build/Linux-x86_64/libmetis -lmetis
#METISLIB = -I /homes/zgongae/PGM-index/include

CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium
MEM = -D WATCH_MEM

all: tdtree

# tdtree: mainTD.cpp treedecomp.h PLF.h Constants.h Segment.h misc.h
# 	$(CXX) $(CFLAGS) $(METISLIB) -o tdtree mainTD.cpp


tdtree: mainTD.cpp Shortcuts_.h PLF.h Constants.h Segment.h misc.h
	$(CXX) $(CFLAGS) -o tdtree mainTD.cpp $(LIBS)


# treedecomp.o: treedecomp.h PLF.h Constants.h Segment.h misc.h 
# 	$(CXX) $(CFLAGS) -c treedecomp.h $(LIBS)
# tdtree: mainTD.cpp treedecomp.o PLF.h Constants.h Segment.h misc.h
# 	$(CXX) $(CFLAGS) -o n3 treedecomp.o mainTD.cpp $(LIBS)






# insertionN2.o: TDGT.h io.h insertionN2.h insertionN2.cpp PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h 
# 	$(CXX) $(CFLAGS) -c insertionN2.cpp $(LIBS)$(METISLIB)
# n2: mainN2.cpp insertionN2.h insertionN2.cpp insertionN2.o TDGT.h io.h PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h TDDijkstra.h 
# 	$(CXX) $(CFLAGS) -o n2 insertionN2.o mainN2.cpp $(LIBS) $(METISLIB)

# insertionN.o: TDGT.h io.h insertionN.h insertionN.cpp PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h 
# 	$(CXX) $(CFLAGS) -c insertionN.cpp $(LIBS)$(METISLIB)
# n: mainN.cpp insertionN.h insertionN.cpp insertionN.o TDGT.h io.h PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h TDDijkstra.h 
# 	$(CXX) $(CFLAGS) -o n insertionN.o mainN.cpp $(LIBS) $(METISLIB)


.PHONY: clean
clean:
		-@rm *.o *.gcno *~ 2> /dev/null || true
		-@rm tdtree 2> /dev/null || true
##############################################################################################################
# ##
# #	Makefile 
# #	Author: zy based on Trasier
# ##
# CC = gcc
# CXX = g++
# #yuxiang's method
# LIBS = -D HELLINGER

# METISLIB = -I /homes/zgongae/metis-5.1.0/include -L /homes/zgongae/metis-5.1.0/build/Linux-x86_64/libmetis -lmetis
# #-L /Users/gongcengyang/Desktop/TestTDGT/metis-5.1.0/build/Darwin-x86_64/libmetis -lmetis

# CFLAGS = --std=c++11 -W -Wall -Wno-sign-compare -O3 -s -pipe -mmmx -msse -msse2 -msse3 # -mcmodel=medium
# MEM = -D WATCH_MEM

# all: TDGT

# insertion.o: TDGT.h io.h insertion.h insertion.cpp PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h 
# 	$(CXX) $(CFLAGS) -c insertion.cpp $(LIBS)$(METISLIB)
# TDGT: main.cpp insertion.h insertion.cpp insertion.o TDGT.h io.h PLF.h Constants.h TGraph.h TNode.h Segment.h misc.h TMatrix.h TDDijkstra.h 
# 	$(CXX) $(CFLAGS) -o TDGT insertion.o main.cpp $(LIBS) $(METISLIB)


# .PHONY: clean
# clean:
# 		-@rm *.o *.gcno *~ 2> /dev/null || true
# 		-@rm TDGT 2> /dev/null || true
