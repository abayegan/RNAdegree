# Makefile for computing number of neighbors, both with and without shift
# P.Clote, Dec 2014

all:	uniformDegree
#all:	nbor0 
#all:	nbor1

uniformDegree: uniformExpNumNbors.o numNborAddRem.o numNborAddRemShift.o aux.o
	gcc uniformExpNumNbors.o numNborAddRem.o numNborAddRemShift.o aux.o -o uniformDegree
	
uniformExpNumNbors.o : uniformExpNumNbors.c
	gcc -c uniformExpNumNbors.c

numNborAddRem.o: numNborAddRem.c
	gcc -c numNborAddRem.c

numNborAddRemShift.o: numNborAddRemShift.c
	gcc -c numNborAddRemShift.c

aux.o:	aux.c
	gcc -c aux.c

clean:
	rm -fr *.o uniformDegree

