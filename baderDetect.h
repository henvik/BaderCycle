#ifndef BADERDETECT_INCLUDED
#define BADERDETECT_INCLUDED

#include"GraphUtilities.h"

//The first phase of the algorithm, finds local cycles, builds the local part of the express graph, and calls the communication step 
void discovery(int num_vert);


//Visits a node. Calls recursively
AdjLst  visit(int v);

//Communicates transArcs and builds the express graphs
void comm_transArcs(EdgeLst *expGraph);

// Builds the buffer for sending the transarcs
void build_transBuffer(int *buffSizes, EdgeLst *trans_buffer);


//Merges the local sub-graphs together into an single graphs, detecting cycles as it goes.
void merge(EdgeLst* expGraph);

int last(int z, int h);

int test(int z, int h);

int clear(int z, int h);
int set(int z, int h);

void RecieveExp(int procNr, EdgeLst* graphBuff);

void SendExp(int procNr, EdgeLst* graphBuff);

void MergeGraphsEx(EdgeLst* expGraph, EdgeLst* graphBuff);

#endif
