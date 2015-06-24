#ifndef BADERDETECT_INCLUDED
#define BADERDETECT_INCLUDED

#include"GraphUtilities.h"

typedef struct RemoveList{
	int capasity;
	int size;
	ExVert **list;

} RemoveList;


//The first phase of the algorithm, finds local cycles, builds the local part of the express graph, and calls the communication step 
ExpGraph* discovery(int num_vert);


//Visits a node. Calls recursively
AdjLst  visit(int v);

//Communicates transArcs and builds the express graphs
void comm_transArcs(ExpGraph *expGraph);

// Builds the buffer for sending the transarcs
void build_transBuffer(int *buffSizes, EdgeLst *trans_buffer, ExpGraph *expGraph);


//Merges the local sub-graphs together into an single graphs, detecting cycles as it goes.
void merge(ExpGraph* expGraph);

int last(int z, int h);

int test(int z, int h);

int clear(int z, int h);
int set(int z, int h);

void completeExpGraph(ExpGraph *G, AdjLst *adjecent);

void recvTransArcs(Edge* list, int size, ExpGraph *expGraph, int procNr);
	

void RecieveExp(int procNr, ExpGraph** graphBuff);

void SendExp(int procNr, ExpGraph* expGraph);

ExpGraph* MergeGraphs(ExpGraph* exp1, ExpGraph* exp2, int origin1, int origin2);

void addToBuffer(int **send_buffer,int toBeAdded ,int *counter, int *size_send_buffer);

int packGraph(ExpGraph *exp,int **send_buffer);

ExpGraph *unPackRecvBuffer(int *recv_buffer);

void printBuffer(int **send_buffer,int counter);

void checkVertex(ExpGraph *exp0,ExVert *xExVert, int origin, RemoveList *toRemove);

bool isExArc(ExVert *terminal ,ExVert *initial);

void transferExArcs(ExpGraph *G, ExVert *globalCurrent, ExVert *init, ExVert *term);

void removeAndReplaceTransArc(ExpGraph *G,ExVert *init, ExVert *term, RemoveList *toRemove );

void removeTransArcs(ExVert *vert,ExVert *remove);

void addToListOfMarked(RemoveList *list,ExVert *add);

RemoveList* new_RemoveList();



#endif
