#ifndef GRAPHUTILS_INCLUDED
#define GRAPHUTILS_INCLUDED
void importGrid(char* file, int **ia, int **ja, int *nv);
void printGraph(int* ia, int *ja, int nv);
void printMappedGraph(int* ia, int *ja,int *map, int nv);
int procNrFromCell(int LocalGridDims[], int CartGridDims[], int CellNr);
#endif