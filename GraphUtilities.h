#ifndef GRAPHUTILS_INCLUDED
#define GRAPHUTILS_INCLUDED
void importGrid(char* file, int **ia, int **ja, int *nv);
void printGraph(int* ia, int *ja, int nv);
void printMappedGraph(int* ia, int *ja,int *map, int nv);
int procNrFromCell(int LocalGridDims[], int CartGridDims[], int CellNr);
void globalCoords2localCoords(int globalCoords[],int localDims[],  int* localCoords );
int  localCoords2localCellNr(int localCoords[], int localDims[]);
void cellNr2cartCoord(int cellNr, int* GlobalDims, int* output);
int cartCoord2cellNr( int x, int y, int* dims);
int globalCellNr2localCellNr(int globalCellNr, int globalDims[],int localDims[]);
#endif
