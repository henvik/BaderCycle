#ifndef GRAPHUTILS_INCLUDED
#define GRAPHUTILS_INCLUDED

typedef struct {
   int    v;
   int    w;
}Edge;

typedef struct {
   int size;
   int* list;
}AdjLst;

typedef struct {
   int size;
   Edge *list;
}EdgeLst;


void importGrid(char* file, int **ia, int **ja, int *nv);
void printGraph(int* ia, int *ja, int nv);
void printMappedGraph(int* ia, int *ja,int *map, int nv);
int procNrFromCell(int LocalGridDims[], int CartGridDims[], int CellNr);
void globalCoords2localCoords(int globalCoords[],int localDims[],  int* localCoords );
int  localCoords2localCellNr(int localCoords[], int localDims[]);
void cellNr2cartCoord(int cellNr, int* GlobalDims, int* output);
int cartCoord2cellNr( int x, int y, int* dims);
int globalCellNr2localCellNr(int globalCellNr, int globalDims[],int localDims[]);


void concatenate(int* A, int* B, int* lengthA, int* lengthB);

Edge new_edge(int v, int w);
void add_trans_arc(int v, int w,EdgeLst *trans_arcs, int *trans_arcs_count,AdjLst* adjecent);
AdjLst new_AdjLst();
AdjLst merge_AdjLst(AdjLst A, AdjLst B);
void add_Edge(AdjLst *A, int edge);


EdgeLst new_EdgeLst(int size);
void add_EdgeToLst(Edge e,int index, EdgeLst *lst);

void add_edges(EdgeLst *ExpGraph, EdgeLst *recv_buffer);
void expand_buffer(EdgeLst *recv_buffer, int recv_size);

#endif
