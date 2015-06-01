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



typedef struct ExVert{
	int vert_num;
	int procNr;
	
	struct TransArc *trans_arcs;
	struct ExArc *exp_arcs;
	
	struct ExVert *next;
}ExVert;



typedef struct TransArc{
	int proc_nr;
	int vert_num;
	struct TransArc *next;
}TransArc;

typedef struct ExArc{
	int  vert_num;
	int proc_nr;
	struct ExArc *next;
}ExArc;

typedef struct ExpGraph{
	int num_vert;
	ExVert* first;
}ExpGraph;



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
void merge_AdjLst(AdjLst *A, AdjLst *B);
void add_Edge(AdjLst *A, int edge);
bool isInAdjLst(AdjLst A,int x);


EdgeLst new_EdgeLst(int size);
void add_EdgeToLst(Edge e,int index, EdgeLst *lst);

void add_edges(EdgeLst *ExpGraph, EdgeLst *recv_buffer);
void expand_buffer(EdgeLst *recv_buffer, int recv_size);


//------------------------ExpGraph

ExpGraph* new_ExpGraph();

ExVert* new_ExVert(int vert_num, int procNr);

void addExVert(ExpGraph *G, ExVert *vert);

void addTransArc(ExVert *vert, TransArc *trans_arc);

void addExArc(ExVert *vert, ExArc *ex_arc);
TransArc* newTransArc(int proc_nr, int vert_num);

ExArc* newExArc(int vert_num, int proc_nr);
bool isInExp_arcs(ExVert *vert,ExArc *ex_arc);


ExVert* FindExVert(ExpGraph *G, int v);

void printfExpGraph(ExpGraph *G);

void addExternalEdge(ExpGraph *G, int internal_vert_num, int external_vert_num, int in_procNr ,int external_procNr, AdjLst *adjecent);

bool isReachable(AdjLst adjecent,int vert_num);

ExpGraph* mergeVertices(ExpGraph *exp1,ExpGraph *exp2);

void removeExVert(ExpGraph *G,ExVert *vert);
void removeExArc(ExVert *vert, ExArc *exarc);

#endif
