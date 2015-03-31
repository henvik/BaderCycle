#ifndef BADERDETECT_INCLUDED
#define BADERDETECT_INCLUDED

typedef struct {
   int    v;
   int    w;
}Edge;

typedef struct {
   int size;
   int* list;
}AdjLst;

void discovery(int num_vert);
AdjLst  visit(int v);
void concatenate(int* A, int* B, int* lengthA, int* lengthB);

Edge new_edge(int v, int w);
void add_trans_arc(int v, int w,Edge* trans_arcs, int trans_arcs_count);
AdjLst new_AdjLst();
AdjLst merge_AdjLst(AdjLst A, AdjLst B);
AdjLst add_Edge(AdjLst A, int edge);
#endif
