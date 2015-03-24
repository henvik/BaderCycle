#ifndef BADERDETECT_INCLUDED
#define BADERDETECT_INCLUDED

typedef struct {
   int    v;
   int    w;
}Edge;

void discovery(int num_vert);
int* visit(int v);
void concatenate(int* A, int* B, int* lengthA, int* lengthB);

Edge new_edge(int v, int w);
void add_trans_arc(int v, int w,Edge* trans_arcs, int trans_arcs_count);
#endif