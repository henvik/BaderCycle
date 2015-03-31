#ifndef DISTGATHER_INCLUDED
#define DISTGATHER_INCLUDED
void localSetup(int nv);
void buildSubGraph(int i, int* local_ia, int* local_ja, int* local_map, int* j_count_out);
void distGraph(int* ia, int* ja);

#endif
