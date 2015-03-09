#ifndef DISTGATHER_INCLUDED
#define DISTGATHER_INCLUDED
void localSetup(int nv);
void cellNr2cartCoord(int cellNr, int** coords, int* dims);
int cartCoord2cellNr( int x, int y, int* dims);
void buildSubGraph(int i, int* local_ia, int* local_ja, int* local_map, int* j_count_out);
void distGraph(int* ia, int* ja);
#endif
