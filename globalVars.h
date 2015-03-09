#ifndef SHAREFILE_INCLUDED
#define SHAREFILE_INCLUDED
#ifdef  MAIN_FILE
//Global variables
int rank,                       // MPI rank
    size,                       // Number of MPI processes
    local_numOfEdges,			//number of edges in the local grid
    dims[2],                    // Dimensions of MPI grid
    coords[2],                  // Coordinate of this rank in MPI grid
    periods[2] = {0,0},         // Periodicity of grid
	gridDims[2],
	local_dims[2],
    north,south,east,west;      // Four neighbouring MPI ranks

MPI_Comm cart_comm; 

int* local_ia;
int* local_ja;
int* local_map;
#else
//Global variables
extern int rank,                       // MPI rank
    size,                       // Number of MPI processes
    local_numOfEdges, 			//number of edges in the local grid
    dims[2],                    // Dimensions of MPI grid
    coords[2],                  // Coordinate of this rank in MPI grid
    periods[2],        // Periodicity of grid
    gridDims[2],
    local_dims[2],
    north,south,east,west;      // Four neighbouring MPI ranks

extern MPI_Comm cart_comm; 
extern int* local_ia;
extern	int* local_ja;
extern	int* local_map;
#endif
#endif


