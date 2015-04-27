#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include"globalVars.h"
#include"baderDetect.h"
// MPI initialization, setting up cartesian communicator
void init_mpi(int argc, char** argv){
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm );
    MPI_Cart_coords( cart_comm, rank, 2, coords );
    
    MPI_Cart_shift( cart_comm, 0, 1, &north, &south );
    MPI_Cart_shift( cart_comm, 1, 1, &west, &east );
}


void def_datatypes(MPI_Datatype *particletype){
    MPI_Datatype oldtypes[1];
    int count,  blockcount[1];
    MPI_Aint  offsets[1];
    count=1;
    blockcount[0]=2;
    oldtypes[0]=MPI_INT;
    offsets[0]=0;
    MPI_Type_struct(count,blockcount,offsets,oldtypes,particletype);
    MPI_Type_commit(particletype);
}
