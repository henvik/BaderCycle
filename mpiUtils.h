#ifndef MPIUTILS_INCLUDED
#define MPIUTILS_INCLUDED
// MPI initialization, setting up cartesian communicator
void init_mpi(int argc, char** argv);
void def_datatypes(MPI_Datatype *particletype);
#endif
