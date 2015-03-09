#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include"DistGather.h"
#include"GraphUtilities.h"
#include"mpiUtils.h"
#include"baderDetect.h"

#define MAIN_FILE
#include "globalVars.h"

int main(int argc,  char* argv[]){

	init_mpi(argc,argv);
	printf("I am process %d of %d\n", rank,size);
	int nv;
	int *ia, *ja;
	if(rank==0){
			int i,k;
			importGrid("Grids/grid1src100.txt", &ia, &ja, &nv);
	//		printGraph(ia,ja,nv);
	}
	//All processors need to know the global problem size. (Would be cool, and probably possible, to avoid this.) 
	MPI_Bcast(&nv, 1, MPI_INT, 0,cart_comm);
	localSetup(nv);
	distGraph(ia,ja);
 	if(rank==0){
 	printMappedGraph(local_ia,local_ja,local_map,nv/size);
 	}

	int *vert= (int *)malloc((50+1)*sizeof(int)); 
	int *comp= (int *)malloc((50+1)*sizeof(int));
	int *ncomp;
	int* work =(int *)malloc(3*50*sizeof(int));

    MPI_Finalize();

	return 0;
}
