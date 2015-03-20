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
	
	int dimsTest[2]={16,16};
	int testCell=cartCoord2cellNr( 6, 7, dimsTest);
	printf("this: %d \n",testCell);
	
	
	if(rank==0){
			importGrid("Grids/grid1src100.txt", &ia, &ja, &nv);
			//printGraph(ia,ja,nv);


	}
	
	//All processors need to know the global problem size. (Would be cool, and probably possible, to avoid this.) 
 	MPI_Bcast(&nv, 1, MPI_INT, 0,cart_comm);
// 	printf("Here\n");
 	localSetup(nv);
	if(rank==0){
 	printf("Local dims is : (%d,%d)\n",local_dims[0],local_dims[1]);
 	int whatIsRank;
	whatIsRank=procNrFromCell(local_dims, gridDims, 75);
	printf("Cell nr 25 resides on processor: %d\n", whatIsRank);
	}
// 	distGraph(ia,ja);
// 	if(rank==0){
 //	printMappedGraph(local_ia,local_ja,local_map,nv/size);
//  	}
	
	
    MPI_Finalize();

	return 0;
}
