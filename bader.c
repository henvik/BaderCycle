#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<stdbool.h>
#include<time.h>

#include"DistGather.h"
#include"GraphUtilities.h"
#include"mpiUtils.h"
#include"baderDetect.h"

#define MAIN_FILE
#include "globalVars.h"




int main(int argc,  char* argv[]){

	init_mpi(argc,argv);
//	printf("I am process %d of %d\n", rank,size);
	int nv;
	int *ia, *ja;
	
	if(rank==0){
			importGrid("Grids/graphSpec512x512.txt", &ia, &ja, &nv);
//			printGraph(ia,ja,nv);
		
	}
				
	//All processors need to know the global problem size. (Would be cool, and probably possible, to avoid this.) 
 	MPI_Bcast(&nv, 1, MPI_INT, 0,cart_comm);
 	localSetup(nv);
	distGraph(ia,ja);
// 	printf("I am nr %d. North=%d, south=%d, west=%d, and east=%d. \n",rank,north,south,west,east);	
	
	begin=clock();
	ExpGraph *expGraph=discovery(nv/size);
	//printfExpGraph(expGraph);
// 	printMappedGraph(local_ia,local_ja,local_map,nv/size);
  
  	merge(expGraph);
	end=clock();
	time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
	printf("For proc %d, BaderDetect took %f \n",rank,time_spent);
    free_local();
    if(rank==0){
	    free(ia);
	    free(ja);
    }
    MPI_Finalize();

	return 0;
}
