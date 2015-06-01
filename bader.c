#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<stdbool.h>

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
			importGrid("Grids/grid1src100.txt", &ia, &ja, &nv);
		//	printGraph(ia,ja,nv);
		
	}
		
// 	 int num=1;
// 	int bit=0;
// // 	
//  	printf("We set the %dth bit of %d to 1: %d\n",bit,num,set(num,bit));
		
	//All processors need to know the global problem size. (Would be cool, and probably possible, to avoid this.) 
 	MPI_Bcast(&nv, 1, MPI_INT, 0,cart_comm);
 	localSetup(nv);
	distGraph(ia,ja);
// 	printf("I am nr %d. North=%d, south=%d, west=%d, and east=%d. \n",rank,north,south,west,east);	

 //	if(rank==3){
	ExpGraph *expGraph=discovery(nv/size);
// 	printMappedGraph(local_ia,local_ja,local_map,nv/size);
  //	}
	merge(expGraph);
    free_local();
    if(rank==0){
	    free(ia);
	    free(ja);
    }
    MPI_Finalize();

	return 0;
}
