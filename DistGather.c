#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include"globalVars.h"


void localSetup(int nv){
	gridDims[0]=sqrt(nv);
	gridDims[1]=sqrt(nv);		
	local_dims[0]=gridDims[0]/dims[0];
	local_dims[1]=gridDims[1]/dims[1];
	
	local_numOfEdges=local_dims[0]*local_dims[1]*4;  //preliminary value
	local_ia= (int *)malloc(local_dims[0]*local_dims[1]*sizeof(int));
	local_ja= (int *)malloc(4*local_dims[0]*local_dims[1]*sizeof(int)); //4 times the number of nodes should be enough (allways for 2D graph). Reallocating later if necassary.
	local_map= (int *)malloc(local_dims[0]*local_dims[1]*sizeof(int));
}

//Makes cartesian 
void cellNr2cartCoord(int cellNr, int** coords, int* dims){
/* Input: 	cellNr  - the cellNr in the grid 
		  	dims    - an array specifying the dimensions of the cartesian grid
	Output: coords  - the cartesian coordinates of the cellNR
*/	
	int* output=(int *)malloc(sizeof(int)*2);
 	output[1]=cellNr%dims[1];
 	output[0]=cellNr/dims[0];	
 	*coords=output;
}

int cartCoord2cellNr( int x, int y, int* dims){
/* Input: 	x,y  - the cartesian coordinates of the cellNR
		  	dims    - an array specifying the dimensions of the cartesian grid
	Output: cellNr  - the cellNr in the grid 
*/	
	int output=x*dims[1]+y;
	return output;
}

void buildSubGraph(int i,int* ia, int* ja, int* local_ia, int* local_ja, int* local_map, int* j_count_out){
/* Input: i - rank of the receiving process
		  ia,ja - pointers to the entire graph
		  j_count - number of elements to be sent.
		  local_ia, local_ja, local_map - pointers to the buffers in which to store the subgraph
*/

		int j,k,l,tmp,count,j_count,
			dest[2],
			xlimits[2],
			ylimits[2];

			MPI_Cart_coords(cart_comm,i, 2, dest);
		//	printf("Dest: (%d,%d), with local_dims: (%d,%d)\n",dest[0],dest[1],local_dims[0],local_dims[1]);
			xlimits[0]=dest[0]*local_dims[0];
			xlimits[1]=xlimits[0]+local_dims[0];
			ylimits[0]=dest[1]*local_dims[1];
			ylimits[1]=ylimits[0]+local_dims[1];
		//	printf("xlims: (%d,%d), ylims: (%d,%d) \n",xlimits[0],xlimits[1],ylimits[0],ylimits[1]);
			count=0;
			j_count=0;
			for(j=xlimits[0];j<xlimits[1];j++){
				for(k=ylimits[0];k<ylimits[1];k++){
		//			printf("gridDims is (%d,%d)\n",gridDims[0],gridDims[1]);
					tmp=cartCoord2cellNr(j,k,gridDims);
		//			printf("tmp is %d for (%d,%d) \n ",tmp,j,k);
					local_map[count]=tmp;
					local_ia[count++]=j_count;				
					for(l=ia[tmp];l<ia[tmp+1];l++){
						if(local_numOfEdges<j_count){
							local_numOfEdges=local_numOfEdges+local_dims[0]*local_dims[1];
							local_ja = (int *) realloc(local_ja,local_numOfEdges);
						}
						local_ja[j_count++]=ja[l];
					}
				}	
			}
	    	local_ia[count]=j_count;
	    	*j_count_out=j_count;
}

//Takes in a sparse grid, including sources, and distributes it to n different processes
void distGraph(int* ia, int* ja){
	if(rank==0){
	
		int i,*j_count;
		
		for(i=1;i<size;i++){
			buildSubGraph(i, ia, ja, local_ia,local_ja,local_map, j_count);
			MPI_Send(local_ia, local_dims[0]*local_dims[1], MPI_INT, i, 0,cart_comm);
			MPI_Send(local_map, local_dims[0]*local_dims[1], MPI_INT, i, 1, cart_comm);
			MPI_Send(local_ja,*j_count , MPI_INT, i, 2, cart_comm);
		}
	}
	else{
	 int size_ia,size_ja,size_map;
	MPI_Status status_ia,status_ja,status_map;
	MPI_Probe(0, 0, cart_comm, &status_ia);
	MPI_Get_count(&status_ia, MPI_INT, &size_ia);
	MPI_Recv(local_ia, size_ia, MPI_INT, 0, 0,cart_comm, MPI_STATUS_IGNORE);
	
	MPI_Probe(0, 1, cart_comm, &status_map);
	MPI_Get_count(&status_map, MPI_INT, &size_map);
	MPI_Recv(local_map, size_map, MPI_INT, 0, 1,cart_comm, MPI_STATUS_IGNORE);
		
	MPI_Probe(0, 2, cart_comm, &status_ja);
	MPI_Get_count(&status_ja, MPI_INT, &size_ja);
	MPI_Recv(local_ja, size_ja, MPI_INT, 0, 2,cart_comm, MPI_STATUS_IGNORE);
}
}
