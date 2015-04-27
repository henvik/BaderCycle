#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include<string.h>

#include"globalVars.h"
#include"baderDetect.h"
#include"GraphUtilities.h"
#include"mpiUtils.h"

#define WHITE 'w'
#define RED 'r'
#define BLACK 'b'
#define GREEN 'g'
#define NEIGHBOURS 4  //north, south, east, west


char * color; //Array for colorcoding

EdgeLst trans_arcs; //Trans arcs
int trans_arcs_size;
int trans_arcs_count;

AdjLst* adjecent;  //Holds the R_v list



// BAderDetect - detection of cycles in a distributed directed graph. 

void discovery(int num_vert){
	//printf("Number of vertices is: %d\n",num_vert);

	color=(char*) malloc(num_vert*sizeof(char)); 
	adjecent=(AdjLst *)malloc(sizeof(AdjLst)*num_vert); //adjecent is a list of pointers. Each pointer points to an array where the R_v's are stored
 	
 	trans_arcs_size=sqrt(num_vert)*4; //same as the size of the "border". Changed later if needed.
 	trans_arcs=new_EdgeLst(trans_arcs_size);
	trans_arcs_count=0;
	
	for(int i=0; i<num_vert;i++){
		color[i]=WHITE;
	} 
	for(int i=num_vert-1; i>=0;i--){
		if(color[i]==WHITE){
			adjecent[i]=visit(i);
		}
	}
// 	printf("Trans_arcs_count: %d\n",trans_arcs_count);
	
	EdgeLst expGraph=new_EdgeLst(0);;
	
	comm_transArcs(&expGraph);
	merge(&expGraph);
	free(color);
	for(int i=0;i<num_vert;i++){
		if(adjecent[i].size!=0){
			free(adjecent[i].list);
		}
	}
	free(adjecent);
	free(trans_arcs.list);
}

AdjLst visit(int v){
//	printf("Visit called on local vertex %d, globCellNr %d.\n",v,local_map[v]);
	adjecent[v]=new_AdjLst();
	color[v]=RED;
	for(int w=local_ia[v]; w<local_ia[v+1]; w++){
		int procNr=procNrFromCell(local_dims, gridDims,local_ja[w]);
		int locCellNr=globalCellNr2localCellNr(local_ja[w],gridDims,local_dims);
		if(procNr==rank){
	//		printf("Internal edge: %d, on proc %d. Local cellNr is %d\n", local_ja[w],rank,globalCellNr2localCellNr(local_ja[w],gridDims,local_dims));
			switch(color[locCellNr]){
				case WHITE:
					adjecent[locCellNr]=visit(locCellNr);
					adjecent[v]=merge_AdjLst(adjecent[v],adjecent[locCellNr]);
					break;
				case BLACK:
					adjecent[v]=merge_AdjLst(adjecent[v],adjecent[locCellNr]);
					break;
				case RED:
					printf("Cycle found %d, im outta here!\n",local_map[v]);
					MPI_Finalize();
					exit(0);
		
			}
		}
		else{
		//	printf("External edge: %d, on proc %d, to %d\n", local_ja[w],rank,procNr);
			add_trans_arc( v,  local_ja[w], &trans_arcs, &trans_arcs_count, adjecent);
		}

	}
	if(adjecent[v].size==0){
		color[v]=GREEN;
	}
	else{
		color[v]=BLACK;
	}
	return adjecent[v];
}

void merge(EdgeLst* expGraph){
	EdgeLst* graphBuff;
	for(int h=0;h<log2(size);h++){
		if(last(rank,h)==0){
			if(test(rank,h)==0){
				//RecieveExp(procNr,graphBuff);
				printf("Rank=%d. Recive expG from processor %d\n",rank, set(rank,h));
				//MergeGraphs(expGraph, graphBuff);			
			}else{
				printf("Rank=%d. Send expG to processor %d \n",rank,clear(rank,h));
				//SendExp(procNr, graphBuff);
			}
		}	
	}

}

int last(int z, int h){
	int res=0;
	for(int i=0;i<h;i++){
		if((z>>i)&1){
		res |=1<<i;
		}
	}
	return res;
}

int test(int z, int h){
	return ((z>>(h))&1);
}

int set(int z, int h){
	return z|=1<<(h);
}

int clear(int z, int h){
	z &= ~(1 << (h));
	return z;
}

void RecieveExp(int procNr, EdgeLst* graphBuff){

}

void SendExp(int procNr, EdgeLst* graphBuff){

}

void MergeGraphsEx(EdgeLst* expGraph, EdgeLst* graphBuff){
	
}

void comm_transArcs(EdgeLst *expGraph){
	//send/recv buffers
	EdgeLst* trans_buffer=(EdgeLst *)malloc(NEIGHBOURS*sizeof(EdgeLst));
	// *expGraph=new_EdgeLst(0);
	printf("OK SO FAR\n");

	int *buffSizes=(int *)calloc(NEIGHBOURS, sizeof(int));
	for(int i=0; i<NEIGHBOURS;i++){
		//pre allocating memory for send/recv buffers. reallocating later if needed.
		trans_buffer[i]=new_EdgeLst(trans_arcs_count*2/NEIGHBOURS);
	}
 	build_transBuffer(buffSizes,trans_buffer);
/*	for(int i=0;i<4;i++)
	printf("Buffsizes for %d is %d\n",i,buffSizes[i]);	
	 	if(west!=-2){
 		for(int i=0;i<buffSizes[2];i++){
 			printf("trans_buffer for %d are %d, %d\n", rank,trans_buffer[2].list[i].v,trans_buffer[2].list[i].w);
 		}
 	}*/
	MPI_Datatype MPI_Edge;
	def_datatypes(&MPI_Edge);

	MPI_Request send_req[NEIGHBOURS],recv_req[NEIGHBOURS];


	if(north!=-2){	
		MPI_Isend(trans_buffer[0].list,buffSizes[0],MPI_Edge,north,0,cart_comm,&send_req[0]);
	}
	if(south!=-2){	
		MPI_Isend(trans_buffer[1].list,buffSizes[1],MPI_Edge,south,1,cart_comm,&send_req[1]);
	}
	if(west!=-2){	
		MPI_Isend(trans_buffer[2].list,buffSizes[2],MPI_Edge,west,2,cart_comm,&send_req[2]);
	}
	if(east!=-2){	
		MPI_Isend(trans_buffer[3].list,buffSizes[3],MPI_Edge,east,3,cart_comm,&send_req[3]);
	}
	int added=0;
	
	if(north!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(north,1,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		expand_buffer(expGraph,recv_size);
		MPI_Recv(expGraph->list+added,recv_size,MPI_Edge,north,1,cart_comm,MPI_STATUS_IGNORE);
		added+=recv_size;

	}
	
	
	if(south!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(south,0,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
	
		expand_buffer(expGraph,recv_size);
		MPI_Recv(expGraph->list+added,recv_size,MPI_Edge,south,0,cart_comm,MPI_STATUS_IGNORE);
		added+=recv_size;
	}
	
	
	if(west!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(west,3,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		expand_buffer(expGraph,recv_size);
		MPI_Recv(expGraph->list+added,recv_size,MPI_Edge,west,3,cart_comm,MPI_STATUS_IGNORE);	
		added+=recv_size;
	}
	
	
	if(east!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(east,2,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		expand_buffer(expGraph,recv_size);
		MPI_Recv(expGraph->list+added,recv_size,MPI_Edge,east,2,cart_comm,MPI_STATUS_IGNORE);
		added+=recv_size;
	}
	
	
	if(north!=-2){	
		MPI_Wait(&send_req[0],MPI_STATUS_IGNORE);
	}
	if(south!=-2){	
		MPI_Wait(&send_req[1],MPI_STATUS_IGNORE);
	}
	if(west!=-2){	
		MPI_Wait(&send_req[2],MPI_STATUS_IGNORE);
	}
	if(east!=-2){	
		MPI_Wait(&send_req[3],MPI_STATUS_IGNORE);
	}
// 	printf("ExpGraph size=%d\n",expGraph.size);
// 	for(int i =0;i<expGraph.size;i++){
// 		printf("expGraph on %d nr %d: %d ->%d\n",rank,i,expGraph.list[i].v,expGraph.list[i].w);
// 	}

	free(buffSizes);
	for(int i=0; i<NEIGHBOURS;i++){
		free(trans_buffer[i].list);
	}
	
}


void build_transBuffer(int * buffSizes, EdgeLst* trans_buffer){
	int procNr;
	for(int i=0; i<trans_arcs_count;i++){
		procNr=procNrFromCell(local_dims,gridDims,trans_arcs.list[i].w);
// 		printf("On processor %d we have %d trans_arcs. This is the %d'th. It goes to %d on %d\n", rank, trans_arcs_count,i,trans_arcs[i].w,procNr);
		//makes a send buffer for each of the (possibly) four neighbours. 
		//The four directions maps to the following indices in the transferbuffer:
		//north=0, south=1, west=2 and east=3,.
		if(procNr==north){
			add_EdgeToLst(trans_arcs.list[i],buffSizes[0],&trans_buffer[0]);
			buffSizes[0]+=1;
		}
		if(procNr==south){
			add_EdgeToLst(trans_arcs.list[i],buffSizes[1],&trans_buffer[1]);
			buffSizes[1]+=1;
		}
		if(procNr==west){
			add_EdgeToLst(trans_arcs.list[i],buffSizes[2],&trans_buffer[2]);
			buffSizes[2]+=1;
		}
		if(procNr==east){
			add_EdgeToLst(trans_arcs.list[i],buffSizes[3],&trans_buffer[3]);
			buffSizes[3]+=1;
		}
	}
}


