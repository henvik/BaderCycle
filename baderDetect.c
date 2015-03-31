#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include<string.h>

#include"globalVars.h"
#include"baderDetect.h"
#include"GraphUtilities.h"

#define WHITE 'w'
#define RED 'r'
#define BLACK 'b'
#define GREEN 'g'



char * color; //Array for colorcoding

Edge* trans_arcs; //Trans arcs
int trans_arcs_size;
int trans_arcs_count=0;

AdjLst* adjecent;  //Holds the R_v list



// BAderDetect - detection of cycles in a distributed directed graph. 

void discovery(int num_vert){
	//printf("Number of vertices is: %d\n",num_vert);

	color=(char*) malloc(num_vert*sizeof(char)); 
	adjecent=(AdjLst *)malloc(sizeof(AdjLst)*num_vert); //adjecent is a list of pointers. Each pointer points to an array where the R_v's are stored
 	
 	trans_arcs_size=sqrt(num_vert)*4; //same as the size of the "border". Changed later if needed.
 	trans_arcs=(Edge *)malloc(sizeof(Edge)*trans_arcs_size);
	
	for(int i=0; i<num_vert;i++){
		color[i]=WHITE;
	} 
	for(int i=num_vert-1; i>=0;i--){
		if(color[i]==WHITE){
			adjecent[i]=visit(i);
		}
	}

	free(color);
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
	//		printf("External edge: %d, on proc %d\n", local_ja[w],procNr);
			add_trans_arc( v,  w, trans_arcs,trans_arcs_count);
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

void concatenate(int* A, int* B,int* lengthA, int* lengthB){
	A=realloc(A,sizeof(int)*(*lengthA+*lengthB));

	memcpy(A+*lengthA,B,*lengthB*sizeof(int));	
	*lengthA=*lengthA+*lengthB;
}
Edge new_edge(int v, int w){
	Edge new={v,w};
	return new;
}

AdjLst new_AdjLst(){
	AdjLst new;
	new.size=0;
	new.list=NULL;
	return new;
}

AdjLst merge_AdjLst(AdjLst A, AdjLst B){
	if(A.size==0 &&  B.size==0){
		return A;
	}else{
	A.list=(int *) realloc(A.list,sizeof(int)*(A.size+B.size));
	memcpy(A.list+A.size,B.list,B.size);
	A.size=A.size+B.size;
	return A;
	}

}

AdjLst add_Edge(AdjLst A, int edge){
	if(A.size==0){
		A.list=(int *)malloc(sizeof(int));
		A.size++;
	}else{
		A.size++;
		A.list=(int *) realloc(A.list,sizeof(int)*(A.size));
	
	}
	A.list[A.size-1]=edge;
	return A;
}

void add_trans_arc(int v, int w,Edge* trans_arcs, int trans_arcs_count){
	if(trans_arcs_count==trans_arcs_size){
		trans_arcs=realloc(trans_arcs, (trans_arcs_size*2));
		trans_arcs_size=trans_arcs_size*2;
	}

	trans_arcs[trans_arcs_count++]=new_edge(v,w);
	adjecent[v]=add_Edge(adjecent[v],w);
}
