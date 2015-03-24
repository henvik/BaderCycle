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

int** adjecent;  //Holds the R_v list
int* adj_size; //Holds the sizes of each R_v list 



// BAderDetect - detection of cycles in a distributed directed graph. 

void discovery(int num_vert){

	color=(char*) malloc(num_vert*sizeof(char)); 
	adjecent=(int**)malloc(sizeof(int*)*num_vert); //adjecent is a list of pointers. Each pointer points to an array where the R_v's are stored
	adj_size=(int*)malloc(sizeof(int)*num_vert); //adjecent is a list of pointers. 
 	
 	trans_arcs_size=sqrt(num_vert)*4; //same as the size of the "border". Changed later if needed.
 	trans_arcs=(Edge *)malloc(sizeof(Edge)*trans_arcs_size);
	
	for(int i=0; i<num_vert;i++){
		color[i]=WHITE;
	} 
	for(int i=num_vert; i>0;i--){
		if(color[i]==WHITE){
			visit(i);
		}
	}


}

int* visit(int v){
	adjecent[v]=NULL;
	adj_size[v]=0;
	color[v]=RED;
	for(int w=local_ia[v]; w<local_ia[v+1]; w++){
		printf("Visiting neighbors of vertex %d, neighbor %d \n",local_map[v],w);

		int procNr=procNrFromCell(local_dims, gridDims,  local_map[w]);
		printf("Rank is: %d. Proc nr for %d is %d \n",rank,local_map[w],procNr);
		if(procNr==rank){
			printf("Internal edge: %d, on proc %d\n", w,rank);
			switch(color[local_ja[w]]){
				case WHITE:
					adjecent[w]=visit(w);
					concatenate(adjecent[v],adjecent[w], &adj_size[v], &adj_size[w]);
				case BLACK:
					concatenate(adjecent[v], adjecent[w], &adj_size[v], &adj_size[w]);	
				case RED:
					printf("Cycle found, im outta here!\n");
					MPI_Finalize();
					exit(0);
		
			}
		}
		else{
			printf("External edge: %d, on proc %d\n", local_map[w],procNr);
			add_trans_arc( v,  w, trans_arcs,trans_arcs_count);
		}

	}
	if(&adj_size[v]==0){
		color[v]=GREEN;
	}
	else{
		color[v]=BLACK;
	}

	return adjecent[v];
}

void concatenate(int* A, int* B,int* lengthA, int* lengthB){
	realloc(A,sizeof(int)*(*lengthA+*lengthB));

	memcpy(A+*lengthA,B,*lengthB*sizeof(int));	
	*lengthA=*lengthA+*lengthB;
}
Edge new_edge(int v, int w){
	Edge new={v,w};
	return new;
}

void add_trans_arc(int v, int w,Edge* trans_arcs, int trans_arcs_count){
	if(trans_arcs_count==trans_arcs_size){
		realloc(trans_arcs, (trans_arcs_size*2));
		trans_arcs_size=trans_arcs_size*2;
	}
	trans_arcs[trans_arcs_count++]=new_edge(v,w);

	realloc(adjecent[v],sizeof(int)*(adj_size[v]+1));
	printf("adj_size[v]: %d",adj_size[v]);
	//adjecent[v][adj_size[v]]=w;
	//memcpy(adjecent[v]+adj_size[v],&w,sizeof(int));	
	//adj_size[v]++;
}