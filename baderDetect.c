#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include<string.h>

#include"globalVars.h"
#include"baderDetect.h"

#define WHITE 'w'
#define RED 'r'
#define BLACK 'b'
#define GREEN 'g'

char * color; //Array for colorcoding
int* trans_arcs; //Trans arcs
int** adjecent;  //Holds the R_v list
int* adj_size; //Holds the sizes of each R_v list 

// BAderDetect - detection of cycles in a distributed directed graph. 

void discovery(int num_vert){

	color=(char*) malloc(num_vert*sizeof(char)); 
	adjecent=(int**)malloc(sizeof(int*)*num_vert); //adjecent is a list of pointers. Each pointer points to an array where the R_v's are stored
	adj_size=(int*)malloc(sizeof(int)*num_vert); //adjecent is a list of pointers. 
 
	
	for(int i=0; i<num_vert;i++){
		color[i]=WHITE;
	} 
	for(int i=0; i<num_vert;i++){
		if(color[i]==WHITE){
		//	visit(i);
		}
	}


}

int* visit(int v){
	adj_size[v]=0;
	color[v]=RED;
	for(int w=local_ia[v]; w<local_ia[v+1]; w++){
		//if(){
			switch(color[local_ja[w]]){
				case WHITE:
					adjecent[w]=visit(w);
					concatenate(adjecent[v],adjecent[w], &adj_size[v], &adj_size[w]);
				case BLACK:
					concatenate(adjecent[v], adjecent[w], &adj_size[v], &adj_size[w]);	
				case RED:
					exit(0);
		
			}
		//}
	}

	return adjecent[v];
}

void concatenate(int* A, int* B,int* lengthA, int* lengthB){
	realloc(A,sizeof(int)*(*lengthA+*lengthB));

	memcpy(A+*lengthA,B,*lengthB*sizeof(int));	
	*lengthA=*lengthA+*lengthB;
}