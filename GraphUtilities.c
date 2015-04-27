#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include<string.h>

#include"globalVars.h"
#include"DistGather.h"
#include"GraphUtilities.h"


void cellNr2cartCoord(int cellNr, int* GlobalDims, int* output){
	/* Input:       cellNr  - the cellNr in the grid 
	 *                         dims    - an array specifying the dimensions of the cartesian grid
	 *                                 Output: coords  - the cartesian coordinates of the cellNR
	 *                                 */
	        output[1]=cellNr/GlobalDims[1];
		        output[0]=cellNr%GlobalDims[0];
}




int cartCoord2cellNr( int x, int y, int* dims){
	/* Input:       x,y  - the cartesian coordinates of the cellNR
	 *                         dims    - an array specifying the dimensions of the cartesian grid
	 *                                 Output: cellNr  - the cellNr in the grid 
	 *                                 */
	        int output=y*dims[1]+x;
		        return output;
}


void globalCoords2localCoords(int globalCoords[],int localDims[],  int* localCoords ){
	/* Input: globalCoords  -  cartesian coordinates in the global grid,
	 *	  localDims 	-  dimension of the local subgrid
	 *Output: localCoords   -  cartesian coordinates in the subgrid	 
	 */
	localCoords[0]=globalCoords[0]%localDims[0];
	localCoords[1]=globalCoords[1]%localDims[1];
}

int  localCoords2localCellNr(int localCoords[],int  localDims[]){
	/*Input: localCoords - local cartesian coordinates
	 * 	 localDims   - dimension of the local subgrid
	 *Outpu: localCellNr - the cellNr in the local adjecency list
	 */

	int localCellNr=localCoords[1]*localDims[0]+localCoords[0];
	return localCellNr;
}
int globalCellNr2localCellNr(int globalCellNr, int globalDims[],int localDims[]){
	/*Input: globalCellNr  - cellNr in the global adjecency list
	 * Output: localCellNr - local cellNr in the local adjecency list
	 */
	int* tmp=(int *)malloc(2*sizeof(int));
	cellNr2cartCoord(globalCellNr,globalDims,tmp);
	globalCoords2localCoords(tmp,localDims,tmp);	
     	int localCellNr = localCoords2localCellNr(tmp,localDims);
 	free(tmp);
	return localCellNr;
}

int procNrFromCell(int LocalGridDims[], int CartGridDims[], int CellNr){ //Finds which processor a certain cell resides on
	int *CartCoords = (int*) malloc(2*sizeof(int));
	cellNr2cartCoord(CellNr, CartGridDims, CartCoords);
//	printf("Cartcoords are (%d,%d)\n",CartCoords[0],CartCoords[1]);
//	printf("LocalGridDims: (%d,%d) \n", LocalGridDims[0],LocalGridDims[1]);
	int procCoord[2];
	procCoord[0]=(CartCoords[0]) /(LocalGridDims[0]);
	procCoord[1]=(CartCoords[1]) / (LocalGridDims[1]);

	//printf("Cellnr: %d, CartCoords, (%d,%d), ProcCoord are (%d, %d)\n",CellNr, CartCoords[0], CartCoords[1], procCoord[0],procCoord[1]);
	int rank;
	MPI_Cart_rank(cart_comm, procCoord, &rank);
	return rank;
}

//Read a .csv file an make a directed graph in sparse matrix format
/*The file has to have the following format:
- The first line consists of a single number, indication the number of vertices in the graph.
- line number i indicates contains a list of vertices with edges from the i'th vertice on the format: k,l,m; 
*/
void importGrid(char* file, int **ia, int **ja, int *nv){
    
    FILE * fp;
    char line[100];
	int len=0;
	
	int l,i,tmp,charLeft;  
    int jaINDEX=0; //Index counter for the array ja.
    int numOfLines; //to store the number of lines

	char* buffer = (char*)malloc(sizeof(char)*32); //Char buffer
	fp=fopen(file,"r");	
	if(fp==NULL){
	   exit(EXIT_FAILURE);
	}
	int line_nr=0;
	//Reads the first line
	fgets(line,60,fp);
	numOfLines=atoi(line);
	*nv=numOfLines;
// Allocating memory
	int *Pia=(int *)malloc((numOfLines+1)*sizeof(int));
	int *Pja=(int *)malloc(numOfLines*4*sizeof(int));	
	//reads the file line by line, starting at the second line
	while(fgets(line,60,fp)){
	     i=0;
	     l=0;
	     charLeft=1;
	     int edgeCount=0;
	   //  As long as there is characters in the line
	     while(charLeft){
	     	//If we find the delimiter signaling the end of an integer
		   if(line[i]==','){
		   		//If this is the first edge from the current vertex
		   	  if(edgeCount==0){
		   	  	edgeCount++;
		   	    Pia[line_nr++]=jaINDEX; 
		   	  }
		   	  //We mark the end of the current integer and convert it to an int and add it to ja
		   	  buffer[l++]='\0';
			  l=0;
		      tmp=atoi(buffer);
		      Pja[jaINDEX++]=tmp;
		      i++;
		  //If we encounter the end of a line
	       }else if(line[i]==';'){
	       	  //If the line was empty we have to mark this
	       	  if(i==0){
	       	  	charLeft=0;
	       	  	Pia[line_nr++]=jaINDEX; 
	       	  	break;
	       	  }
	         charLeft=0;
	       }
	       //If non of the above we just read the next char in the line
	        else{
	   		  buffer[l++]=line[i++];
	   	  }
	   }
	}
	if(Pia[line_nr-1]==jaINDEX){
			Pia[line_nr]=jaINDEX;
	}
	else{
		Pia[line_nr]=jaINDEX;
	}
	//Set the pointers(couldn't do it any other way?)
	*ia=Pia;
	*ja=Pja;
	fclose(fp);
}


void printGraph(int* ia, int *ja, int nv){
	int i,j;
	for(i=0;i<nv;i++){
		printf("Vertex %d has edges to: ",i);
		for(j=ia[i];j<ia[i+1];j++){
			printf("%d, ",ja[j]);
		}
		printf("\n");
	}
}

void printMappedGraph(int* ia, int *ja,int *map, int nv){
	int i,j;
	for(i=0;i<nv;i++){
		printf("Vertex %d has edges to: ",map[i]);
		for(j=ia[i];j<ia[i+1];j++){
			printf("%d, ",ja[j]);
		}
		printf("\n");
	}
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

EdgeLst new_EdgeLst(int size){
	EdgeLst new;
	new.size=size;
	if(size==0){
		new.list=NULL;
	}else{
		new.list=(Edge *)malloc(size*sizeof(Edge));
	}
	return new;
}

void add_EdgeToLst(Edge e,int index, EdgeLst* lst){
	if(index>=lst->size){
		lst->size=lst->size+index*2;
		lst->list=realloc(lst->list,lst->size*sizeof(Edge));
	}
	lst->list[index]=e;
}


AdjLst merge_AdjLst(AdjLst A, AdjLst B){
	if(A.size==0 &&  B.size==0){
		return A;
	}else{
	A.list=(int *) realloc(A.list,sizeof(int)*(A.size+B.size));
	memcpy(A.list+A.size,B.list,B.size*sizeof(int));
	A.size=A.size+B.size;
	return A;
	}

}

void add_Edge(AdjLst *A, int edge){
	if(A->size==0){
		A->list=(int *)malloc(sizeof(int));
		A->size++;
	}else{
		A->size++;
		A->list=(int *) realloc(A->list,sizeof(int)*(A->size));
	
	}
	A->list[A->size-1]=edge;
	
}

void add_trans_arc(int v, int w,EdgeLst *trans_arcs, int *trans_arcs_count, AdjLst* adjecent){
	if(*trans_arcs_count==trans_arcs->size){
		trans_arcs->list=realloc(trans_arcs->list,sizeof(Edge)*(trans_arcs->size*2));
		trans_arcs->size=trans_arcs->size*2;
	}

	trans_arcs->list[*trans_arcs_count]=new_edge(local_map[v],w);
	*trans_arcs_count+=1;
	add_Edge(&adjecent[v],w);
}

void expand_buffer(EdgeLst *recv_buffer, int recv_size){
	if(recv_buffer->size==0){
		recv_buffer->list=(Edge *)malloc(sizeof(Edge)*recv_size);
		recv_buffer->size+=recv_size;
	}else{
		recv_buffer->list=realloc(recv_buffer->list,(recv_buffer->size + recv_size)*sizeof(Edge));
		recv_buffer->size+=recv_size;
	}
}
