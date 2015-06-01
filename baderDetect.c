#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include<string.h>
#include<stdbool.h>

#include"globalVars.h"
#include"baderDetect.h"
#include"GraphUtilities.h"
#include"mpiUtils.h"

#define WHITE 'w'
#define RED 'r'
#define BLACK 'b'
#define GREEN 'g'
#define NEIGHBOURS 4  //north, south, east, west
#define INITSIZE 5

char * color; //Array for colorcoding

EdgeLst trans_arcs; //Trans arcs
int trans_arcs_size;
int trans_arcs_count;

AdjLst* adjecent;  //Holds the R_v list

ExpGraph *expGraph;


// BAderDetect - detection of cycles in a distributed directed graph. 

ExpGraph* discovery(int num_vert){
	//printf("Number of vertices is: %d\n",num_vert);

	color=(char*) malloc(num_vert*sizeof(char)); 
	adjecent=(AdjLst *)malloc(sizeof(AdjLst)*num_vert); //adjecent is a list of pointers. Each pointer points to an array where the R_v's are stored
 	
 	trans_arcs_size=sqrt(num_vert)*4; //same as the size of the "border". Changed later if needed.
 	trans_arcs=new_EdgeLst(trans_arcs_size);
	trans_arcs_count=0;
	
	expGraph=new_ExpGraph();
	
	for(int i=0; i<num_vert;i++){
		color[i]=WHITE;
	} 
	for(int i=num_vert-1; i>=0;i--){
		if(color[i]==WHITE){
			adjecent[i]=visit(i);
		}
	}
// 	printf("Trans_arcs_count: %d\n",trans_arcs_count);
	if(rank==2){	
		for(int i=0;i<num_vert;i++){
		
			printf(" %d:",local_map[i]);
			for(int j=0;j<adjecent[i].size;j++){
				printf(" %d,",adjecent[i].list[j]);
			}
			printf("\n");
		}
	}
	//EdgeLst expGraph=new_EdgeLst(0);;
	
	comm_transArcs(expGraph);
	completeExpGraph(expGraph,adjecent);
	if(rank==2)
	printfExpGraph(expGraph);

//	merge(&expGraph);
	free(color);

	for(int i=0;i<num_vert;i++){
		if(adjecent[i].size!=0){
			free(adjecent[i].list);
		}
	}
	free(adjecent);
	free(trans_arcs.list);
	return expGraph;
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
					merge_AdjLst(&adjecent[v],&adjecent[locCellNr]);
					break;
				case BLACK:
					merge_AdjLst(&adjecent[v],&adjecent[locCellNr]);
					break;
				case RED:
					printf("Cycle found %d, im outta here!\n",local_map[v]);
					MPI_Finalize();
					exit(0);
			}
		}
		else{
		//	printf("External edge: %d, on proc %d, to %d\n", local_ja[w],rank,procNr);
			//add_trans_arc( v,  local_ja[w], &trans_arcs, &trans_arcs_count, adjecent);
			addExternalEdge(expGraph, local_map[v], local_ja[w],rank,procNr, &adjecent[v]);
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


void completeExpGraph(ExpGraph *G,AdjLst *adjecent){
	//Add express acrs
	ExVert *current=G->first;
	while(current){
		if(current->procNr!=rank){
			TransArc *cTarc=current->trans_arcs;
			while(cTarc){
				int index=globalCellNr2localCellNr(cTarc->vert_num,gridDims,local_dims);
				for(int i=0; i<adjecent[index].size;i++){
					addExArc(current,newExArc(adjecent[index].list[i], cTarc->proc_nr));
				}
				cTarc=cTarc->next;
			}
			current=current->next;
			continue;
		}
		ExVert *iterater=G->first;
		while(iterater){
			if(iterater==current){
				iterater=iterater->next;
				continue;
			}
			TransArc *cTarc=iterater->trans_arcs;	
			while(cTarc){
				int index=globalCellNr2localCellNr(current->vert_num,gridDims,local_dims);
				if(isReachable(adjecent[index],cTarc->vert_num)){
					addExArc(current,newExArc(cTarc->vert_num, cTarc->proc_nr));
				}
				cTarc=cTarc->next;
			}
			iterater=iterater->next;
	}
	current=current->next;	
	}
	//Pack for communication
	
}




void merge(ExpGraph* expGraph){
	EdgeLst* graphBuff;
	for(int h=0;h<log2(size);h++){
		if(last(rank,h)==0){
			if(test(rank,h)==0){
		//		printf("Rank=%d. Recive expG from processor %d\n",rank, set(rank,h));
				
				int procNr=set(rank,h);
				ExpGraph *expRecieved;
				RecieveExp(procNr,&expRecieved);
								
				MergeGraphs(expGraph, expRecieved,rank,procNr);	
						
			}else{
		//		printf("Rank=%d. Send expG to processor %d \n",rank,clear(rank,h));
				 int procNr=clear(rank,h);
				SendExp(procNr, expGraph);
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

void SendExp(int procNr, ExpGraph* exp){
		 
		int *send_buffer=(int *)malloc(sizeof(int)*exp->num_vert*12); //magic number 12 is semi-arbitrary to allocate enough memory for each vertex.
		int send_count=packGraph(exp, &send_buffer);
		
	//	printf("we are sending %d ints to  procNr %d\n",send_count,procNr);
		//printBuffer(&send_buffer,send_count);
		MPI_Send(send_buffer,send_count,MPI_INT,procNr,0,cart_comm);
		free(send_buffer);
}





//Packs the Express graph in the following manner:
//All the vertices are stored sequentially with a -1 signaling the end of one vertice and the start of the next.
//Within each vertice the info is stored as following. The first element contains the vertice number, the second contians the procNr.
//The trans arc are stored in pairs of two integers. The first is the vertice, the second the proc nr. 
//A -2 signals the end of transArcs and the start of express arcs. These are stored in pairs of two, like the trans arcs.
//A -9 signals the end of the express graph.
int packGraph(ExpGraph *exp, int **send_buffer){
	int size_send_buffer=exp->num_vert*12;
	
	int counter=0;
	ExVert *current=exp->first;
	while(current){
		addToBuffer(send_buffer,current->vert_num,&counter,&size_send_buffer);
		addToBuffer(send_buffer,current->procNr,&counter,&size_send_buffer);

		TransArc *cTarc=current->trans_arcs;
		while(cTarc){
			
			addToBuffer(send_buffer,cTarc->vert_num,&counter,&size_send_buffer);
			addToBuffer(send_buffer,cTarc->proc_nr,&counter,&size_send_buffer);
			cTarc=cTarc->next;
		}
		addToBuffer(send_buffer,-2,&counter,&size_send_buffer);
		
		ExArc *cEarc=current->exp_arcs;
		while(cEarc){	
			addToBuffer(send_buffer,cEarc->vert_num,&counter,&size_send_buffer);
			addToBuffer(send_buffer,cEarc->proc_nr,&counter,&size_send_buffer);
			cEarc=cEarc->next;
		}

		addToBuffer(send_buffer,-1,&counter,&size_send_buffer);
		current=current->next;
	}
	addToBuffer(send_buffer,-9,&counter,&size_send_buffer);
	return counter;
}

void addToBuffer(int **send_buffer,int toBeAdded ,int *counter, int *size_send_buffer){
	if((*counter)>=(*size_send_buffer)){
		*send_buffer = (int *) realloc(*send_buffer, sizeof(int)*(*size_send_buffer)*2);
		(*size_send_buffer)*=2;
	}
//	printf("Now counter are: %d \n",(*counter));
	(*send_buffer)[(*counter)++]=toBeAdded;
}

void printBuffer(int **send_buffer,int counter){
	int* current=*send_buffer;
	int testcount=0;
	printf("%d, %d, %d\n",*current,*(current+1),(*current+2));
	while(*current!=-9){
		printf("Vertex %d  ",*current);
		printf("on processor %d has transArcs: \n     ",*(current+1) );
		current+=2;
		testcount+=2;
		while(*current!=-2){
			printf("%d, ",*current);
			current++;
			testcount++;
			
		}
		printf("\n      ");
		current++;
		testcount++;
		printf("and express arcs ");
		while(*current!=-1){
			printf("%d, ",*current);
			current++;
			testcount++;
		}
		printf("\n");
		current++;
		testcount++;
	}
		printf("Total of %d elements\n",testcount);
}

void RecieveExp(int procNr, ExpGraph** exp){
	MPI_Status recv_status;
	int recv_size;
	MPI_Probe(procNr,0,cart_comm,&recv_status);
	MPI_Get_count(&recv_status,MPI_INT,&recv_size);
		
	int *recv_buff =malloc(recv_size*sizeof(int));
	MPI_Recv(recv_buff,recv_size,MPI_INT,procNr,0,cart_comm,MPI_STATUS_IGNORE);
	//printBuffer(&recv_buff,recv_size);
	*exp=unPackRecvBuffer(recv_buff);
}

ExpGraph *unPackRecvBuffer(int *recv_buffer){
	ExpGraph *newExpGraph=new_ExpGraph();
	//printf("Got here\n");	
	int index=0;
	while(recv_buffer[index]!=-9){
		ExVert *newExVert=new_ExVert(recv_buffer[index],recv_buffer[index+1]);
		index+=2;

		while(recv_buffer[index]!=-2){
			addTransArc(newExVert, newTransArc(recv_buffer[index+1], recv_buffer[index]));
			//NB: Note that we here reverse the input of vert_num and procNr because of the definition of newTransArc.
			index+=2;
		}
		index++;
		
		while(recv_buffer[index]!=-1){
			addExArc(newExVert,  newExArc(recv_buffer[index], recv_buffer[index+1] ) );
			index+=2;
	
		}
		addExVert(newExpGraph,newExVert);
		index++;
	}
	return newExpGraph;
}


void MergeGraphs(ExpGraph* exp1, ExpGraph* exp2, int origin1, int origin2){
	ExpGraph *exp0=mergeVertices(exp1,exp2);
	
	//Iterate through trans_arcs in exp1 and exp2 
	ExVert *cExVert=exp1->first;
	while(cExVert){
		if(cExVert->procNr==origin1){//check whether this is an entrance vertex
			cExVert=cExVert->next;
			continue;
		}	
		TransArc *cTarc=cExVert->trans_arcs;
		while(cTarc){
			ExVert *vert=FindExVert(exp0,cTarc->vert_num);
			if(!vert){
				addTransArc(FindExVert(exp0,cExVert->vert_num),newTransArc(cTarc->proc_nr,cTarc->vert_num));
			}else{
				if(isExArc(vert,cExVert)){
					printf("cycle found\n");
					MPI_Finalize();
					exit(0);
				}
				removeAndReplaceTransArc(exp0,cExVert,vert);

			}
			cTarc=cTarc->next;
		}		
		cExVert=cExVert->next;
	}
}

bool isExArc(ExVert *terminal ,ExVert *initial){
	ExArc *current=terminal->exp_arcs;
	while(current){
		if(current->vert_num==initial->vert_num){
			return true;
		}
		current=current->next;
	}
	return false;
}

void removeAndReplaceTransArc(ExpGraph *G,ExVert *init, ExVert *term ){
	ExVert *globalCurrent=G->first;
	while(globalCurrent){
	
		ExArc *localCurrent = globalCurrent->exp_arcs;
		while(localCurrent){
			if(localCurrent->vert_num==init->vert_num){
				transferExArcs(globalCurrent,init,term);
				removeExArc(globalCurrent,localCurrent);

			}
			localCurrent=localCurrent->next;
		}
	globalCurrent=globalCurrent->next;
	}
	removeExVert(G,init);
	removeExVert(G,term);
}

void transferExArcs(ExVert *globalCurrent, ExVert *init, ExVert *term){
	ExArc *currentArc =term->exp_arcs;
	while(currentArc){
		addExArc(globalCurrent,newExArc(currentArc->vert_num, rank));
		currentArc=currentArc->next;
	} 
}

void comm_transArcs(ExpGraph *expGraph){
	//send/recv buffers
	EdgeLst* trans_buffer=(EdgeLst *)malloc(NEIGHBOURS*sizeof(EdgeLst));
	//*expGraph=new_EdgeLst(0);
	printf("OK SO FAR\n");

	int *buffSizes=(int *)calloc(NEIGHBOURS, sizeof(int));
	for(int i=0; i<NEIGHBOURS;i++){
		//pre allocating memory for send/recv buffers. reallocating later if needed.
		trans_buffer[i]=new_EdgeLst(INITSIZE);
	}
 	build_transBuffer(buffSizes,trans_buffer, expGraph);
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
	
	Edge *recv_buff;
	if(north!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(north,1,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		recv_buff =malloc(recv_size*sizeof(Edge));
		
		MPI_Recv(recv_buff,recv_size,MPI_Edge,north,1,cart_comm,MPI_STATUS_IGNORE);
		recvTransArcs(recv_buff, recv_size, expGraph, north);	
		free(recv_buff);
		
	}
	
	
	if(south!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(south,0,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
	
		recv_buff =malloc(recv_size*sizeof(Edge));
		
		MPI_Recv(recv_buff,recv_size,MPI_Edge,south,0,cart_comm,MPI_STATUS_IGNORE);
		
		recvTransArcs(recv_buff, recv_size, expGraph, south);	
		free(recv_buff);
	}
	
	
	if(west!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(west,3,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		recv_buff =malloc(recv_size*sizeof(Edge));
		
		MPI_Recv(recv_buff,recv_size,MPI_Edge,west,3,cart_comm,MPI_STATUS_IGNORE);
		
		recvTransArcs(recv_buff, recv_size, expGraph, west);
		free(recv_buff);	
	}
	
	
	if(east!=-2){	
		MPI_Status recv_status;
		int recv_size;
		MPI_Probe(east,2,cart_comm,&recv_status);
		MPI_Get_count(&recv_status,MPI_Edge,&recv_size);
		
		recv_buff =malloc(recv_size*sizeof(Edge));
		
		MPI_Recv(recv_buff,recv_size,MPI_Edge,east,2,cart_comm,MPI_STATUS_IGNORE);
		
		recvTransArcs(recv_buff, recv_size, expGraph, east);
		free(recv_buff);		
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


void build_transBuffer(int * buffSizes, EdgeLst* trans_buffer, ExpGraph *expGraph){	
	ExVert* current=expGraph->first;
	
	//iterate through every trans arc of every  vertex.
	while(current){
		TransArc *currTransArc=current->trans_arcs;
		while(currTransArc){
// 			printf("On processor %d we have %d trans_arcs. This is the %d'th. It goes to %d on %d\n", rank, trans_arcs_count,i,trans_arcs[i].w,procNr);
			//makes a send buffer for each of the (possibly) four neighbours. 
			//The four directions maps to the following indices in the transferbuffer:
			//north=0, south=1, west=2 and east=3,.
			if(currTransArc->proc_nr==north){
				add_EdgeToLst(new_edge(current->vert_num,currTransArc->vert_num),buffSizes[0],&trans_buffer[0]);
				buffSizes[0]+=1;
			}
			if(currTransArc->proc_nr==south){
				add_EdgeToLst(new_edge(current->vert_num,currTransArc->vert_num),buffSizes[1],&trans_buffer[1]);
				buffSizes[1]+=1;
			}
			if(currTransArc->proc_nr==west){
				add_EdgeToLst(new_edge(current->vert_num,currTransArc->vert_num),buffSizes[2],&trans_buffer[2]);
				buffSizes[2]+=1;
			}
			if(currTransArc->proc_nr==east){
			
				add_EdgeToLst(new_edge(current->vert_num,currTransArc->vert_num),buffSizes[3],&trans_buffer[3]);
				buffSizes[3]+=1;
			}
		currTransArc=currTransArc->next;
		}
	current=current->next;
	}
	
}

void recvTransArcs(Edge* list, int size, ExpGraph *expGraph, int procNr){
	
	for(int i=0; i<size;i++){
		ExVert *new=new_ExVert(list[i].v,procNr);
		addExVert(expGraph,new);
		addTransArc(new, newTransArc( rank, list[i].w));
	}

}

