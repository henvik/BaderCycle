#ifndef GRAPHUTILS_INCLUDED
#define GRAPHUTILS_INCLUDED

typedef struct {
   int    v;
   int    w;
}Edge;

typedef struct {
   int size;
   int* list;
}AdjLst;

typedef struct {
   int size;
   Edge *list;
}EdgeLst;



typedef struct ExVert{
	int vert_num;
	int procNr;
	
	struct TransArc *trans_arcs;
	struct ExArc *exp_arcs;
	
	struct ExVert *next;
}ExVert;



typedef struct TransArc{
	ExVert *head;
	struct TransArc *next;
}TransArc;

typedef struct ExArc{
	ExVert *head;
	struct ExArc *next;
}ExArc;

typedef struct ExpGraph{
	int num_vert;
	ExVert* first;
}ExpGraph;


/*Reads the grid from a .csv file and saves it as an adjecancy list in ia and ja.
	Input: 
		- file : char string containing the path
				 to the csv file
	Output: ia: pointer to a string of integers.
				 The i'th element of the lists contains the 
				index in ja of the first edge for node i,
			ja: pointer to a string of integer. 
				The ia[i]'th element of ja conatins the first
				edge of the i'th vertex.
*/
void importGrid(char* file, int **ia, int **ja, int *nv);

/*Prints an adjecancy graph 
	Input: 
		ia: pointer to a string of integers. The i'th element of the lists contains the 
			index in ja of the first edge for node i,
		ja: pointer to a string of integer. The ia[i]'th element of ja conatins the first
			edge of the i'th vertex.
*/

void printGraph(int* ia, int *ja, int nv);

/*Prints an adjecancy graph, where the indices in ia have been maped to global coordinates contained in map
	Input:
		ia: pointer to a string of integers. The i'th element of the lists contains the 
			index in ja of the first edge for node map[i],
		ja: pointer to a string of integer. The ia[i]'th element of ja conatins the first
			edge of the map[i]'th vertex.
 */
void printMappedGraph(int* ia, int *ja,int *map, int nv);

/*Calculates the rank of the processor owning a cell with a certain cell number
	Input:
		LocalGridDims: array specifying the dimensions of the local grids, distributed 
					   across all of the processors in the network
		CartGridDims:  array specifying the dimensions of the global grid.
		CellNr:		   the global cell number of the cell
	Output: 
		Integer specifying the rank of the processor which owns the cell.
*/
int procNrFromCell(int LocalGridDims[], int CartGridDims[], int CellNr);


/* Takes the global coordinates of a cell and returns its local coordinates  
	Input: 
		globalCoords: the coordinates of the cell in reference to the global grid
		localDims: the dimension of the local grid
	Output: 
		localCoords: the coordinates of the cell with respect to the local grid
*/
void globalCoords2localCoords(int globalCoords[],int localDims[],  int* localCoords );

/*Converts local coordinates to local cell nr
	Input:
		loocalCoords: the local coordinates
		localDims: dimension to the local cartesian grid
	Output:
		local cellnr  
 */
int  localCoords2localCellNr(int localCoords[], int localDims[]);

/*Converts global cell nr to global cartesian coordinates
	Input:
		cellNr: the global cell number
		GlobalDims: array containing the dimensions of the global cartesian grid
	Output:
		output: pointer to an array containing the cartesian coordinates in the global grid
*/
void cellNr2cartCoord(int cellNr, int* GlobalDims, int* output);

/*Converts cartesian coordinates to global cell nr
	Input:
		x: x-coordinate in the global grid
		y: y-coordinate in the global grid
		dims: dimensions of the global grid
	Output:
		the cellNr 
*/
int cartCoord2cellNr( int x, int y, int* dims);

/*Converts global cell nr to local cell nr
	Input: 
		globalCellNr: the cell nr of the cell in the global grid
		globalDims: the dimensions of the global grid
		localDims: the dimension of the local grids
	output:
		the cell nr in the local grid
*/
int globalCellNr2localCellNr(int globalCellNr, int globalDims[],int localDims[]);

/*Concatenates the two arrays
	Input: 
		A: array of integer
		B: array of integers
		lengthA: number of elements in A
		lengthB: number of elements in B
	output: 
		A: contains the concatenated array containing the elements of A followed by the elements of A
		lengthA: pointer to integer specifying the number of elements in the new array
*/
void concatenate(int* A, int* B, int* lengthA, int* lengthB);


Edge new_edge(int v, int w);
void add_trans_arc(int v, int w,EdgeLst *trans_arcs, int *trans_arcs_count,AdjLst* adjecent);

/*Creates a new instance of the structure AdjLst 
	Input: 
		NONE
	Output:
		new AdjLst allocated on the heap
*/
AdjLst new_AdjLst();

/*Merges two AdjLsts into one.
	Input: 
		A: AdjLst 
		B: AdjLst
	Output:
		A: New AdjLst now containing the elments of A and B.

*/
void merge_AdjLst(AdjLst *A, AdjLst *B);

/*Adds an edge to an AdjLSt
	Input: 
		A: the AdjLst that we want to add an edge to
		edge: the cellNr to the edge that we want to add
	Output: 
		A: the AdjLst now containing the new edge
*/
void add_Edge(AdjLst *A, int edge);

/* Searches an AdjLst for a certain edge
	Input: 
		A: the AdjLst that we want to search
		x: the cellNr of the edge that we are searching for
	Output:
		boolean value. True if x is in A. False otherwise
*/
bool isInAdjLst(AdjLst A,int x);


/*No longer in use*/
EdgeLst new_EdgeLst(int size);
void add_EdgeToLst(Edge e,int index, EdgeLst *lst);
void add_edges(EdgeLst *ExpGraph, EdgeLst *recv_buffer);
void expand_buffer(EdgeLst *recv_buffer, int recv_size);


//------------------------ExpGraph
/*Creates a new instance of the structure EXpGraph
	Input:
		takes no input
	Output:
		new pointer to a new ExpGraph element
*/
ExpGraph* new_ExpGraph();

/*Creates a new instance of the structure ExVert
	Input:
		vert_num: unique vertex number identifying the vertex
		procNr: rank of the processor on which the vertex currently resides.
	Output:
		new instance of an ExVert
*/
ExVert* new_ExVert(int vert_num, int procNr);

/*Adds a vertex two an express graph
	Input: 
		G: express graph 
		vert: vertex 
	Output:
		G: express graph now containing the new vertex
*/
void addExVert(ExpGraph *G, ExVert *vert);

/*Adds a trans arc to an express vertex
	Input: 
		vert: express vertex 
		trans_arc: the trans arc 
	Output:
		vert: express vertex now containing new trans arc
 */
void addTransArc(ExVert *vert, TransArc *trans_arc);

/*Adds Express arc to express vertex
	Input: 
		vert: Express vertex
		ex_arc: express arc
	Output: 
		vert: express vertex now containing new express arc
 */
void addExArc(ExVert *vert, ExArc *ex_arc);

/*Creates new instance of the structure TransArc
	Input:
		proc_nr: rank of the processor owning the terminal vertex of the trans arc
		vert_num: the vertex number of the terminal vertex. 
	Output:
		pointer to a new TransArc
*/
TransArc* newTransArc(ExVert *head);


/*Creates a new instance of the structure ExArc
	Input: 
		vert_num: vertex number of the terminal vertex.
		proc_nr: rank of the processor owning the terminal vertex of the express arc.
	Output:
		pointer to a new ExArc
*/
ExArc* newExArc(ExVert *head);

/*Searches through an ExVerts express arcs looking for certain express arc
	Input:
		vert: the vertex to be searched
		ex_arc: express arc
	Output:
		Returns true if ex_arc is found in vert. False otherwise.
*/
bool isInExp_arcs(ExVert *vert,ExArc *ex_arc);

/*Searches an express graph for a vertex
	Input: 
		G: express graph
		v: vertex number of the sought vertex
	Output:
		Pointer to vertex. NULL if vertex is not found.
*/
ExVert* FindExVert(ExpGraph *G, int v);

/*Prints the express graph in a formatted matter
	Input: 
		G: express graph
*/
void printfExpGraph(ExpGraph *G);


/*Adds a trans arc to the express graph G, and a new edge to the adjacency list of exit the vertex	
	Input:
		G: Express graph
		internal_vert_num: vertex number of the exit vertex
		external_vert_num: vertex number of the entrance vertex
		in_procNr: the rank of the processor owning the exit vertex
		external_procNr: rank of the processor owning the entrance vertex
		adjacent: adjacency list of the exit vertex.
	Output:
		G: express graph with the new trans arc added
		adjacent: adjacency list with the new edge added
*/
void addExternalEdge(ExpGraph *G, int internal_vert_num, int external_vert_num, int in_procNr ,int external_procNr, AdjLst *adjecent);

/*Checks if a certain vertex is reachable from another
	Input:
		adjecent: Adjacency list of the initial vertex
		vert_num: vertex number of vertex in question
	Output:
		True if vert_num is reachable. False otherwise.
*/
bool isReachable(AdjLst adjecent,int vert_num);

/*Merges the vertices of exp1 and exp2 into one ExpGraph
	Input:
		exp1: express graph
		exp2: express graph
	Output:
		new express graph which vertices are the union of the vertices in exp1 and exp2
*/
ExpGraph* mergeVertices(ExpGraph *exp1,ExpGraph *exp2, int origin2);

/*Removes vertex from express graph
	Input:
		G: express graph
		vert: vertex to be removed
	Output:
		G: express graph with vert removed
*/

void removeExVert(ExpGraph *G,ExVert *vert);

/*Removes express arc from express vertex
	Input: 
		vert: express vertex 
		exarc: express arc to be removed
	Output:
		vert: express vertex with exarc removed
*/
void removeExArc(ExVert *vert, ExArc *exarc);

#endif
