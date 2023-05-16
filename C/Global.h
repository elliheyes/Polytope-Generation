/*  ======================================================================  */
/*  ==========	     			   	  	       	==========  */
/*  ==========                    G L O B A L                   ==========  */
/*  ==========						        ==========  */
/*  ======================================================================  */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>



/*  ============           	P A R A M E T E R S              ============  */

#define MIN -3		  /* lower bound on integers considered for points */
#define MAXNVRTS 14  	  /* maximum number of vertices */
#define POLYDIM 3	  /* dimension of polytopes */
#define BINLEN 3	  /* maximum length of binary number */
#define POPSIZE 350       /* population size */
#define NUMGEN 500        /* number of generations */
#define NUMCUTS 1         /* number of cuts in a crossing */
#define RANKING 0         /* macro for ranking method to select breeding pairs */
#define ROULETTE 1        /* macro for roulette method to select breading pairs */
#define MUTRATE 0.005     /* mutation rate */
#define ALPHA 3.0         /* value of alpha parameter in selection probability */
#define KEEPFITEST 1      /* copy fitest individual into next population */
#define DONTKEEPFITEST 0  /* don't copy over fitest individual */
#define MONITORON 1       /* terminal monitor on */
#define MONITOROFF 0      /* terminal monitor off */

/** fitness weights **/
#define DIST_WEIGHT 1     /* distance of the facets from the origin */
#define IP_WEIGHT 1       /* IP property */
#define NVERTS_WEIGHT 0   /* number of vertices */
#define NPTS_WEIGHT 0     /* number of points */
#define H11_WEIGHT 0      /* h11 hodge number */
#define H12_WEIGHT 0      /* h12 hodge number */
#define H13_WEIGHT 0      /* h13 hodge number */
#define H22_WEIGHT 0      /* h22 hodge number */
#define EULER_WEIGHT 0    /* euler number */
#define NVERTS 6    
#define NPTS 7
#define H11 1
#define H12 1
#define H13 1
#define H22 1
#define EULER 1

/* PALP limits */
#define POINT_Nmax 2000000         /* maximum number of points */
#define VERT_Nmax 128	           /* maximum number of vertices */
#define FACE_Nmax 10000	           /* maximum number of faces */
#define SYM_Nmax 46080             /* maximum number of symmetries */
#define EQUA_Nmax VERT_Nmax        /* maximum number of facets */		
#define AMBI_Dmax (5 * POLYDIM)    /* maximum dimension of the ambient space */
#define FIB_Nmax 3000              /* maximum number of allowed weight relations among points that define the IP simplices */

/* data types */
#define	Long long long    
#define LLong long long
#define GL_Long	Long 

#define INT_Nbits 32
#define LONG_LONG_Nbits 64
/*
These numbers should be set to the actual numbers of bits occupied by the
structures "unsigned int" and "unsigned long long" in your version of C.
If they are set to lower values, everything still works but may be
considerably slowed down.
*/

#if (VERT_Nmax <= INT_Nbits)
typedef		        unsigned int            INCI;
#elif (VERT_Nmax <= LONG_LONG_Nbits)
typedef		        unsigned long long	INCI;
#else
#define I_NUI     ((VERT_Nmax-1)/INT_Nbits+1)
typedef struct {unsigned int ui[I_NUI];}   INCI;
#endif
/*
An INCI encodes the incidence relations between a face and a list of
vertices as a bit pattern (1 if a vertex lies on the face, 0 otherwise).
Depending on the allowed number VERT_Nmax of vertices, a single "unsigned int"
or "unsigned long long" may be sufficient.
If VERT_Nmax is larger than the number of bits in a "long long integer", an
array of unsigned integers is used to simulate an integer type of the required
size.
*/

#if (VERT_Nmax <= LONG_LONG_Nbits)
#define INCI_M2(x)     ((x) % 2)              /* value of first bit      */
#define	INCI_AND(x,y)  ((x) & (y))            /* bitwise logical and     */
#define	INCI_OR(x,y)   ((x) | (y))            /* bitwise logical or      */
#define	INCI_XOR(x,y)  ((x) ^ (y))            /* bitwise exclusive or    */
#define	INCI_EQ(x,y)   ((x) == (y))           /* check on equality       */
#define INCI_LE(x,y)   INCI_EQ(INCI_OR(x,y),y)/* bitwise less or equal */
#define INCI_EQ_0(x)   INCI_EQ(x,INCI_0())    /* check if all bits = 0   */
#define INCI_0()       (0)                    /* set all bits to 0       */
#define INCI_1()       (1)                    /* set only first bit to 1 */
#define INCI_D2(x)     ((x) / 2)              /* shift by one bit        */
#define INCI_PN(x,y)   (2 * (x) + !(y))       /* shift and set first bit */
/*
For an INCI defined as a single unsigned (long long) integer whose bits are
regarded as representing incidences, these are useful definitions.
INCI_PN is particularly useful when a new vertex is added: if x represents
an equation E w.r.t. some vertex list and y is the result of evaluating E
on some new vertex V, then INCI_PN(x,y) represents x w.r.t. the vertex list
enhanced by V.
*/

#else
#define INCI_M2(x)      ((x).ui[0] % 2)
INCI INCI_AND(INCI x, INCI y);
INCI INCI_OR(INCI x, INCI y);
INCI INCI_XOR(INCI x, INCI y);
int  INCI_EQ(INCI x, INCI y);
int  INCI_LE(INCI x, INCI y);
int  INCI_EQ_0(INCI x);
INCI INCI_0();
INCI INCI_1();
INCI INCI_D2(INCI x);
INCI INCI_PN(INCI x, Long y);
#endif
/*
If we need more bits than can be represented by a single unsigned long long,
these routines are designed to simulate the above definitions.
*/

/*  ============           	S T R U C T U R E S              ============  */

typedef struct {int list[BINLEN];} binary;
/*
A binary number as a list.
*/

typedef struct {int len; int points[MAXNVRTS][POLYDIM];} pointlist;
/*
A polytope state described by the convex hull of a list of points.
pl.len is the number of points and pl.points is the actual list of points.
*/

typedef struct {int len, terminal; float fitness; int bits[MAXNVRTS*POLYDIM*BINLEN];} bitlist;
/*
A bitlist state describing a polytope.
bl.len is the length of the bitlist.
bl.bits is the actual bitlist of 0s and 1s.
bl.fitness is the fitness of the bitlist.
bl.terminal is 1 if the bitlist is a terminal state and 0 otherwise.

*/

typedef struct {int size, nterm; float maxfitness, avfitness; bitlist bl[POPSIZE];} population;
/*
A population of bitlists of size pop.size.
pop.maxfitness and pop.avfitness are the maximum and average fitness of the population respectively.
pop.nterm is the number of terminal states in the population.
pop.bl is the actual list of bitlists.
*/  

typedef struct {int nv; Long x[POLYDIM][VERT_Nmax];} NormalForm;
/* 
The normal form of the vertex matrix.
NF.nv is the number of vertices.
NF.x is the normal form matrix.
*/

typedef struct {int n, np; Long x[POINT_Nmax][POLYDIM];} PolyPointList;
/*
A list of lattice points of a polytope.
P.x[i][j] is the j'th coordinate of the i'th lattice point.
P.n is the dimension of the polytope and P.np the number of points in the list.
*/

typedef struct {int v[VERT_Nmax]; int nv;} VertexNumList;
/*
The list of vertices of a polytope, referring to some PolyPointList P.
The j'th coordinate of the i'th vertex is then given by P.x[V.v[i]][j].
V.nv is the number of vertices of P.
*/

typedef struct {Long a[POLYDIM], c;} Equation;
/*
An equation of the type ax+c=0
*/

typedef struct {int ne; Equation e[EQUA_Nmax];} EqList;
/*
A list of equations. EL.ne is the number of equations in the list.
*/

typedef Long PairMat[EQUA_Nmax][VERT_Nmax];
/*
The matrix whose entries are the pairings av+c between the vertices v and the equations (a,c).
*/

typedef struct {int mp, mv, np, nv, n, pic, cor, h22, h1[POLYDIM-1];} BaHo;                                                                   
/*
This structure is related to Batyrev's formulas for Hodge numbers.
n     ... dimension of the polytope
pic   ... Picard number
cor   ... sum of correction terms
h1[i] ... Hodge number h_{1i}
h22   ... Hodge number h_{22} (if n = 5)
mp, mv, np, nv denote the numbers of points/vertices in the M and N lattices,
repectively.
*/

typedef struct {
    int nf[POLYDIM+1];			                /* #(faces)[dim]  */
 	INCI v[POLYDIM+1][FACE_Nmax]; 		        /*  vertex info   */
 	INCI f[POLYDIM+1][FACE_Nmax]; 		        /* V-on-dual info */
 	Long nip[POLYDIM+1][FACE_Nmax];		        /* #IPs on face  */
 	Long dip[POLYDIM+1][FACE_Nmax];} 	FaceInfo;   /* #IPs on dual  */
/*
nf[i] denotes the number of faces of dimension i
   (the number of faces of dimension n-i-1 of the dual polytope).
v[i][j] encodes the incidence relation of the j'th dim-i face with the vertices
nip[i][j] is the number of interior points of the j'th dim-i face.
f[i][j] and dip[i][j] give the same informations for the dual (n-i-1
   dimensional) faces, with f[i][j] referring to the dual vertices.
*/

typedef	struct 	{int p[SYM_Nmax][VERT_Nmax];}		VPermList;
/*
The vertex permutation list.
*/


/*  ============              B I T L I S T                ============  */


bitlist pts2bts(pointlist pl);
/* 
Convert a pointlist into a bitlist.
*/

pointlist bts2pts(bitlist bl);
/* 
Convert a bitlist into a pointlist.
*/

int randomint(int min, int max);
/*
Generate a random integer in the range from min to max.
*/

int randomchoice(float p[POPSIZE], int len);
/*
Generate a random choice for integers 0,...,len-1 for a probability distribution p.
*/

void fprintbitlist(FILE *fp, bitlist bl);
/*
Write bitlist bl to a file with file pointer fp.
*/

bitlist * freadbitlist(char *filename, int len);
/*
Read pointlists from a file and convert to bitlists.
*/

void flipbit(bitlist *bl, int pos);
/*
Flip bit in position pos for a bitlist bl.
*/

void crossbitlists(bitlist *bl1, bitlist *bl2, int numcuts, int cuts[NUMCUTS]);
/*
Cross bitlists bl1 and bl2, with numcuts number of cuts at positions specified in array cuts.
*/

bitlist randomstate();
/*
Generate a random bitlist state.
*/

int compbitlist(const void *p1, const void *p2);
/*
Compare two bistlist by their fitness. 
*/

NormalForm normalform(bitlist bl);
/*
Compute the vertex matrix normal form of a polytope from its bitlist. 
*/

int NFsequal(NormalForm NF1, NormalForm NF2);
/*
Determine if two normal forms are equal.
*/

int bitlistsequal(bitlist bl1, bitlist bl2);
/*
Determine if two bitlists are equal.
*/

int bitlistsequiv(bitlist bl1, bitlist bl2);
/*
Determine if two bitlists are equivalent.
*/

  
/*  ============              P O P U L A T I O N                ============  */


population randompop(int popsize);
/*
Generate a random population of size popsize.
*/


void nextpop(population pop, population *newpop, int meth, int numcuts, int keepfitest, 
			float mutrate, float alpha);
/* 
Generate the next population from a given one by performing selection, crossover and mutation. 
*/


/*  ============              E V O L U T I O N                ============  */


population * evolvepop(population initialpop, int numgen, int meth, int numcuts,
			      int keepfitest, float mutrate, float alpha, int monitor);
/*
Evolve a given population over numgen generations.
*/

void monitorevol(int gen, population *pop);
/*
Monitor the evolution of a population.
*/

bitlist * termstates(population *evol, int numgen, int *numterm);
/*
Extracts the terminal states from a population.
*/

void removeredundancy(bitlist *bl, int *len);
/*
Removes the redundancy in a list of bitlists bl of length len.
The number of reduced states is assigned to len and the bl is sorted such that the 
first len bitlists comprise the reduced list.
*/

bitlist * termstatesred(population *evol, int numgen, int *numterm);
/*
Extracts the terminal states from a population and removes the redundancy.
*/   
			   
bitlist * searchenv(int numrun, int numevol, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, int *numterm);
/*
Repeatedly evolves random initial populations extracts terminal states and saves them to a file.
Evolves numrun*numevol populations of size popsize over numgen generations.
*/


/*  ============             F I T N E S S               ============  */


void fitness(bitlist * bl);
/*
Computes the fitness of the bitlist bl and assigns the result to bl.fitness.
*/


/*  ============             V E R T E X                ============  */


void Sort_VL(VertexNumList *V);
/*
Sorts the entries _V->v[i] in ascending order.
*/

Long Eval_Eq_on_V(Equation *E, Long *V, int n);
/*
Evaluates E on V, i.e. calculates \sum_{i=0}^{n-1} E->a[i] * V[i] + E->c.
*/

int  Vec_Greater_Than(Long *X, Long *Y, int n);
/*
Returns 1 if *X > *Y in the sense that X[i] > Y[i] for the first i where
X[i] and Y[i] differ, returns 0 if *X < *Y and gives an error message if
X[i] equals Y[i] for all i in {0,...n-1}.
*/

Equation EEV_To_Equation(Equation *E1, Equation *E2, Long *V, int n);
/*
Returns the equation describing the span of the vector V and the intersection
of the hyperplanes corresponding to E1 and E2; n is the dimension.
*/

void Make_VEPM(PolyPointList *P, VertexNumList *VNL, EqList *EL, PairMat PM);
/*
Calculates the matrix of pairings between the vertices in VNL and the
equations in EL.
*/

void EL_to_PPL(EqList *EL, PolyPointList *DP, int *n);
/*
Converts *EL to the incomplete PolyPointList *DP corresponding to the dual
polytope; *n is the dimension. 
*/

int  Find_Equations(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
For the polytope determined by P, *VNL and *EL are calculated.
*VNL is the complete list of vertices of P.
*EL is the complete list of equations determining the facets of P.
Find_Equations returns 1 if P has IP property (i.e., it has the
origin in its interior) and 0 otherwise.
*/

void Complete_Poly(Long VPM[][VERT_Nmax],EqList *E,int nv,PolyPointList *P);
/*
Given the vertex pairing matrix VPM, the EqList *E and the number nv of
vertices, the complete list of lattice points *P is determined.
*/

INCI Eq_To_INCI(Equation *E, PolyPointList *P, VertexNumList *VNL);
/*
Converts *E to an INCI.
*/

void Make_Incidence(PolyPointList *P, VertexNumList *VNL, EqList *EL, FaceInfo *FI);
/*
Creates the structure FaceInfo *FI from *P, *VNL and *EL.
*/

int QuickAnalysis(PolyPointList *_P, BaHo *_BH, FaceInfo *_FI);
/*
Fast computation of FaceInfo and Hodge numbers.
*/


/*  ============             N O R M A L   F O R M                ============  */


void swap(int *i,int *j);
/*
Swaps *i and *j.
*/

void  Make_Poly_Sym_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F, Long NF[POLYDIM][VERT_Nmax]);
/*
Given *P, *VNL and *EL the normal form coordinates NF of the vertice is computed.
*/
