/*  ======================================================================  */
/*  ==========	     			   	  	   	==========  */
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


#define	Long long long    
#define LLong long long


#define MIN -15		  /* lower bound on integers considered for points */
#define MAXNVRTS 6  	  /* maximum number of vertices */
#define POLYDIM 5	  /* dimension of polytopes */
#define BINLEN 5	  /* maximum length of binary number */
#define POPSIZE 500       /* population size */
#define NUMGEN 500        /* number of generations */
#define NUMCUTS 1         /* number of cuts in a crossing */
#define METHOD 1          /* macro for ranking (0) or roulette (1) method to select breeding pairs */
#define MUTRATE 0.005     /* mutation rate */
#define ALPHA 3.0         /* value of alpha parameter in selection probability */
#define KEEPFITEST 1      /* copy fitest individual into next population (1) or don't (0) */
#define MONITOR 1         /* terminal monitor on (1) or off (0) */

/* fitness weights */
#define DIST_WEIGHT 1 
#define IP_WEIGHT 1 
#define NVERTS_WEIGHT 0
#define NVERTS 6
#define NPTS_WEIGHT 0 /* if >0 slows down computation */
#define NPTS 7
#define H11_WEIGHT 0 /* if >0 slows down computation */
#define H11 1
#define H12_WEIGHT 0 /* if >0 slows down computation */
#define H12 1
#define H13_WEIGHT 0 /* if >0 slows down computation */
#define H13 1
#define H22_WEIGHT 0 /* if >0 slows down computation */
#define H22 1
#define EULER_WEIGHT 0 /* if >0 slows down computation */
#define EULER 1


#define POINT_Nmax 2000000         /* maximum number of points */
#define VERT_Nmax 128	           /* maximum number of vertices */
#define FACE_Nmax 10000	           /* maximum number of faces */
#define SYM_Nmax 46080             /* maximum number of symmetries */
#define EQUA_Nmax VERT_Nmax        /* maximum number of facets */		
#define AMBI_Dmax (5 * POLYDIM)    /* maximum dimension of the ambient space */

#define FIB_Nmax 3000 /* maximum number of allowed weight relations among points that define the IP simplices */

#define GL_Long	Long /* uses W_to_GLZ like in Rat.c */

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
Binary list object.
*/

typedef struct {int len; int points[MAXNVRTS][POLYDIM];} pointlist;
/*
Point list of polytope of length len.
*/

typedef struct {int len; int bits[MAXNVRTS*POLYDIM*BINLEN]; float fitness; int terminal;} bitlist;
/*
Bitlist data object.
bl.bits is the actual bitlist of length bl.len.
bl.fitness is the fitness of the corresponding polytope.
bl.terminal is 1 is the state is terminal and 0 otherwise.
*/

typedef struct {int size; float maxfitness; float avfitness; int nterm; bitlist bl[POPSIZE];} population;  
/* 
Population data object.
pop.bl is the list of bitlists of length pop.size.
pop.maxfitness and pop.avfitness is the maximum and average fitness respectively.
pop.nterm is the number of terminal states in the population.
*/

typedef struct {int nv; Long x[POLYDIM][VERT_Nmax];} NormalForm; 
/* 
Data object for the normal form of a polytope.
NF.x is the d x Nv normal form matrix and NF.nv is the number of vertices, i.e. Nv.
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
convert a points list into a bit list after adding max 
*/

pointlist bts2pts(bitlist bl);
/* 
convert a bit list into a points list and subtract max 
*/

int randomint(int min, int max);
/* 
generate a random integer in a range from min to max 
*/

int randomchoice(float p[POPSIZE], int len);
/* 
generate a random choice for integers 0,...,len-1 for a probability distribution p 
*/

void fprintbitlist(FILE *fp, bitlist bl);
/* 
write a bitlist bl to a file with file pointer fp 
*/

bitlist * freadbitlist(char *filename, int len);
/* 
read bitlists from a file 
*/

void flipbit(bitlist *bl, int pos);
/* 
flips bit in position pos for a bitlist bl 
*/

void crossbitlists(bitlist *bl1, bitlist *bl2, int numcuts, int cuts[NUMCUTS]);
/* 
crosses bitlists bl1 and bl2, with numcuts number of cuts at positions specified in array cuts 
*/

bitlist randomstate();
/* 
the randomstate function which generates a random state 
*/

int compbitlist(const void *p1, const void *p2);
/* 
compare two bistlist by their fitness
*/

NormalForm normalform(bitlist bl);
/* 
compute the normal form of a polytope from its bitlist 
*/

int NFsequal(NormalForm NF1, NormalForm NF2);
/* 
decide if two normal forms are equal 
*/

int bitlistsequal(bitlist bl1, bitlist bl2);
/* 
decide if two bistlist are identical 
*/

int bitlistsequiv(bitlist bl1, bitlist bl2);
/* 
decide if two bistlist describe equivalent polytopes 
*/
  
  
/*  ============              P O P U L A T I O N                ============  */


population randompop(int popsize);
/* 
generate a random population of size popsize 
*/

float avfitness(population *pop);
/* 
find average fitness of a population 
*/

void nextpop(population pop, population *newpop, int meth, int numcuts, int keepfitest,
	     float mutrate, float alpha);
/* 
generate the next population from a given one 
*/


/*  ============              E V O L U T I O N                ============  */


population * evolvepop(population initialpop, int numgen, int meth, int numcuts,
			      int keepfitest, float mutrate, float alpha, int monitor);
/* 
generate the next population from a given one 
*/

bitlist * termstates(population *evol, int numgen, int *numterm);
/* 
select terminal states from a population 
*/

void removeredundancy(bitlist *bl, int *len);
/* 
remove redundqncy in list of bitlists 
*/

bitlist * termstatesred(population *evol, int numgen, int *numterm);
/* 
select terminal states from a population and remove redundancy 
*/
			   
bitlist * searchenv(int numrun, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, int *numterm);
/* 
repeated evolution of a random initial population, extracting terminal states 
*/
			   

/*  ============             F I T N E S S               ============  */


void fitness(bitlist * bl);
/* 
the fitness function which updates the fitness of a bitlist 
*/


/*  ============             V E R T E X                ============  */


Long Eval_Eq_on_V(Equation *E, Long *V, int n);
/*
Evaluates E on V, i.e. calculates \sum_{i=0}^{n-1} E->a[i] * V[i] + E->c.
*/

void Make_VEPM(PolyPointList *P, VertexNumList *VNL, EqList *EL, PairMat PM);
/*
Calculates the matrix of pairings between the vertices in VNL and the
equations in EL.
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

int QuickAnalysis(PolyPointList *_P, BaHo *_BH, FaceInfo *_FI);
/*
Fast computation of FaceInfo and Hodge numbers.
*/


/*  ============             N O R M A L   F O R M                ============  */


void  Make_Poly_Sym_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F, Long NF[POLYDIM][VERT_Nmax]);
/*
Given *P, *VNL and *EL, the normal form coordinates NF of the vertices is computed.
*/
