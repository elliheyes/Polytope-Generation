/*  ======================================================================  */
/*  ==========	     			                     	==========  */
/*  ==========                    G L O B A L                   ==========  */
/*  ==========				                        ==========  */
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


#define MIN -10	          /* lower bound on integers considered for points */
#define MAXNVRTS 6  	  /* maximum number of vertices */
#define POLYDIM 5         /* dimension of polytopes */
#define BINLEN 4          /* maximum length of binary number */
#define POPSIZE 150       /* population size */
#define NUMGEN 300        /* number of generations */
#define NUMCUTS 1         /* number of cuts in a crossing */
#define RANKING 0         /* macro for ranking method to select breeding pairs */
#define ROULETTE 1        /* macro for roulette method to select breading pairs */
#define MUTRATE 0.005     /* mutation rate */
#define ALPHA 3.0         /* value of alpha parameter in selection probability */
#define KEEPFITEST 1      /* copy fitest individual into next population */
#define DONTKEEPFITEST 0  /* don't copy over fitest individual */
#define MONITORON 1       /* terminal monitor on */
#define MONITOROFF 0      /* terminal monitor off */

/* fitness weights */
#define DIST_WEIGHT 1 
#define IP_WEIGHT 1 
#define INTERIOR_WEIGHT 0 /* if >0 slows down computation */
#define NVERTS_WEIGHT 0
#define NVERTS 6
#define NPTS_WEIGHT 1
#define NPTS 7
#define H11_WEIGHT 0
#define H11 1
#define H12_WEIGHT 0
#define H12 1
#define H13_WEIGHT 0
#define H13 1
#define H22_WEIGHT 0
#define H22 1
#define EULER_WEIGHT 0
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

struct binary
{
  int list[BINLEN];  /* the actual binary list */
};

struct pointlist
{
  int len;                    /* the number of points */
  int points[MAXNPTS][POLYDIM];  /* the actual point list */
};

struct bitlist
{
  int len;                        /* number of bits in bitlist */
  int bits[MAXNPTS*POLYDIM*BINLEN];  /* the actual bitlist */
  float fitness;                  /* fitness of bitlist */
  int terminal;                   /* terminal or not */
};

struct population
{
  int size;                    /* size of population */
  float maxfitness;            /* maximal fitness in population */
  float avfitness;             /* average fitness in population */
  int nterm;                   /* number of terminal states in population */
  struct bitlist bl[POPSIZE];  /* the actual population */
};  

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

typedef struct {int mp, mv, np, nv, n, pic, cor, h22, h1[POLYDIM-1];} BaHo;                                                                    BaHo;
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



/* convert a decimal number into a binary list */
struct binary decimal2binary(int num);

/* convert a binary list into a decimal number */
int binaryToDecimal(struct binary bin);

/* convert a points list into a bit list after adding max */
struct bitlist pts2bts(struct pointlist pl);

/* convert a bit list into a points list and subtract max */
struct pointlist bts2pts(struct bitlist bl);

/* generate a random bitlist of lenght len */
struct bitlist randombitlist(int len);

/* generate a random integer in a range from min to max */
int randomint(int min, int max);

/* generate a random choice for integers 0,...,len-1 for a probability distribution p */
int randomchoice(float p[POPSIZE], int len);

/* write a bitlist bl to a file with file pointer fp */
void fprintbitlist(FILE *fp, struct bitlist bl);

/* flips bit in position pos for a bitlist bl */
void flipbit(struct bitlist *bl, int pos);

/* copies the bits in bitlist bl from position pos1 to pos2-1 into a new bitlist */
struct bitlist copybitlist(struct bitlist bl, int pos1, int pos2);

/* pastes bitlist blpaste into bitlist bl at position pos, overwriting previous content in bl */
void pastebitlist(struct bitlist *bl, struct bitlist blpaste, int pos);

/* a simple insertion sort of an integer array into ascending order */
void isort(int arr[], int len);

/* crosses bitlists bl1 and bl2, with numcuts number of cuts at positions specified in array cuts */
void crossbitlists(struct bitlist *bl1, struct bitlist *bl2, int numcuts, int cuts[NUMCUTS]);

/* the randomstate function which generates a random state */
struct bitlist randomstate();

/* compare two bistlist by their fitness */
int compbitlist(const void *p1, const void *p2);

/* decide if two bistlist are identical */
int bitlistsequal(struct bitlist bl1, struct bitlist bl2);

/* decide if two bistlist describe equivalent polytopes */
int bitlistsequiv(struct bitlist bl1, struct bitlist bl2);
  
  
  
/*  ============              P O P U L A T I O N                ============  */



/* generate a random population of size popsize */
struct population randompop(int popsize);

/* sort a population by fitness */
void sortpop(struct population *pop);

/* mutate a population */
void mutatepop(struct population *pop, float mutrate);

/* find fitest in a population */
struct bitlist fitestinpop(struct population *pop);

/* find average fitness of a population */
float avfitness(struct population *pop);

/* find number of terminal state in population */
int nterminalpop(struct population *pop);

/* generate the next population from a given one */
void nextpop(struct population pop, struct population *newpop, int meth, int numcuts, int keepfitest,
	     float mutrate, float alpha);



/*  ============              E V O L U T I O N                ============  */



/* genetically evolve a population */
struct population * evolvepop(struct population initialpop, int numgen, int meth, int numcuts,
			      int keepfitest, float mutrate, float alpha, int monitor);

/* monitor evolution of a population */
void monitorevol(int gen, struct population *pop);

/* select terminal states from a population */
struct bitlist * termstates(struct population *evol, int numgen, int *numterm);

/* remove redundqncy in list of bitlists */
void removeredundancy(struct bitlist *bl, int *len);

/* select terminal states from a population and remove redundancy */
struct bitlist * termstatesred(struct population *evol, int numgen, int *numterm);

/* repeated evolution of a random initial population, extracting terminal states */
struct bitlist * searchenv(int numevol, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, FILE * fp, int *numterm);



/*  ============             F I T N E S S               ============  */



/* the fitness function which returns the fitness of a bitlist */
void fitness(struct bitlist * bl);



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

void Make_Incidence(PolyPointList *P, VertexNumList *VNL, EqList *EL,
                    FaceInfo *FI);
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

int  Make_Poly_Sym_NF(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		      int *SymNum, int V_perm[][VERT_Nmax], Long NF[POLYDIM][VERT_Nmax]);
/*
Given *P, *VNL and *EL, the following objects are determined:
the number *SymNum of GL(n,Z)-symmetries of the polytope,
the *SymNum vertex permutations V_perm realising these symmetries,
the normal form coordinates NF of the vertices,
the number of symmetries of the vertex pairing matrix
    (this number is the return value of Make_Poly_Sym_NF).
*/
