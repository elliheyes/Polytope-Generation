/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ==========         F I T N E S S   F U N C T I O N          ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

#include "Global.h"

/* compute the fitness of a polytope from its bitlist */

void fitness(struct bitlist *bl)
{
  float penalties[2];
  float score;
  int IP;
  int i,j,k;
  int interior,numInterior;
  int totalDist;
  
  /* compute the list of points from the bitlist */
  struct pointlist pl = bts2pts(*bl); 
  
  /* define the data objects */
  VertexNumList V;
  EqList *E = (EqList *) malloc(sizeof(EqList));
  PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList));
  PairMat *PM = (PairMat *) malloc(sizeof(PairMat));
  
  /* define the number of points and polytope dimension */
  _P->np=NPTS; 
  _P->n=POLYDIM; 
  
  /* define the points */
  for(i=0; i<NPTS; i++){
    for(j=0; j<POLYDIM; j++){
      _P->x[i][j]=pl.points[i][j];
    }
  } 
	
  /* find the bounding hyperplane equations of the polytope */
  IP=Find_Equations(_P,&V,E); 
  
  /* find the vertex pairing matrix */
  Make_VEPM(_P,&V,E,*PM); 

  /* find the complete list of points */
  Complete_Poly(*PM,E,V.nv,_P); 

  /* determine the interior points */
  numInterior = 0;
  for(i=0; i<_P->np; i++){
    interior = 1;
    for(j=0; j<E->ne; j++){
      if(Eval_Eq_on_V(&E->e[j],_P->x[i],_P->n)==0) interior = 0;
    }
    if(interior) numInterior += 1;
  }
  
  /* define penalty for the number of interior points */
  if(numInterior == 1) penalties[0] = 0;
  else if(numInterior == 0) penalties[0] = -1;
  else penalties[0] = -numInterior;
  
  /* determine the total distance of the faces from the origin */
  totalDist = 0.;
  for(i=0; i<E->ne; i++){
  	totalDist += llabs(llabs(E->e[i].c)-1);
  }

  /* define penalty for the distance of the faces */
  penalties[1] = -totalDist;
  
  /* compute the total fitness score */
  score = penalties[0] + penalties[1];
  
  /* update bitlist fitness and terminal */
  bl->fitness = score;
  
  if(bl->fitness == 0) bl->terminal = 1;
  else bl->terminal = 0;
  
}




















