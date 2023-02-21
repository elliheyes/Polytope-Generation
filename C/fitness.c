/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ==========         F I T N E S S   F U N C T I O N          ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

#include "Global.h"

/* compute the fitness of a polytope from its bitlist */
void fitness(struct bitlist *bl)
{
  float score;
  int IP;
  int i,j,k;
  int range=pow(2,BINLEN)-1;
  int interior;
  float numInterior;
  float totalDist, avDist;
  int h11, h12, h13, h22, euler;
  
  /* compute the list of points from the bitlist */
  struct pointlist pl = bts2pts(*bl); 

  if(pl.len < POLYDIM) score = -100;
  else{
    /* define the data objects */
  	VertexNumList V;
  	EqList *E = (EqList *) malloc(sizeof(EqList));
  	PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList));
  	PairMat *PM = (PairMat *) malloc(sizeof(PairMat));
  	BaHo BH;
  	FaceInfo *FI = (FaceInfo *) malloc(sizeof(FaceInfo));
  
  	/* define the polytope dimension */
  	_P->n=POLYDIM; 
  	
  	/* define the number of points */
  	_P->np=pl.len; 
  
  	/* define the points */
  	for(i=0; i<pl.len; i++) for(j=0; j<POLYDIM; j++) _P->x[i][j]=pl.points[i][j];

    /* find the bounding hyperplane equations of the polytope */
    IP=Find_Equations(_P,&V,E); 
    
    /* find the vertex pairing matrices */
  	Make_VEPM(_P,&V,E,*PM); 

  	/* find the complete lists of points */
  	Complete_Poly(*PM,E,V.nv,_P,&V); 
  	
  	/* compute the fitness score */
  	score = 0;
  	
  	/* penalty for the IP property */
  	if(IP_WEIGHT>0) score += IP_WEIGHT*(IP-1);
  	
  	/* penalty for number of interior points */
  	if(INTERIOR_WEIGHT>0){
  		numInterior = 0.;
  		for(i=0; i<_P->np; i++){
  			interior = 1;
  			for(j=0; j<E->ne; j++){
  				if(Eval_Eq_on_V(&E->e[j],_P->x[i],_P->n) == 0) interior=0;
  			}
  			if(interior) numInterior++;
  		}
  		score += -INTERIOR_WEIGHT*abs(1-numInterior)/(range*POLYDIM);
  	}
  	
  	/* penalty for the distance of facets from the origin */
  	if(DIST_WEIGHT>0){
    	totalDist = 0.;
    	for(i=0; i<E->ne; i++) totalDist += llabs(E->e[i].c-1);
    	avDist = totalDist/E->ne;
    	score += -DIST_WEIGHT*avDist/(range*POLYDIM);
    }

	/* penalty for the number of vertices */
	if(NVERTS_WEIGHT>0) score += -NVERTS_WEIGHT*abs(V.nv-NVERTS);
	
	/* penalty for the hodge numbers */
	if(H11_WEIGHT > 0 || H12_WEIGHT > 0 || H13_WEIGHT > 0 || H22_WEIGHT > 0 || EULER_WEIGHT > 0){
		if(score == 0){
			QuickAnalysis(_P, &BH, FI);
			if(H11_WEIGHT > 0) score += -H11_WEIGHT*abs(BH.h1[1]-H11);
			if(H12_WEIGHT > 0) score += -H12_WEIGHT*abs(BH.h1[2]-H11);
			if(H13_WEIGHT > 0) score += -H13_WEIGHT*abs(BH.h1[3]-H11);
			if(H22_WEIGHT > 0) score += -H22_WEIGHT*abs(BH.h22-H22);
			if(EULER_WEIGHT > 0) score += -EULER_WEIGHT*abs(6*(8+BH.h1[1]+BH.h1[3]-BH.h1[2])-EULER);
		}
	}
  
    /* free allocated memory */
    free(E);free(_P);free(PM);free(FI);
  }
  
  /* update bitlist fitness and terminal */
  bl->fitness = score;
  
  if(bl->fitness == 0) bl->terminal = 1;
  else bl->terminal = 0;

}


