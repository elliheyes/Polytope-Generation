/*  ======================================================================  */
/*  ==========	     			   	  	        ==========  */
/*  ==========         F I T N E S S   F U N C T I O N          ==========  */
/*  ==========						        ==========  */
/*  ======================================================================  */

#include "Global.h"

/* exchange two rows of a matrix */
void swapRows(int R, int C, int mat[R][C], int row1, int row2, int col)
{
	int i;
    for (i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}
 
/* compute the rank of a matrix */
int rankOfMatrix(int R, int C, int mat[R][C])
{
    int row, col, i, rank = C;
 
    for (row = 0; row < rank; row++)
    { 
        if (mat[row][row])
        {
           for (col = 0; col < R; col++)
           {
               if (col != row)
               {
                 double mult = (double)mat[col][row] / mat[row][row];
                 for (i = 0; i < rank; i++) mat[col][i] -= mult * mat[row][i];
              }
           }
        }
        else
        {
            int reduce = 1;
 
            for (i = row + 1; i < R;  i++)
            {
                if (mat[i][row])
                {
                    swapRows(R, C, mat, row, i, rank);
                    reduce = 0;
                    break ;
                }
            }
 
            if (reduce)
            {
                rank--;
 
                for (i = 0; i < R; i ++)
                    mat[i][row] = mat[i][rank];
            }
            row--;
        }
    }
    return rank;
}

/* compute the dimension of a polytope from its pointlist */
int dimensionOfPolytope(pointlist pl)
{
    int i, j, mat[pl.len][POLYDIM];
    
    for(i=0; i<pl.len; i++){
        for(j=0; j<POLYDIM; j++){
            mat[i][j] = pl.points[i][j]-pl.points[0][j];
        }
    }
    
    int dim = rankOfMatrix(pl.len, POLYDIM, mat);
    
    return dim;
}

/* compute the fitness of a polytope from its bitlist */
void fitness(bitlist *bl)
{
  int dim;
  float score;
  pointlist pl;
  
  /* compute the list of points from the bitlist */
  pl = bts2pts(*bl); 
  
  /* compute the dimension of the polytope */
  dim = dimensionOfPolytope(pl);
  
  /* if the dimension of the polytope is less than POLYDIM then assign a large penalty */
  if(dim < POLYDIM) score = -10;
  else{
    int i, j, k, IP, interior, h11, h12, h13, h22, euler, range=pow(2,BINLEN)-1;
  	float numInterior, totalDist, avDist;
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
  	
  	/* compute the fitness score */
  	score = 0;
  	
  	/* penalty for the IP property */
  	if(IP_WEIGHT>0) score += IP_WEIGHT*(IP-1);
  	
  	/* penalty for the distance of facets from the origin */
  	if(DIST_WEIGHT>0){
    	totalDist = 0.;
    	for(i=0; i<E->ne; i++) totalDist += llabs(E->e[i].c-1);
    	avDist = totalDist/E->ne;
    	score += -DIST_WEIGHT*avDist/(range*POLYDIM);
    }

	/* penalty for the number of vertices */
	if(NVERTS_WEIGHT>0) score += -NVERTS_WEIGHT*abs(V.nv-NVERTS);
	
	if(score==0 && (NPTS_WEIGHT > 0 || EULER_WEIGHT > 0 || H11_WEIGHT > 0 || 
		H12_WEIGHT > 0 || H13_WEIGHT > 0 || H22_WEIGHT > 0 )){
		/* find the vertex pairing matrices */
  		Make_VEPM(_P,&V,E,*PM); 
  	
  		/* find the complete lists of points */
  		Complete_Poly(*PM,E,V.nv,_P); 
  		
  		/* penalty for the number of points */
		if(NPTS_WEIGHT>0) score += -NPTS_WEIGHT*abs(_P->np-NPTS);
	
		/* penalties for topological data */
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
