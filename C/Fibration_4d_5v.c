/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ==========       F I B R A T I O N   F U N C T I O N S      ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

#include "Global_4d_5v.h"

/* find the largest element in a list of integers */
Long max(int len, Long x[len]){
  int i;
  Long maximum = x[0];
  
  for (i=1; i<len; i++){
    if (x[i] > maximum) maximum = x[i];
  }
  
  return maximum;
}

/* compute the dot product */
Long dotprod(int dim, Long v1[dim], Long v2[dim])
{
  int k;
  Long sum = 0;
  
  for (k=0; k<dim; k++){
    sum += v1[k]*v2[k];
  }
  
  return sum;
}

/* identify points with small inner product with all vertices */
void findshort(int p, int v, int dim, Long points[p][dim], Long vertices[v][dim], struct slist * sList, struct shortlist * shortList)
{
  int i,j,k,m;
  Long dots[v];
  
  for (i=0; i<p; i++){
    for (j=0; j<v; j++){
      dots[j] = dotprod(dim, points[i], vertices[j]);
    }
  	m = max(v,dots);
  	if (m<4 && m!=0){
  	  sList->list[m-1] += 1;
  	  if (m == 1){
  	     for (k=0; k<dim; k++) shortList->S1[sList->list[m-1]][k] = points[i][k];
  	  }
  	  else if (m == 2){
  	    for (k=0; k<dim; k++) shortList->S2[sList->list[m-1]][k] = points[i][k];
  	  }
  	  else{
  	    for (k=0; k<dim; k++) shortList->S3[sList->list[m-1]][k] = points[i][k];
  	  }
  	}
  }
} 

/* check if a vector exists in a list */
int member(int dim1, int dim2, Long v[dim2], Long list[dim1][dim2]){
  int i, j, test;
 
  for (i=0; i<dim1; i++){
  	test = 1;
    for (j=0; j<dim2; j++){
      if (v[j] != list[i][j]) test = 0;
    }
    if (test){
      return 1;
    }
  }
  return 0;
}

/* check if two vectors are equal */
int vecequal(int dim, Long v1[dim], Long v2[dim]){
  int i;
 
  for (i=0; i<dim; i++){
  	if (v1[i] != v2[i]) return 0;
  }
  	
  return 1;
}

/* look for all types of fibres */
int fiber(struct shortlist * shortList, struct slist * sList, int dim)
{
  int i,j,k;
  Long x[dim], y[dim], z[dim];
  
  for (i=0; i<sList->list[0]; i++){
    for (k=0; k<dim; k++) x[j] = -1*shortList->S1[i][k];
    if (member(sList->list[0], dim, x, shortList->S1)){
      for (j=i+1; j<sList->list[0]; j++){
      for (k=0; k<dim; k++) y[k] = -1*shortList->S1[j][k];
        if (member(sList->list[0], dim, y, shortList->S1)){
          if (vecequal(dim, shortList->S1[i], y)) return 1;
        }
      }
    } 
  }
  
  for (i=0; i<sList->list[1]; i++){
    for (j=i+1; j<sList->list[1]; j++){
      for (k=0; k<dim; k++) z[k] = -1*shortList->S2[i][k] -1*shortList->S2[j][k];
      if (member(sList->list[0], dim, z, shortList->S1) || 
          member(sList->list[1], dim, z, shortList->S2)) return 2;
    }
  }
  
  for (i=0; i<sList->list[2]; i++){
    for (j=0; j<sList->list[2]; j++){
      for (k=0; k<dim; k++) z[k] = (-1*shortList->S3[i][k] -1*shortList->S3[j][k])/2;
      if (member(sList->list[0], dim, z, shortList->S1)) return 3; 
    }
  }
  
  return 0;
}

/* determine if there exists an elliptic fibration */
int fibration(PolyPointList *_P, VertexNumList *_V)
{
  int i,j;
  Long points[_P->np][POLYDIM];
  Long vertices[_V->nv][POLYDIM];
  
  for(i=0; i<_P->np; i++) for(j=0; j<POLYDIM; j++) points[i][j] = (_P->x)[i][j];
  for(i=0; i<_V->nv; i++) for(j=0; j<POLYDIM; j++) vertices[i][j] = (_P->x)[_V->v[i]][j];

  struct shortlist * shortList = (struct shortlist *) malloc(sizeof(struct shortlist));
  struct slist * sList = (struct slist *) malloc(sizeof(struct slist));
  int x;

  findshort(_P->np, _V->nv, POLYDIM, points, vertices, sList, shortList);

  x = fiber(shortList, sList, _P->n); 
  
  free(shortList);free(sList);
  
  if (x > 0) return 1; 
  
  return 0;
}

