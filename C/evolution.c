/*  ======================================================================  */
/*  ==========	     			   	  	     	==========  */
/*  ==========       E V O L U T I O N   F U N C T I O N S      ==========  */
/*  ==========						        ==========  */
/*  ======================================================================  */

#include "Global.h"


/* genetically evolve a population */
struct population * evolvepop(struct population initialpop, int numgen, int meth, int numcuts,
			      int keepfitest, float mutrate, float alpha, int monitor)
{
  int gen, nterm;
  struct population *evol;

   /* check validity of numgen */
  if (numgen>NUMGEN) numgen=NUMGEN;
  
  /* allocate memory for evolution */
  evol = calloc(numgen,sizeof(struct population));
  if (evol == NULL) {
    printf("evolvepop: memory allocation failed");
    exit(0);
  }  

  /* load initial population */
  evol[0]=initialpop;
  if (monitor) monitorevol(0,&(evol[0])); 

  /* main loop over generations */
  nterm=initialpop.nterm;
  for (gen=0; gen<numgen-1; gen++) {
    nextpop(evol[gen],&(evol[gen+1]),meth,numcuts,keepfitest,mutrate,alpha);  /* determine next generation */
    nterm=nterm+evol[gen+1].nterm;  /* count number of terminal states */
    if (monitor) monitorevol(gen+1,&(evol[gen+1]));  /* monitor */ 
  }

  if (monitor)printf("\nNumber of terminal states: %i\n",nterm);
  
  return evol;
}  


/* monitor evolution of a population */
void monitorevol(int gen, struct population *pop)
{
  if (!(gen%10)) printf("Gen    AvFit     MaxFit    #Term\n");
  printf("%3i    %2.4f    %2.4f    %3i\n",gen,pop->avfitness,pop->maxfitness,pop->nterm);
  fflush(stdout);
}  


/* select terminal states from a population */
struct pointlist * termstates(struct population *evol, int numgen, int *numterm)
{
  int i,j, k;
  struct pointlist *pl;
  
  /* find the number of terminal states */
  *numterm=0;
  for (k=0; k<numgen; k++)
    for (i=0; i<evol[k].size; i++) 
      if ((evol[k].bl)[i].terminal){
        (*numterm)++;
      }
 
   /* allocate memory for terminal states */
  pl = calloc(*numterm,sizeof(struct pointlist)); 

  /* store terminal states in allocated space */
  j=0;
  for (k=0; k<numgen; k++)
    for (i=0; i<evol[k].size; i++) 
      if ((evol[k].bl)[i].terminal){
	    pl[j]=bts2pts((evol[k].bl)[i]); 
	    j++;
      }
				 
   return pl;
}


/* remove equality and equivalence redundancy in list of bitlists */
void removeredundancy(struct pointlist *pl, int *len)
{
  /* remove equality redundancy */
  int cnonred, cactive, red, k; 
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (pointlistsequal(pl[k],pl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    pl[cnonred]=pl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(pl,cnonred,sizeof(struct pointlist),compbitlist);
  }
  
  /* remove equivalence redundancy */
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (pointlistsequiv(pl[k],pl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    pl[cnonred]=pl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(pl,cnonred,sizeof(struct pointlist),compbitlist);
  }
}


/* select terminal states from a population add mirror duals and remove redundancy */
struct pointlist * termstatesred(struct population *evol, int numgen, int *numterm)
{
  int i,j,k,IP;
  struct pointlist pl, *pl1, *pl2;

  /* extract terminal states */
  pl1=termstates(evol,numgen,numterm);
  
  /* add mirror duals */
  pl2 = calloc(2*(*numterm),sizeof(struct pointlist));
  for(i=0; i<*numterm; i++){
  	pl2[i*2] = pl1[i];
  
    VertexNumList V;
    EqList *E = (EqList *) malloc(sizeof(EqList));
    PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList));

    _P->n=POLYDIM; 
    _P->np=pl1[i].len; 
    for(j=0; j<pl1[i].len; j++) for(k=0; k<POLYDIM; k++) _P->x[j][k]=pl1[i].points[j][k];
      
    IP=Find_Equations(_P,&V,E);
      
    pl.len=E->ne;
    for(j=0; j<E->ne; j++) for(k=0; k<POLYDIM; k++) pl.points[j][k]=E->e[j].a[k];
      
    pl2[i*2+1]=pl;
      
    free(E);free(_P);
  } 
  
  *numterm = 2*(*numterm);
  
  free(pl1);
    
  /* remove redundancy in the list of terminal states */
  removeredundancy(pl2,numterm);
  
  return pl2;
}  


/* repeated evolution of a random initial population, extracting terminal states */
struct pointlist * searchenv(int numrun, int numevol, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, FILE * fp, int *numterm)
{
  int i, j, k, IP, crun, cevol, n1, n2, nterm1, nterm2;
  struct population *evol;
  struct pointlist *pl, *plterm1, *plterm2, *pltermOld1, *pltermOld2;

  /* main loop over runs */
  nterm1=0;
  for(crun=0; crun<numrun; crun++){
  
    /* loop over evolutions */
    nterm2=0;
    for(cevol=0; cevol<numevol; cevol++){
    
      /* evolve random population */
      evol=evolvepop(randompop(popsize),numgen,meth,numcuts,keepfitest,mutrate,alpha,0);
      
      /* extract terminal states add mirror duals and remove redundancy */
      pl=termstatesred(evol,numgen,&n2);
      nterm2=nterm2+n2;
      
      /* allocate or re-allocate memory for terminal states */
      if(cevol==0) {
        plterm2=calloc(nterm2,sizeof(struct pointlist));
      }
      else{
        if(n2!=0){
            pltermOld2=calloc(nterm2 - n2,sizeof(struct pointlist));;
            for (i=0; i< nterm2 - n2; i++) pltermOld2[i]=plterm2[i];
            free(plterm2);
            plterm2=calloc(nterm2,sizeof(struct pointlist));
            for (i=0; i< nterm2 - n2; i++) plterm2[i]=pltermOld2[i];
            free(pltermOld2);
        }
      }
    
      /* load in new terminal states */
      for (i=0; i<n2; i++) plterm2[nterm2-n2+i]=pl[i];

	  /* free allocated memory */ 
      free(pl);free(evol);
    }
    n1=nterm2;
    nterm1=nterm1+n1;
    
    /* allocate or re-allocate memory for terminal states */
    if(crun==0){
      plterm1=calloc(nterm1,sizeof(struct pointlist));
    }
    else{
      if(n1!=0){
        pltermOld1=calloc(nterm1 - n1,sizeof(struct pointlist));;
        for (i=0; i< nterm1 - n1; i++) pltermOld1[i]=plterm1[i];
        free(plterm1);
        plterm1=calloc(nterm1,sizeof(struct pointlist));
        for (i=0; i< nterm1 - n1; i++) plterm1[i]=pltermOld1[i];
        free(pltermOld1);
      }
    }
    
    /* load in new terminal states */
    for (i=0; i<n1; i++) plterm1[(nterm1)-n1+i]=plterm2[i];
      
    /* remove redundancy on set */
    removeredundancy(plterm1,&nterm1);
    
    /* monitor */
    if (monitor) {
      if (!(crun%10)) printf("   Run     #Term     #AllTerm\n");
      printf("%6i    %6i    %6i\n",crun,n1,nterm1);
      fflush(stdout);
    }
    
    /* print to file */
    fprintf(fp,"%d %d\n",nterm2,nterm1);
    fflush(fp);

  }
  
  *numterm = nterm1;

  return plterm1; 
}
