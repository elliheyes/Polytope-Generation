/*  ======================================================================  */
/*  ==========	     			                      	==========  */
/*  ==========       E V O L U T I O N   F U N C T I O N S      ==========  */
/*  ==========				                        ==========  */
/*  ======================================================================  */

#include "Global_5d_7p.h"

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
struct bitlist * termstates(struct population *evol, int numgen, int *numterm)
{
  int i,j, k;
  struct bitlist *bl;
  
  /* find the number of terminal states */
  *numterm=0;
  for (k=0; k<numgen; k++)
    for (i=0; i<evol[k].size; i++) 
      if ((evol[k].bl)[i].terminal){
        (*numterm)++;
      }
 
   /* allocate memory for terminal states */
  bl = calloc(*numterm,sizeof(struct bitlist)); 

  /* store terminal states in allocated space */
  j=0;
  for (k=0; k<numgen; k++)
    for (i=0; i<evol[k].size; i++) 
      if ((evol[k].bl)[i].terminal){
	    bl[j]=(evol[k].bl)[i]; j++;
      }
				 
   return bl;
}


/* remove equality and equivalence redundancy in list of bitlists */
void removeredundancy(struct bitlist *bl, int *len)
{
  /* remove equality redundancy */
  int cnonred, cactive, red, k; 
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (bitlistsequal(bl[k],bl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    bl[cnonred]=bl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(bl,cnonred,sizeof(struct bitlist),compbitlist);
  }
  
  /* remove equivalence redundancy */
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (bitlistsequiv(bl[k],bl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    bl[cnonred]=bl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(bl,cnonred,sizeof(struct bitlist),compbitlist);
  }
}

/* remove equality redundancy in list of bitlists */
void removeequality(struct bitlist *bl, int *len)
{
  int cnonred, cactive, red, k; 
  
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (bitlistsequal(bl[k],bl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    bl[cnonred]=bl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(bl,cnonred,sizeof(struct bitlist),compbitlist);
  }
  
  
}

/* remove equivalence redundancy in list of bitlists */
void removeequiv(struct bitlist *bl, int *len)
{
  int cnonred, cactive, red, k; 
  
  if (*len>1) {
    cnonred=1; cactive=1;
    while (cactive < *len) {
      red=0; k=0;
      while (!red && k<cnonred) {
        if (bitlistsequiv(bl[k],bl[cactive])) red=1;
        k++;
      }
      if (!red) {
	    bl[cnonred]=bl[cactive];
	    cnonred++;
	  }
      cactive++;
    }
    *len=cnonred;
    qsort(bl,cnonred,sizeof(struct bitlist),compbitlist);
  }
}


/* select terminal states from a population and remove redundancy */
struct bitlist * termstatesred(struct population *evol, int numgen, int *numterm)
{
  struct bitlist *bl;

  /* extract terminal states */
  bl=termstates(evol,numgen,numterm);
  
  /* remove redundancy in the list of terminal states */
  removeredundancy(bl,numterm);
  
  return bl;
}  


/* repeated evolution of a random initial population, extracting terminal states */
struct bitlist * searchenv(int numevol, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, FILE * fp, int *numterm)
{
  int cevol, n, i, nterm;
  struct population *evol;
  struct bitlist *bl, *blterm, *bltermOld;

  /* main loop over evolutions */
  nterm=0;
  for (cevol=0; cevol<numevol; cevol++){
    
    /* evolve random population */
    evol=evolvepop(randompop(popsize),numgen,meth,numcuts,keepfitest,mutrate,alpha,0);
    
    /* extract terminal states and remove redundancy */
    bl=termstatesred(evol,numgen,&n);
    nterm=nterm+n;
    
    /* allocate or re-allocate memory for terminal states */
    if (cevol==0) {
        blterm=calloc(nterm,sizeof(struct bitlist));
    }
    else {
        if(n!=0){
            bltermOld=calloc(nterm - n,sizeof(struct bitlist));;
            for (i=0; i< nterm - n; i++) bltermOld[i]=blterm[i];
            free(blterm);
            blterm=calloc(nterm,sizeof(struct bitlist));
            for (i=0; i< nterm - n; i++) blterm[i]=bltermOld[i];
            free(bltermOld);
        }
        /*blterm=realloc(blterm,(*nterm)*sizeof(struct bitlist));*/
    }
    
    if (blterm==NULL) {
      printf("serchenv: memory allocation failed");
      exit(0);
    }
    
    /* load in new terminal states */
    for (i=0; i<n; i++) blterm[(nterm)-n+i]=bl[i];
    
    /* remove redundancy on entire set */
    removeredundancy(blterm,&nterm);
    
    /* monitor */
    if (monitor) {
      if (!(cevol%10)) printf("   Run     #Term     #AllTerm\n");
      printf("%6i    %6i    %6i\n",cevol,n,nterm);
      fflush(stdout);
    }
    
    /* print to file */
    fprintf(fp,"%d %d\n",n,nterm);
    fflush(fp);
    
    /* free allocated memory */ 
    free(bl);free(evol);
  }
  
  *numterm = nterm;

  return blterm; 
}
