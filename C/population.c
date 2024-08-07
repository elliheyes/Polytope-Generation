/*  ======================================================================  */
/*  ==========	     			   	  	    	==========  */
/*  ==========      P O P U L A T I O N   F U N C T I O N S     ==========  */
/*  ==========						        ==========  */
/*  ======================================================================  */

#include "Global.h"


/* find fittest in a population */
bitlist fitestinpop(population *pop)
{
  float fmax;
  int pos, i;

  pos=0;
  if (pop->size>0) {
    fmax=((pop->bl)[0]).fitness; 
    for (i=1; i<pop->size; i++)
      if (((pop->bl)[i]).fitness>fmax) {pos=i; fmax=((pop->bl)[i]).fitness;}
  }

  return (pop->bl)[pos];
}  
  
      
/* find average fitness of a population */
float avfitness(population *pop)
{
  float fsum=0.;
  int i;

  for (i=0; i<pop->size; i++) fsum=fsum+((pop->bl)[i]).fitness;

  return fsum/pop->size;
} 


/* find number of terminal state in population */
int nterminalpop(population *pop)
{
  int i,nterm;

  nterm=0;
  for (i=0; i<pop->size; i++) if (((pop->bl)[i]).terminal) nterm++;

  return nterm;
}   


/* generate a random population of size popsize */
population randompop(int popsize)
{
  int size, i;
  population pop;
  bitlist blfitest;

  size=((popsize>POPSIZE) ? POPSIZE : popsize);
  for (i=0; i<size; i++) pop.bl[i]=randomstate();
  pop.size=size;

  /* compute average and maximal fitness for random population */
  pop.avfitness=avfitness(&pop);    /* average fitness in population */
  blfitest=fitestinpop(&pop);		/* fittest in population */
  pop.maxfitness=blfitest.fitness;  /* maximal fitness in population */
  pop.nterm=nterminalpop(&pop);     /* number of terminal state in mutated population */
  
  return pop;
}  


/* sort a population by fitness */
void sortpop(population *pop)
{
  qsort(pop->bl,pop->size,sizeof(bitlist),compbitlist);
}


/* mutate a population */
void mutatepop(population *pop, float mutrate)
{
  int len, nbits, nmut, i, k, pos, ipos, bpos;
  bitlist blfitest;

  len=((pop->bl)[0]).len;     /* length of bitlist, assume they are the same across population */
  nbits=len*(pop->size);      /* total number of bits in population */
  nmut=round(nbits*mutrate);  /* number of bits to be mutated */

  /* run over mutations */
  for (i=0; i<nmut; i++) {
    pos=randomint(0,nbits);			    /* random bit position in entire population */
    ipos=pos/len;       	     	    /* individual where mutation occurs */
    bpos=pos%len;          			    /* bit position in individual */
    flipbit(&((pop->bl)[ipos]),bpos);	/* flip bit */
  }
  
  for (k=0; k<pop->size; k++) fitness(&((pop->bl)[k])); /* update fitness for mutated population */
  pop->avfitness=avfitness(pop);                        /* average fitness in mutated population */
  blfitest=fitestinpop(pop);		                     /* fittest in population */
  pop->maxfitness=blfitest.fitness;                     /* maximal fitness in mutated population */
  pop->nterm=nterminalpop(pop);                         /* number of terminal state in mutated population */
  
}  


/* generate the next population from a given one by selection, crossing and mutation */
void nextpop(population pop, population *newpop, int meth, int numcuts, int keepfitest,
	     float mutrate, float alpha)
{
  bitlist blfitest, bl1, bl2;
  float fmax, fav, df, p[POPSIZE];
  int n, len, k, i, i1, i2, j, cutpos[NUMCUTS];

   n=pop.size;             /* size of population */
   len=((pop.bl)[0]).len;  /* length of bitlist - assume the same across population */
  
   if (n>1){

     /* find fitest in population */
     blfitest=fitestinpop(&pop);
  
    /* fix parameters if out of range */
    if (numcuts<1 || numcuts>NUMCUTS) numcuts=1;
    if (alpha<=1.) alpha=2.;
 
    /* compute selection probabilities p */
    
    /* ranking method */
    if (meth==0) {
      sortpop(&pop);
      for (k=0; k<n; k++) p[k]=2*(1+(n-k-1)/(n-1)*(alpha-1))/(1+alpha)/n;
    }  
    else {
      /* roulette method */
      fmax=blfitest.fitness;
      fav=avfitness(&pop);
      df=fmax-fav;
      if (df<=0) for (k=0; k<n; k++) p[k]=1./n;
      else for (k=0; k<n; k++) p[k]=((alpha-1)*(((pop.bl)[k]).fitness-fav)+df)/df/n;
    }

    /* select and cross pairs of individuals */
    i=0;
    while (i<n){
	  /* select the positions of two individuals */
      i1=randomchoice(p,n); i2=randomchoice(p,n);             
      
      /* copy the two selected individuals */
      bl1=pop.bl[i1]; bl2=pop.bl[i2];                       
      
      /* determine the random cut positions */
      for (j=0; j<numcuts; j++) cutpos[j]=randomint(0,len-1); 
      
      /* cross the two individuals */
      crossbitlists(&bl1,&bl2,numcuts,cutpos);    
      
      /* copy the crossed individuals into the new population */
      (newpop->bl)[i]=bl1; (newpop->bl)[i+1]=bl2;       
           
      i=i+2;
    }
    
    newpop->size=n;

    /* mutate new population */
    mutatepop(newpop,mutrate);

    /* if fitest individual is kept copy to new population */
    if (keepfitest) (newpop->bl)[0]=blfitest;
 }  
} 
