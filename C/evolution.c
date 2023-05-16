/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ==========       E V O L U T I O N   F U N C T I O N S      ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

#include "Global.h"


/* genetically evolve a population */
population * evolvepop(population initialpop, int numgen, int meth, int numcuts, int keepfitest, 
						float mutrate, float alpha, int monitor)
{
  int gen, nterm;
  population *evol;

   /* check validity of numgen */
  if (numgen>NUMGEN) numgen=NUMGEN;
  
  /* allocate memory for evolution */
  evol = calloc(numgen,sizeof(population));
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
void monitorevol(int gen, population *pop)
{
  if (!(gen%10)) printf("Gen    AvFit     MaxFit    #Term\n");
  printf("%3i    %2.4f    %2.4f    %3i\n",gen,pop->avfitness,pop->maxfitness,pop->nterm);
  fflush(stdout);
}  


/* select terminal states from a population */
bitlist * termstates(population *evol, int numgen, int *numterm)
{
  int i,j, k;
  bitlist *bl;
  
  /* find the number of terminal states */
  *numterm=0;
  for (k=0; k<numgen; k++)
    for (i=0; i<evol[k].size; i++) 
      if ((evol[k].bl)[i].terminal){
        (*numterm)++;
      }
 
   /* allocate memory for terminal states */
  bl = calloc(*numterm,sizeof(bitlist)); 

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
void removeredundancy(bitlist *bl, int *len)
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
    qsort(bl,cnonred,sizeof(bitlist),compbitlist);
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
    qsort(bl,cnonred,sizeof(bitlist),compbitlist);
  }
}


/* select terminal states from a population and remove redundancy */
bitlist * termstatesred(population *evol, int numgen, int *numterm)
{
  bitlist *bl;

  /* extract terminal states */
  bl=termstates(evol,numgen,numterm);
  
  /* remove redundancy in the list of terminal states */
  removeredundancy(bl,numterm);
  
  return bl;
} 


/* select new terminal states from a generated list */
void newtermstates(NormalForm * NFsOld, bitlist * blNew, int numtermOld, int len, int *numtermNew)
{
  int cnew, cactive, old, k;
  NormalForm NF;
  
  cnew=0; cactive=0;
  while (cactive < len){
    /* compute the normal form of the new polytope */
    NF = normalform(blNew[cactive]);
  
    /* check whether normal form is in existing list */
    old=0; k=0; 
    while (!old && k<numtermOld){ 
      if (NFsequal(NFsOld[k],NF)) old=1;
      k++;
    }
    
    /* if normal form is not in existing list then bring it to the front of the list */
    if (!old) {
	  blNew[cnew]=blNew[cactive];
	  
	  /* update index */
	  cnew++;
	}	
	/* update index */
    cactive++;
  }
  
  *numtermNew=cnew;
  qsort(blNew,cnew,sizeof(bitlist),compbitlist);
} 


/* repeated evolution of a random initial population, extracting terminal states */
bitlist * searchenv(int numrun, int numevol, int numgen, int popsize, int meth, int numcuts,
			   int keepfitest, float mutrate, float alpha, int monitor, int *numterm)
{
  int i, j, k, IP, crun, cevol, n1, n2, nterm1, nterm2;
  population *evol;
  bitlist *bl, *blterm1, *blterm2, *bltermOld1, *bltermOld2;
  NormalForm *NFs, *NFsOld;
  
  FILE * fp1 = fopen("Num_Terminal_States.txt","w");
  FILE * fp2 = fopen("Terminal_States.txt","w");

  /* main loop over runs */
  nterm1=0;
  for(crun=0; crun<numrun; crun++){
  
    /* loop over evolutions */
    nterm2=0;
    for(cevol=0; cevol<numevol; cevol++){
    
      /* evolve random population */
      evol=evolvepop(randompop(popsize),numgen,meth,numcuts,keepfitest,mutrate,alpha,0);
      
      /* extract terminal states and remove redundancy */
      bl=termstatesred(evol,numgen,&n2);
      nterm2=nterm2+n2;
      
      /* allocate or re-allocate memory for terminal states */
      if(cevol==0) {
        blterm2=calloc(nterm2,sizeof(bitlist));
      }
      else{
        if(n2!=0){
            bltermOld2=calloc(nterm2-n2,sizeof(bitlist));;
            for (i=0; i< nterm2-n2; i++) bltermOld2[i]=blterm2[i];
            free(blterm2);
            blterm2=calloc(nterm2,sizeof(bitlist));
            for (i=0; i< nterm2-n2; i++) blterm2[i]=bltermOld2[i];
            free(bltermOld2);
        }
      }
    
      /* load in new terminal states */
      for (i=0; i<n2; i++) blterm2[nterm2-n2+i]=bl[i];

	  /* free allocated memory */ 
      free(bl);free(evol);
    }
    
    if(crun==0){
      /* remove redundancy */
      removeredundancy(blterm2, &nterm2);
      n1=nterm2;
      nterm1=n1;
      
      /* allocate memory for terminal states */
      blterm1=calloc(nterm1,sizeof(bitlist));
      
      /* allocate memory for normal forms */
      NFs=calloc(nterm1,sizeof(NormalForm));
    }
    else{
      /* remove redundancy */
      removeredundancy(blterm2, &nterm2);
      
      /* select new terminal states */
      newtermstates(NFs,blterm2,nterm1,nterm2,&n1);
      nterm1=nterm1+n1;
      
      /* re-allocate memory for terminal states and normal forms */
      if(n1!=0){
        bltermOld1=calloc(nterm1-n1,sizeof(bitlist));
        NFsOld=calloc(nterm1-n1,sizeof(NormalForm));
        for (i=0; i< nterm1-n1; i++){
          bltermOld1[i]=blterm1[i];
          NFsOld[i]=NFs[i];
        };
        free(blterm1);free(NFs);
        blterm1=calloc(nterm1,sizeof(bitlist));
        NFs=calloc(nterm1,sizeof(NormalForm));
        for (i=0; i< nterm1-n1; i++){
          blterm1[i]=bltermOld1[i];
          NFs[i]=NFsOld[i];
        };
        free(bltermOld1);free(NFsOld);
      }
    }
    
    /* load in new terminal states and normal forms */
    for (i=0; i<n1; i++){
      blterm1[(nterm1)-n1+i] = blterm2[i];
      NFs[(nterm1)-n1+i] = normalform(blterm2[i]);
    }; 
    
    /* monitor */
    if (monitor) {
      if (!(crun%10)) printf("   Run     #Term     #AllTerm\n");
      printf("%6i    %6i    %6i\n",crun,nterm2,nterm1);
      fflush(stdout);
    }
    
    /* print number of terminal states to file */
    fprintf(fp1,"%d %d\n",nterm2,nterm1);
    fflush(fp1);
    
    /* print terminal states to file */
    for (i=0; i<n1; i++) fprintbitlist(fp2, blterm2[i]); 
  }
  
  /* close files */
  fclose(fp1); fclose(fp2); 
  
  /* assign total number of reduced terminal states */
  *numterm = nterm1;

  /* return list of reduced terminal states */
  return blterm1; 
}

