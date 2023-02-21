/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========        P O L Y T O P E   G E N E T I C           ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"

int main (int narg, char* fn[]){

  /* initialise the population */
  struct population initialpop = randompop(POPSIZE);
  
  /* evolve the population over generations */
  struct population * evol = evolvepop(initialpop, NUMGEN, ROULETTE, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA, MONITORON); 
  
  /* extract reduced list of terminal states */
  int n; 
  struct bitlist *bl;
  bl = termstatesred(evol, NUMGEN, &n);
  printf("Number of reduced terminal states: %d\n",n);

  /* write terminal states to a file */
  FILE * fp = fopen("Terminal_States.txt","w");  
  for(int i=0; i<n; i++) fprintbitlist(fp, bl[i]); 

  return 0;
} 




