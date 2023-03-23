/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========        P O L Y T O P E   G E N E T I C           ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"

int main (int narg, char* fn[]){

  int i, n;
  
  /* repeated evolution of a random initial population, extracting terminal states */
  FILE * fp1 = fopen("4d_5v_Num_Terminal_States.txt","w");
  struct pointlist * pl = searchenv(10, 1, NUMGEN, POPSIZE, ROULETTE, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA, MONITORON, fp1, &n); 
  fclose(fp1); 

  /* write terminal states to a file */
  FILE * fp2 = fopen("4d_5v_Terminal_States.txt","w");  
  for(i=0; i<n; i++) fprintpointlist(fp2, pl[i]); 
  fclose(fp2); 

  return 0;
} 
