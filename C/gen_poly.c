/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========        P O L Y T O P E   G E N E T I C           ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"

int main (int narg, char* fn[]){

  int n;
  
  /* repeated evolution of a random initial population, extracting terminal states */
  struct pointlist * pl = searchenv(1000, 1, NUMGEN, POPSIZE, ROULETTE, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA, MONITORON, &n); 

  return 0;
} 
