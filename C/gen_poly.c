/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========        P O L Y T O P E   G E N E T I C           ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"

int main (int narg, char* fn[]){

  int n;
  
  /* set random number generator seed */
  srand(clock()); 
  
  /* repeated evolution of a random initial population, extracting terminal states */
  bitlist * bl = searchenv(1000, NUMGEN, POPSIZE, ROULETTE, NUMCUTS, KEEPFITEST, MUTRATE, ALPHA, MONITOR, &n); 

  return 0;
} 
