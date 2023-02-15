/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========        P O L Y T O P E   G E N E T I C           ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"

int main (int narg, char* fn[]){

  struct population initialpop = randompop(POPSIZE);
  
  struct population * pop = evolvepop(initialpop, NUMGEN, ROULETTE, NUMCUTS, KEEPFITEST, STDMUTRATE, STDALPHA, MONITORON); 
  
  /* struct bitlist * blterm = searchenv(10, NUMGEN, POPSIZE, ROULETTE, NUMCUTS, KEEPFITEST, STDMUTRATE, STDALPHA, MONITORON); */
 
  return 0;
} 
