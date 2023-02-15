/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ==========         B I T L I S T   F U N C T I O N S        ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

#include <math.h>
#include "Global.h"


/* convert decimal into binary */
struct binary decimal2binary(int num)
{
  struct binary bin;
  int binlist[BINLEN];
  int i,j,k;
  
  if(num == 0){
    for(i=0; i < BINLEN; i++){
      bin.list[i] = 0;
    }
    return bin;
  }
  
  for(j=0; num > 0; j++){
    binlist[j] = num % 2;
    num /= 2;
  }
  
  for(k=0; k < BINLEN-j; k++){
    bin.list[k] = 0;
  }
    
  int l=BINLEN-j;
  int m;
  for(m=j-1; m >= 0; m--){
    bin.list[l++] = binlist[m];
  }
  
  return bin;
}


/* convert bit list to integer */
int binaryToDecimal(struct binary bin)
{
  int i, num;

  num=0;
  for (i=0; i<BINLEN; i++)
    num=num+bin.list[i]*pow(2,BINLEN-i-1);

  return num;
}  


/* convert a points list into a bit list after adding max */
struct bitlist pts2bts(struct pointlist pl)
{
  struct bitlist bl;
  struct bitlist * blp = &bl;
  int i,j,k;
  int num;
  struct binary bin;
    
  bl.len = NPTS*POLYDIM*BINLEN;

  for(i=0; i < NPTS; i++){
    for(j=0; j < POLYDIM; j++){
      num = pl.points[i][j];
      bin = decimal2binary(num - MIN);
      for(k=0; k < BINLEN; k++){
        bl.bits[(i*POLYDIM*BINLEN)+(j*BINLEN)+k] = bin.list[k];
      }
    }
  }
  
  fitness(blp);
  if(bl.fitness == 1) bl.terminal = 1;
  else bl.terminal = 0;
    
  return bl;
}


/* convert a bit list into a points list and subtract max */
struct pointlist bts2pts(struct bitlist bl)
{
  struct pointlist pl;
  int i,j,k;
  struct binary bin;
  
    for(i=0; i < NPTS; i++){
      for(j=0; j < POLYDIM; j++){
        for(k=0; k < BINLEN; k++){
          bin.list[k] = bl.bits[(i*POLYDIM*BINLEN)+(j*BINLEN)+k];
        }
        pl.points[i][j] = binaryToDecimal(bin) + MIN;
      }
    }
  
  return pl;
}


/* generate a random integer in a range from min to max */
int randomint(int min, int max)
{
  return round((1.* rand())/RAND_MAX*(max-min))+min;
}  


/* generate a random state */
struct bitlist randomstate()
{
  struct bitlist bl;
  struct pointlist pl;
  int i,j;
  
  for(i=0; i<NPTS; i++){
  	for(j=0; j<POLYDIM; j++){
  	  pl.points[i][j] = randomint(MIN,MIN-1+pow(2,BINLEN));
  	}
  }
  
  bl = pts2bts(pl);
  
  return bl;
} 


/* generate a random choice for integers 0,...,len-1 for a probability distribution p */
int randomchoice(float p[POPSIZE], int len)
{
  float psum, P, fran;
  int i;

  /* determine normalisation */
  psum=0;
  for (i=0; i<len;i++ ) psum=psum+p[i];

  /* generate random float in range [0,1] */
  fran=(1.*rand())/RAND_MAX;

  P=0.; i=0;
  while (P<fran) {P=P+p[i]/psum; i++;}

  return i-1;
}  

 
/* flips bit in position pos for a bitlist bl */
void flipbit(struct bitlist *bl, int pos)
{
  if ((pos>=0) && (pos<(bl->len))){
    (bl->bits)[pos]=((bl->bits)[pos]+1) % 2;
  } 
}  


/* copies the bits in bitlist bl from position pos1 to pos2-1 into a new bitlist */
struct bitlist copybitlist(struct bitlist bl, int pos1, int pos2)
{
  struct bitlist blcopy;
  int i, p1, p2;

  /* make sure positions are viable */
  if (pos1<=pos2) {p1=pos1; p2=pos2;}
  else {p1=pos2; p2=pos1;}
  if (p1<0) p1=0;
  if (p2<0) p2=0;
  if (p1>bl.len) p1=bl.len;
  if (p2>bl.len) p2=bl.len;

  /* copy section of bit list */
  if (p1<p2) {
    for (i=0; i<(p2-p1); i++) blcopy.bits[i]=bl.bits[p1+i];
    blcopy.len=p2-p1;
  }
  else blcopy.len=0;
  blcopy.fitness=0.; blcopy.terminal=0;

  return blcopy;
}  


/* pastes bitlist blpaste into bitlist bl at position pos, overwriting previous content in bl */
void pastebitlist(struct bitlist *bl, struct bitlist blpaste, int pos)
{
  int len, i;

  /* make sure position and length are ok */
  if (pos<0) pos=0;
  if (pos>=bl->len) pos=bl->len-1;
  len=((pos+blpaste.len<=bl->len) ? blpaste.len : bl->len-pos);

  /*copy b1paste into bl */
  if (len>0) {
    for (i=0; i<len; i++) (bl->bits)[pos+i]=blpaste.bits[i];
  }
}

									  
/* a simple insertion sort of an integer array into ascending order */
void isort(int arr[], int len)
{
    int i, key, j;
    for (i = 1; i < len; i++)
    {
        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}
									  

/* crosses bitlists bl1 and bl2, with numcuts number of cuts at positions specified in array cuts */
void crossbitlists(struct bitlist *bl1, struct bitlist *bl2, int numcuts, int cuts[NUMCUTS])
{

  int i, minlen, start, end, cuts1[NUMCUTS+2];
  struct bitlist blpart1, blpart2;

  /* set lengths of bitlists to common minimum */
  minlen=((bl1->len < bl2->len) ? bl1->len : bl2->len);
  bl1->len=minlen; bl2->len=minlen;

  /* check range of cut positions */
  cuts1[0]=0; cuts1[numcuts+1]=minlen;
  for (i=0; i<numcuts; i++) {
    cuts1[i+1]=cuts[i];
    if (cuts[i]<0) cuts1[i+1]=0;
    if (cuts[i]>=minlen) cuts1[i+1]=minlen-1;
  }  
  
  /* sort cut positions in ascending order */
  isort(cuts1,numcuts+2);
  /* for (i=0; i<numcuts+2; i++) printf("%i,",cuts1[i]);
     printf("\n"); */

  /* swap the relevant parts of the two bitlists */
  for (i=0; i<=numcuts; i=i+2) {
    start=cuts1[i]; end=cuts1[i+1];
    blpart1=copybitlist(*bl1,start,end);
    blpart2=copybitlist(*bl2,start,end);
    pastebitlist(bl1,blpart2,start);
    pastebitlist(bl2,blpart1,start);
  }  
}  

/* compare two bistlist by their fitness */
int compbitlist(const void *p1, const void *p2)
{
  const struct bitlist *bl1=p1, *bl2=p2;

  if (bl1->fitness < bl2->fitness) return 1;
  else if  (bl1->fitness > bl2->fitness) return -1;
  else return 0;
}  


/* decide if two bistlist are identical */
int bitlistsequal(struct bitlist bl1, struct bitlist bl2)
{
  int i, equal;
  
  if (bl1.len != bl2.len) return 0;
  else {
    equal=1; i=0;
    while (equal && i<bl1.len) {
      equal = equal && (bl1.bits[i]==bl2.bits[i]);
      i++;
    }

    return equal;
  }
}


