/*  ======================================================================  */
/*  ==========	     			   	  	        ==========  */
/*  ==========         B I T L I S T   F U N C T I O N S        ==========  */
/*  ==========						        ==========  */
/*  ======================================================================  */

#include "Global.h"

/* compare two arrays, return 0 if identical and 1 otherwise */
char compareArray(int a[],int b[],int size){
	int i;
	for(i=0;i<size;i++) if(a[i]!=b[i]) return 1;
	return 0;
}


/* convert decimal number into binary number */
binary decimal2binary(int num)
{
  binary bin;
  int binlist[BINLEN];
  int i,j,k;
  
  if(num == 0){
    for(i=0; i < BINLEN; i++) bin.list[i] = 0;
    return bin;
  }
  else{
    for(j=0; num > 0; j++){
      binlist[j] = num % 2;
      num /= 2;
    }
  
    for(k=0; k < BINLEN-j; k++) bin.list[k] = 0;
    
    int m,l=BINLEN-j;
    for(m=j-1; m >= 0; m--) bin.list[l++] = binlist[m];
  
    return bin;
  }
}


/* convert bitlist to an integer */
int binaryToDecimal(binary bin)
{
  int i,num=0;
  for (i=0; i<BINLEN; i++) num=num+bin.list[i]*pow(2,BINLEN-i-1);
  return num;
}  


/* convert a pointlist to a bitlist after adding max */
bitlist pts2bts(pointlist pl)
{
  int i,j,k,num;
  bitlist bl, *blp=&bl;;
  binary zerobin, bin;
  
  /* define bitlist length */
  bl.len = MAXNVRTS*POLYDIM*BINLEN;
  
  /* define binary number corresponding to 0 */
  zerobin = decimal2binary(-MIN);

  /* convert decimal points to binary numbers and add to bitlist */
  for(i=0; i<MAXNVRTS; i++){
    if (i<pl.len) {
      for(j=0; j < POLYDIM; j++){
        num = pl.points[i][j];
        bin = decimal2binary(num - MIN);
        for(k=0; k < BINLEN; k++) bl.bits[(i*POLYDIM*BINLEN)+(j*BINLEN)+k] = bin.list[k];
      }
    }
    else {
      for(j=0; j < POLYDIM; j++) 
        for(k=0; k < BINLEN; k++) bl.bits[(i*POLYDIM*BINLEN)+(j*BINLEN)+k] = zerobin.list[k];
    }
  }
  
  /* compute the fitness */
  fitness(blp);
  
  return bl;
}


/* convert a bitlist into a reduced pointlist and subtract max */
pointlist bts2pts(bitlist bl)
{
  int i,j,k,count=0,new;
  int point[POLYDIM];
  pointlist pl;
  binary bin;
  
  for(i=0; i<MAXNVRTS; i++){
    for(j=0; j<POLYDIM; j++){
      for(k=0; k<BINLEN; k++){
        bin.list[k] = bl.bits[(i*POLYDIM*BINLEN)+(j*BINLEN)+k];
      }
      point[j] = binaryToDecimal(bin) + MIN;  
    }
    
    /* check if new point already exists in the list */
    new=1;
    for(j=0; j<count; j++) if(!compareArray(point, pl.points[j], POLYDIM)) new=0;
    
    if(new){
      for(j=0; j<POLYDIM; j++) pl.points[count][j] = point[j]; 
      count++;
    }
    
  }
  pl.len = count;

  return pl;
}


/* generate a random integer in a range from min to max */
int randomint(int min, int max)
{
  return round((1.* rand())/RAND_MAX*(max-min))+min;
}  


/* generate a random bitlist state */
bitlist randomstate()
{
  int i,j;
  bitlist bl;
  pointlist pl;
  
  /* randomly generate points */
  for(i=0; i<MAXNVRTS; i++) for(j=0; j<POLYDIM; j++) pl.points[i][j] = randomint(MIN,MIN-1+pow(2,BINLEN));
  pl.len = MAXNVRTS;
  
  /* convert pointlist to bitlist */
  bl = pts2bts(pl);
  
  return bl;
} 


/* generate a random choice for integers 0,...,len-1 for a probability distribution p */
int randomchoice(float p[POPSIZE], int len)
{
  float psum,P,fran;
  int i;

  /* determine normalisation */
  psum=0;
  for (i=0; i<len;i++ ) psum=psum+p[i];

  /* generate random float in range [0,1] */
  fran=(1.*rand())/RAND_MAX;

  P=0.; i=0;
  while (P<fran && i<len) {P=P+p[i]/psum; i++;}

  return i-1;
}  

 
/* flips bit in position pos for a bitlist bl */
void flipbit(bitlist *bl, int pos)
{
  if ((pos>=0) && (pos<(bl->len))) (bl->bits)[pos]=((bl->bits)[pos]+1) % 2;
}  


/* copies the bits in bitlist bl from position pos1 to pos2-1 into a new bitlist */
bitlist copybitlist(bitlist bl, int pos1, int pos2)
{
  bitlist blcopy;
  int i,p1,p2;

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
void pastebitlist(bitlist *bl, bitlist blpaste, int pos)
{
  int len,i;

  /* make sure position and length are ok */
  if (pos<0) pos=0;
  if (pos>=bl->len) pos=bl->len-1;
  len=((pos+blpaste.len<=bl->len) ? blpaste.len : bl->len-pos);

  /*copy b1paste into bl */
  if (len>0) for (i=0; i<len; i++) (bl->bits)[pos+i]=blpaste.bits[i];
}

									  
/* a simple insertion sort of an integer array into ascending order */
void isort(int arr[], int len)
{
    int i,key,j;
    for (i=1; i<len; i++)
    {
        key=arr[i];
        j=i-1;
        while (j>=0 && arr[j]>key)
        {
            arr[j+1]=arr[j];
            j=j-1;
        }
        arr[j+1]=key;
    }
}
									  

/* crosses bitlists bl1 and bl2, with numcuts number of cuts at positions specified in array cuts */
void crossbitlists(bitlist *bl1, bitlist *bl2, int numcuts, int cuts[NUMCUTS])
{

  int i,minlen,start,end,cuts1[NUMCUTS+2];
  bitlist blpart1, blpart2;

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
  const bitlist *bl1=p1, *bl2=p2;

  if (bl1->fitness < bl2->fitness) return 1;
  else if  (bl1->fitness > bl2->fitness) return -1;
  else return 0;
}  

/* compute the normal form of a polytope from its bitlist */
NormalForm normalform(bitlist bl)
{
  int i,j,IP;
  pointlist pl;
  PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList));
  VertexNumList V;
  EqList *E = (EqList *) malloc(sizeof(EqList));
  VPermList *VP = (VPermList*) malloc(sizeof(VPermList));
  Long NF[POLYDIM][VERT_Nmax];
  NormalForm NF0;
  
  /* transform the bitlists into point lists */
  pl = bts2pts(bl);
  
  /* define the polytope dimension and the number of points */
  _P->n=POLYDIM;
  _P->np=pl.len; 
  
  /* fill in the points */
  for(i=0; i<pl.len; i++) for(j=0; j<POLYDIM; j++) _P->x[i][j]=pl.points[i][j];

  /* find the bounding hyperplane equations of the polytope */
  IP=Find_Equations(_P,&V,E); 

  /* compute normal form */
  Make_Poly_Sym_NF(_P, &V, E, NF);
  NF0.nv=V.nv; 
  for (i=0; i<POLYDIM; i++) for(j=0; j<VERT_Nmax; j++) NF0.x[i][j]=NF[i][j];
  
  free(_P);free(E);free(VP);
  
  return NF0;
}


/* decide if two bistlist are identical */
int bitlistsequal(bitlist bl1, bitlist bl2)
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


/* decide if two normal forms are equal */
int NFsequal(NormalForm NF1, NormalForm NF2)
{
  int i, j, equal;
  
  /* check if the number of vertices match */
  if (NF1.nv != NF2.nv) return 0;
  
  /* if number of vertices match, compare the normal form matrices */
  else{
    equal=1; 
    for(i=0; i<POLYDIM; i++){
      for(j=0; j<NF1.nv; j++){
        if (NF1.x[i][j]!=NF2.x[i][j]){
          equal = 0;
          break;
        } 
      }
      if (!equal) break;
    }
    return equal;
  }
}


/* decide if two bistlist are equivalent by comparing their normal forms */
int bitlistsequiv(bitlist bl1, bitlist bl2)
{
  int i, j, equal;
  pointlist pl1, pl2;
  VertexNumList V01, V02;
  EqList *E01 = (EqList *) malloc(sizeof(EqList)),
         *E02 = (EqList *) malloc(sizeof(EqList));
  PolyPointList *_P01 = (PolyPointList *) malloc(sizeof(PolyPointList)),
   				*_P02 = (PolyPointList *) malloc(sizeof(PolyPointList));
  int IP1, IP2;
  
  /* transform the bitlists into point lists */
  pl1 = bts2pts(bl1);
  pl2 = bts2pts(bl2);
  
  /* define the polytope dimension and number of points */
  _P01->n=POLYDIM; _P02->n=POLYDIM; 
  _P01->np=pl1.len; _P02->np=pl2.len;
    
  /* define the points */
  for(i=0; i<pl1.len; i++) for(j=0; j<POLYDIM; j++) _P01->x[i][j]=pl1.points[i][j];
  for(i=0; i<pl2.len; i++) for(j=0; j<POLYDIM; j++) _P02->x[i][j]=pl2.points[i][j];
    
  /* find the bounding hyperplane equations of the polytopes */
  IP1 = Find_Equations(_P01,&V01,E01); 
  IP2 = Find_Equations(_P02,&V02,E02); 
  
  /* check if the number of vertices match */
  if (V01.nv != V02.nv){
     /* free allocated memory */
  	 free(E01);free(E02);free(_P01);free(_P02);
  
     /* if the number of vertices don't match then return 0 */
  	 return 0;
  }
  else {
  	VPermList *VP01 = (VPermList*) malloc(sizeof(VPermList)),
              *VP02 = (VPermList*) malloc(sizeof(VPermList)); 
  	Long NF1[POLYDIM][VERT_Nmax], NF2[POLYDIM][VERT_Nmax];
  	NormalForm NF01, NF02;
  
    /* compute the normal forms of the polytopes */
	Make_Poly_Sym_NF(_P01, &V01, E01, NF1);
    Make_Poly_Sym_NF(_P02, &V02, E02, NF2);   
    NF01.nv = V01.nv; NF02.nv = V02.nv; 
    for (i=0; i<POLYDIM; i++){
      for(j=0; j<VERT_Nmax; j++){
         NF01.x[i][j]=NF1[i][j];
         NF02.x[i][j]=NF2[i][j];
      }
    }
        
   	/* compare the normal forms of the two polytopes */
    equal = NFsequal(NF01, NF02); 
    
    /* free allocated memory */
    free(E01);free(E02);free(_P01);free(_P02);free(VP01);free(VP02); 
    
    return equal;
  }
  
}


/* write bitlist to a file in the format of the list of vertices  */
void fprintbitlist(FILE * fp, bitlist bl)
{
  int i,j,IP;
  pointlist pl;
  VertexNumList V;
  EqList *E = (EqList *) malloc(sizeof(EqList));
  PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList));
  
  pl = bts2pts(bl);
  
  _P->n=POLYDIM; 
  _P->np=pl.len; 
  
  for(i=0; i<pl.len; i++) for(j=0; j<POLYDIM; j++) _P->x[i][j]=pl.points[i][j];
    
  IP=Find_Equations(_P,&V,E); 
  
  fprintf(fp,"[");
  for(i=0; i<V.nv; i++){
    fprintf(fp,"[");
    for(j=0; j<POLYDIM; j++){
      fprintf(fp,"%lld",_P->x[V.v[i]][j]);
      if(j!=POLYDIM-1) fprintf(fp,",");
    }
    fprintf(fp,"]");
  }
  fprintf(fp,"]\n");
  
  free(E);free(_P);
}


/* read bitlists from a file */
bitlist * freadbitlist(char *filename, int len)
{
  char * line;
  int i, j, k, c, d;
  size_t size = 0;
  int arr[100];
  pointlist pl;
  bitlist * bl = calloc(len,sizeof(bitlist)); 
  
  FILE * fp = fopen(filename, "r");
  
  if (fp == NULL) {
    printf("Error opening the file\n");
    exit(1);
  }
  
  /* read file in line by line */
  line = NULL;
  i=0; 
  while (getline(&line, &size, fp ) != -1) {  
    /* determine the number of vertices */
    j=0; c=0; d=0;
    while (j<strlen(line)) {
      if (line[j] == '[' || line[j] == ',' || line[j]=='\n') j++;
      else if (line[j] == ']') {
        j++; c++;
      }
      else if (line[j] == '-') {
        j++; j++; d++;
      }
      else {
        j++; d++;
      }
    }
    
    /* write vertices to a pointslist */
    pl.len = c-1;
    j=0; c=0; d=0;
    while (j<strlen(line)) {
      
      if (line[j] == '[' || line[j] == ',' || line[j]=='\n') j++;
      else if (line[j] == ']') {
        j++; c++; d=0;
      }
      else if (line[j] == '-') {
        pl.points[c][d] = -((int)(line[j+1])-(int)('0'));
        j++; j++; d++;
      }
      else {
        pl.points[c][d] = (int)(line[j])-(int)('0');
        j++; d++;
      }
    }

    /* convert to bitlist and add to list */
    bl[i] = pts2bts(pl);
    
    /* update index */
    i++;
  }
  
  /* close file */ 
  fclose(fp);

  /* return list of bitlists */  
  return bl;
}
