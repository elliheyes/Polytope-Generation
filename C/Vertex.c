#include "Global.h"
#include "Rat.h"


#ifndef	CEQ_Nmax	
#define CEQ_Nmax        EQUA_Nmax
#endif

typedef struct {int ne; Equation e[CEQ_Nmax];}              CEqList;

/*  ======================================================================  */
/*  ==========		     			  	                     	==========  */
/*  ==========	   I N C I D E N C E S (as bit patterns)    	==========  */
/*  ==========						                         	==========  */ 
/*  ======================================================================  */

#if (VERT_Nmax > LONG_LONG_Nbits)
INCI INCI_AND(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])&(y.ui[i]); return z;}
INCI INCI_OR(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])|(y.ui[i]); return z;}
INCI INCI_XOR(INCI x, INCI y){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=(x.ui[i])^(y.ui[i]); return z;}
int  INCI_EQ(INCI x, INCI y){
  int i; for (i=0;i<I_NUI;i++) if ((x.ui[i])!=(y.ui[i])) return 0; return 1;}
int  INCI_EQ_0(INCI x){
  int i; for (i=0;i<I_NUI;i++) if (x.ui[i]) return 0; return 1;}
int  INCI_LE(INCI x, INCI y){
  int i;
  /* for (i=0;i<I_NUI;i++) if ((x.ui[i]|y.ui[i])!=(y.ui[i])) return 0; */
  unsigned int *X=x.ui, *Y=y.ui;
  for (i=0;i<I_NUI;i++) if ((X[i]|Y[i])!=(Y[i])) return 0;
  return 1;}
INCI INCI_0(){
  INCI z; int i; for (i=0;i<I_NUI;i++) z.ui[i]=0; return z;}
INCI INCI_1(){
  INCI z; int i; z.ui[0]=1; for (i=1;i<I_NUI;i++) z.ui[i]=0; return z;}
INCI INCI_PN(INCI x, Long y){
  INCI z; int i; 
  z.ui[0]=(x.ui[0]<<1)|(!y); 
  for (i=1;i<I_NUI;i++) z.ui[i]=(x.ui[i]<<1)|(x.ui[i-1]>>(INT_Nbits-1)); 
  return z;}
INCI INCI_D2(INCI x){
  INCI z; int i; 
  for (i=0;i<I_NUI-1;i++) z.ui[i]=(x.ui[i]>>1)|(x.ui[i+1]<<(INT_Nbits-1));
  z.ui[I_NUI-1]=x.ui[I_NUI-1]>>1;
  return z;}
int  INCI_lex_GT(INCI *x, INCI *y){
  int i=I_NUI; while(i--) if(x->ui[i]>y->ui[i]) return 1; 
  else if(x->ui[i]<y->ui[i]) return 0; return 0; }
int  INCI_LmR(INCI *x, INCI *y){ puts("Implement INCI_LmR"); exit(0); }
#else
int  INCI_lex_GT(INCI *x, INCI *y){ return (*x > *y) ? 1 : 0 ; }
int  INCI_LmR(INCI *x, INCI *y){ return (*x>*y) ? 1 : (*x<*y) ? -1 : 0; }
/* int  Lead_Vert(INCI x){int i=0; while(!(x%2)) {i++; x/=2;} return i;} */
#endif

int INCI_abs(INCI X){
  int abs=0; while(!INCI_EQ_0(X)) {abs+=INCI_M2(X); X=INCI_D2(X);} return abs;
}

INCI Eq_To_INCI(Equation *_Eq, PolyPointList *_P, VertexNumList *_V){
  int j;  INCI X=INCI_0();
  for (j=0;j<_V->nv;j++) X=INCI_PN(X,Eval_Eq_on_V(_Eq,_P->x[_V->v[j]],_P->n));
  return X;
}

/*  ======================================================================  */
/*  ==========	     			   	  	                    	==========  */
/*  ========== G E N E R A L   P U R P O S E   R O U T I N E S  ==========  */
/*  ==========						                            ==========  */
/*  ======================================================================  */

void Make_VEPM(PolyPointList *_P, VertexNumList *_V, EqList *_E, PairMat PM){
  int i, j;
  for (i=0;i<_E->ne;i++) for (j=0;j<_V->nv;j++) 
    PM[i][j]=Eval_Eq_on_V(&_E->e[i],_P->x[_V->v[j]],_P->n);
}

#define	 LLong_EEV		(1)    /* 1 @ [4662 4 20 333 422 1554 2329] */
#define  TEST_EEV	      	(0)	       /* compare Long to LLong EEV */

Equation EEV_To_Equation(Equation *_E1, Equation *_E2, Long *_V, int n){
  /* Calculate the equation spanned by _V and the intersection of _E1, _E2  */
  int i; Long l, m, g; Equation Eq;
  l=Eval_Eq_on_V(_E2,_V,n);
  m=Eval_Eq_on_V(_E1,_V,n);
  g=NNgcd(l,m); assert(g); l/=g; m/=g;
#if ((!(LLong_EEV))||(TEST_EEV))			    /* Long version */
  for(i=0;i<n;i++) Eq.a[i]=l*_E1->a[i]-m*_E2->a[i];
  { int gcd=Eq.c=l*_E1->c-m*_E2->c;
    for(i=0;i<n;i++) gcd=NNgcd(gcd,Eq.a[i]); assert(gcd);
    if (gcd!=1) { for(i=0;i<n;i++) Eq.a[i]/=gcd; Eq.c/=gcd;}}
#endif
#if ((LLong_EEV)||(TEST_EEV))				   /* LLong version */
  { LLong A[POLY_Dmax], C, G; for(i=0;i<n;i++) 
	A[i]=((LLong) l)*((LLong)_E1->a[i])-((LLong)m)*((LLong)_E2->a[i]);
    G=C=((LLong) l)*((LLong)_E1->c)-((LLong)m)*((LLong)_E2->c);
    for(i=0;i<n;i++) G=LNNgcd(G,A[i]); assert(G);
    if(G!=1) {C/=G; for(i=0;i<n;i++) A[i]/=G;}
#if	(TEST_EEV)						 /* Compare */
    {	int e=(Eq.c!=C); for(i=0;i<n;i++) if(Eq.a[i]!=A[i]) e=1;     
	if(e) { printf("Error in EEV: l=%d m=%d g=%d\n",l,m,g);
	for(i=0;i<n;i++)printf("%d ",_E1->a[i]);printf("  %d = E1\n",_E1->c);
	for(i=0;i<n;i++)printf("%d ",_E2->a[i]);printf("  %d = E2\n",_E2->c);
	for(i=0;i<n;i++)printf("%d ",Eq.a[i]);printf("  %d = Eq\n",Eq.c);
	for(i=0;i<n;i++)printf("%d ",A[i]);printf("  %d = LL_Eq\n",C);
	exit(0); }
    }
#else
    Eq.c=C; for(i=0;i<n;i++) Eq.a[i]=A[i];
#endif
  }
#endif
  return Eq;
}

Long Eval_Eq_on_V(Equation *E, Long *V, int i){    
  Long p=E->c; while(i--) p+=V[i]*E->a[i];
  return p;
}

int Vec_Greater_Than(Long *X, Long *Y, int i){	    /* return 1 iff `X > Y' */
  while(i--) {if(X[i]>Y[i]) return 1; if(X[i]<Y[i]) return 0;} 
  return 0;
}

int Vec_Equal(Long *X, Long *Y, int i){	    /* return 1 iff `X == Y' */
  while(i--) if(X[i]!=Y[i]) return 0;
  return 1;
}

int  IsGoodCEq(Equation *_E, PolyPointList *_P, VertexNumList *_V){
  int i=_V->nv; 
  Long s;
  while(!(s=Eval_Eq_on_V(_E, _P->x[_V->v[--i]], _P->n))); 
  if(s < 0) { int j=_P->n; while(j--) _E->a[j]=-_E->a[j]; _E->c=-_E->c; }
  while(i) if(Eval_Eq_on_V(_E, _P->x[_V->v[--i]], _P->n) < 0) return 0;
  return 1;
}

int  Search_New_Vertex(Equation *_E, PolyPointList *_P){
  int i, v=0; 
  Long *X=_P->x[0], x=Eval_Eq_on_V(_E,X,(_P->n));
  for(i=1;i<_P->np;i++)      {  
    Long *Y=_P->x[i], y=Eval_Eq_on_V(_E,Y,(_P->n));
    if(y>x) continue;
    if(y==x) if(Vec_Greater_Than(X,Y,_P->n)) continue;
    v=i; X=Y; x=y;     }
  return v;
}

int EL_to_PPL(EqList *_E, PolyPointList *_P, int *n){
  int i,j;
  for (i=0;i<_E->ne;i++){
    if (_E->e[i].c!=1) {_P->np=0; return 0;}
    for (j=0;j<*n;j++) _P->x[i][j]=_E->e[i].a[j];}
  _P->np=_E->ne;
  _P->n=*n;
  return 1;
}

/*  ======================================================================  */
/*  ==========		     			  	                     	==========  */
/*  ==========          S T A R T -- S I M P L E X              ==========  */
/*  ==========		     			  	                    	==========  */
/*  ======================================================================  */

/*	return 0 <=> max.dim., E.ne==P.n+1, made Simplex of Vertices of P;  *
 *      return (P.n-E.ne) == codim. > 0  <=> E.ne defining equations on E;  */

#define  VERT_WITH_MAX_DISTANCE (0)    /* 0 @ [1845 2 15 97 247 610 874]    */
#define	 LONG_EQ_FIRST		(0)    /* 0 @ [3425 2 7 137 429 1141 1709]  */
#define	 TEST_GLZ_EQ		(0)		 /* trace StartSimplex EQs  */

Long VZ_to_Base(Long *V,int *d,Long M[POLY_Dmax][POLY_Dmax])  /* 0 iff V=0 */
{    int p[POLY_Dmax], i, j, J=0; Long g=0, W[POLY_Dmax], *G[POLY_Dmax]; 
     for(i=0;i<*d;i++) 	if(V[i]) {W[J]=V[i]; G[J]=M[i]; p[J++]=i;}
			else for(j=0;j<*d;j++) M[i][j]=(i==j);
     if(J) if(p[0]) { G[0]=M[0]; for(j=0;j<*d;j++) M[p[0]][j]=(j==0);}
     if(J>1) g=W_to_GLZ(W,&J,G); else if(J){g=*W; M[0][0]=0; M[0][p[0]]=1;}
     if(J>1)
     {  for(i=0;i<J;i++) { int I=J; 
	for(j=*d-1;j>=0;j--) G[i][j] = (V[j]) ? G[i][--I] : 0; assert(I==0);}
     }	return g;
}

int  OrthBase_red_by_V(Long *V, int *d, Long A[][POLY_Dmax], int *r,
	Long B[][POLY_Dmax])
{    int i, j, k; Long W[POLY_Dmax], G[POLY_Dmax][POLY_Dmax];
     for(i=0;i<*r;i++) {int j; W[i]=0; for(j=0;j<*d;j++) W[i]+=A[i][j]*V[j];}
     assert( VZ_to_Base(W,r,G) );
     for(i=0;i<*r-1;i++) for(k=0;k<*d;k++)
     {	B[i][k]=0; for(j=0;j<*r;j++) B[i][k]+=G[i+1][j]*A[j][k];
     }
#if	(TEST_GLZ_EQ)
	printf("A -> B ... V = "); for(k=0;k<*d;k++) printf(" %5d",V[k]);
	printf("  W=");for(k=0;k<*r;k++)printf(" %5d",W[k]);puts(""); {int 
	a,b; for(a=0;a<*r-1;a++){for(b=0;b<*d;b++)printf(" %5d",A[a][b]);
	printf("  =A  B=  ");for(b=0;b<*d;b++)printf(" %5d",B[a][b]);
	puts("");}for(b=0;b<*d;b++)printf(" %5d",A[a][b]);printf("  =A\n");}
#endif
     	return (*r)--;
}

int  New_Start_Vertex(Long *V0,Long *Ea, PolyPointList *P,int *v) /* P.x[v] */
{    Equation E; 
     int i, n=0, p=0; 
     Long d, dn=0, dp=0, *Xn=P->x[0], *Xp=Xn; 
     
     for(i=0;i<P->n;i++) E.a[i]=Ea[i];
     
     E.c=0; 
     E.c=-Eval_Eq_on_V(&E,V0,P->n);
     
     d=Eval_Eq_on_V(&E,P->x[0],P->n); 
     
     if(d>0) dp=d; 
     if(d<0) dn=d;
     
     for(i=1;i<P->np;i++)
     {	d=Eval_Eq_on_V(&E,P->x[i],P->n); 
        if(d==0) continue; 
	    if(d==dp) if(Vec_Greater_Than(P->x[i],Xp,P->n)) Xp=P->x[p=i];
        if(d>dp)  {dp=d; Xp=P->x[p=i];}
	    if(d==dn) if(Vec_Greater_Than(P->x[i],Xn,P->n)) Xn=P->x[n=i];
        if(d<dn)  {dn=d; Xn=P->x[n=i];}
     }
     
     if(dp) 
     	if(dn) 				 /* points on both sides */
			#if	(VERT_WITH_MAX_DISTANCE)
				{if(dp+dn>0) *v=p; else *v=n;}
			#else
				{if(dp+dn>0) *v=n; else *v=p;}
			#endif
	 	else   *v=p;				/* d >=0 */
     else if(dn) *v=n;				/* d <=0 */
          else return 0;

     return 1;
}
	
int  GLZ_Start_Simplex(PolyPointList *_P, VertexNumList *_V, CEqList *_C)
{    
	/* returns codimension */
    int i, x=0, y=0, *VN=_V->v, *d=&_P->n, r=*d, b[POLY_Dmax]; 
    
    Long *X=_P->x[x], *Y=_P->x[y], XX=0, YY=0, B[(POLY_Dmax*(POLY_Dmax+1))/2][POLY_Dmax], W[POLY_Dmax]; 
	
	if(_P->np<2) 
     {	for(x=0;x<_P->n;x++) for(y=0;y<_P->n;y++) _C->e[x].a[y]=(x==y);
	assert(_P->np>0); for(x=0;x<_P->n;x++) _C->e[x].c=-_P->x[0][x];
	return _C->ne=_P->n;
     }
     
    for(i=1; i<_P->np; i++) 
    {	Long *Z=_P->x[i]; 	
		if(Vec_Greater_Than(X,Z,_P->n)) X=_P->x[x=i];	/* (x_n)-max: VN[0] */
		if(Vec_Greater_Than(Z,Y,_P->n)) Y=_P->x[y=i];	/* (x_n)-min: VN[1] */
    }	assert(x!=y);	     /* at this point I need two different vertices */\
     
    for(i=0;i<*d;i++) 
    {	Long Xi=(X[i]>0) ? X[i]: -X[i],
		Yi=(Y[i]>0) ? Y[i]: -Y[i]; 
		if(Xi>XX)XX=Xi; 
		if(Yi>YY)YY=Yi;
	}
	
    if(YY<XX) {VN[0]=y;VN[1]=x;} else {VN[0]=x;VN[1]=y;} 
    
    _V->nv=2; 
    y=VN[1];
    X=_P->x[VN[0]];
    Y=_P->x[VN[1]]; 
    
    for(i=0;i<*d;i++) b[i]=(i*(2*(*d)-i+1))/2;
    
    for(x=0;x<*d;x++) for(i=0;i<*d;i++) B[x][i]=(x==i); /* b[i+1]-b[i]=d-i */
    
    for(x=1;x<*d;x++)
    {	for(i=0;i<*d;i++) W[i]=_P->x[y][i]-X[i];
		OrthBase_red_by_V(W,d,&B[b[x-1]],&r,&B[b[x]]); 
		for(i=0;i<r;i++) 
			#if	(LONG_EQ_FIRST)
				if(New_Start_Vertex(X,B[b[x]+r-i-1],_P,&y)) break;
			#else
				if(New_Start_Vertex(X,B[b[x]+i],_P,&y)) break;
			#endif
			if(i==r) break;
		_V->v[_V->nv++]=y;	       /* x = dim(span) < d */
    }
     
    if(x<*d)
    {	for(y=0;y<r;y++)
		{	Equation *E=&_C->e[y]; 
			Long *Z=B[b[x]+y]; 
	    	E->c=0; 
	    	for(i=0;i<*d;i++) E->a[i]=Z[i]; 
	    	E->c=-Eval_Eq_on_V(E,X,_P->n); 
		}   
		return _C->ne=r;
    }
    else
    {	Equation *E=_C->e; 
     	Long *Z=B[b[*d-1]]; 
     	E->c=0; 
     	_C->ne=2;
		for(i=0;i<*d;i++) E->a[i]=Z[i];
		E->c=-Eval_Eq_on_V(E,X,_P->n); 
		if(Eval_Eq_on_V(E,_P->x[_V->v[*d]],_P->n)<0) {
			for(i=0;i<*d;i++)E->a[i]=-Z[i];
	  		E->c*=-1;}
		X=_P->x[_V->v[r=*d]]; 
        for(x=1;x<*d;x++)			/* now the 2nd equation */
		{	Y=_P->x[_V->v[x-1]]; 
			for(i=0;i<*d;i++) W[i]=X[i]-Y[i];
	    	OrthBase_red_by_V(W,d,&B[b[x-1]],&r,&B[b[x]]);
		}
		E=&_C->e[1]; E->c=0;
		for(i=0;i<*d;i++) E->a[i]=Z[i];
		E->c=-Eval_Eq_on_V(E,X,_P->n); 
		assert(XX=Eval_Eq_on_V(E,_P->x[_V->v[*d-1]],_P->n)); 
		if(XX<0) {for(i=0;i<*d;i++) E->a[i]=-Z[i]; E->c*=-1;}
        for(x=*d-2;x>=0;x--)			/* omit vertex #x */
		{	r=*d-x; 
			for(y=x+1;y<*d;y++)
	    	{	Y=_P->x[_V->v[y]]; 
	    		for(i=0;i<*d;i++) W[i]=X[i]-Y[i];
	    		OrthBase_red_by_V(W,d,&B[b[y-1]],&r,&B[b[y]]);
	    	}
			E=&_C->e[(_C->ne)++]; E->c=0;
			for(i=0;i<*d;i++) E->a[i]=Z[i];
			E->c=-Eval_Eq_on_V(E,X,_P->n); 
			assert(XX=Eval_Eq_on_V(E,_P->x[_V->v[x]],_P->n)); 
			if(XX<0) {for(i=0;i<*d;i++)E->a[i]=-Z[i]; E->c*=-1;}
		}
    } 
     
    assert(*d+1==_C->ne); 
     
    for(x=0;x<_C->ne;x++) 
     	for(i=0;i<=*d;i++)
     		assert((x==i)==(0!=Eval_Eq_on_V(&_C->e[x],_P->x[_V->v[*d-i]],_P->n)));
     
    return 0;
}


/*  ======================================================================  */
/*  ==========		     			  		                    ==========  */
/*  ==========	     P O L Y H E D R O N   A N A L Y S I S  	==========  */
/*  ==========					                        		==========  */ 
/*  ======================================================================  */

void Make_New_CEqs(PolyPointList *_P, VertexNumList *_V, CEqList *_C, 
			EqList *_F, INCI *CEq_I, INCI *F_I){
  int i,j, Old_C_ne=_C->ne;
  static CEqList Bad_C;
  static INCI Bad_C_I[CEQ_Nmax];

  Bad_C.ne=_C->ne=0;
  for (i=0;i<Old_C_ne;i++){
    Long dist = Eval_Eq_on_V(&_C->e[i],_P->x[_V->v[_V->nv-1]],_P->n);
    CEq_I[i]=INCI_PN(CEq_I[i],dist);
    if (dist<0) {Bad_C.e[Bad_C.ne]=_C->e[i]; Bad_C_I[Bad_C.ne++]=CEq_I[i];}
    else {_C->e[_C->ne]=_C->e[i]; CEq_I[_C->ne++]=CEq_I[i];}}

  Old_C_ne=_C->ne;
  for (i=0;i<_F->ne;i++) F_I[i]=
	INCI_PN(F_I[i],Eval_Eq_on_V(&_F->e[i],_P->x[_V->v[_V->nv-1]],_P->n));
  for (j=0;j<_F->ne;j++) if (!INCI_M2(F_I[j]))
    for (i=0;i<Bad_C.ne;i++){
      INCI New_Face=INCI_AND(Bad_C_I[i],F_I[j]);
      int k;
      if (INCI_abs(New_Face)<_P->n-1) continue;
      for (k=0;k<Bad_C.ne;k++) if (INCI_LE(New_Face,Bad_C_I[k])) 
	if (k!=i) break;
      if (k!=Bad_C.ne) continue;
      for (k=0;k<Old_C_ne;k++) if (INCI_LE(New_Face,CEq_I[k])) break;
      if (k!=Old_C_ne) continue;
      for (k=0;k<_F->ne;k++) if (INCI_LE(New_Face,F_I[k])) if (k!=j) break;
      if (k!=_F->ne) continue; 
      assert(_C->ne<CEQ_Nmax);
      CEq_I[_C->ne]=INCI_PN(INCI_D2(New_Face),0);
      _C->e[_C->ne]=EEV_To_Equation(&(Bad_C.e[i]),&(_F->e[j]),
				    _P->x[_V->v[_V->nv-1]],_P->n);
      assert(IsGoodCEq(&(_C->e[_C->ne++]),_P,_V));}
  for (j=0;j<Old_C_ne;j++) if (!INCI_M2(CEq_I[j])) 
    for (i=Bad_C.ne-1;i>=0;i--){
      INCI New_Face=INCI_AND(Bad_C_I[i],CEq_I[j]);
      int k;
      if (INCI_abs(New_Face)<_P->n-1) continue;
      for (k=0;k<Bad_C.ne;k++) if (INCI_LE(New_Face,Bad_C_I[k])) 
	if (k!=i) break;
      if (k!=Bad_C.ne) continue;
      for (k=0;k<Old_C_ne;k++) if (INCI_LE(New_Face,CEq_I[k]))
	if (k!=j) break;
      if (k!=Old_C_ne) continue;
      for (k=0;k<_F->ne;k++) if (INCI_LE(New_Face,F_I[k])) break;
      if (k!=_F->ne) continue;
      assert(_C->ne<CEQ_Nmax);
      CEq_I[_C->ne]=INCI_PN(INCI_D2(New_Face),0);
      _C->e[_C->ne]=EEV_To_Equation(&(Bad_C.e[i]),&(_C->e[j]),
				    _P->x[_V->v[_V->nv-1]],_P->n);
      assert(IsGoodCEq(&(_C->e[_C->ne++]),_P,_V));}
}

int  FE_Search_Bad_Eq(CEqList *_C, EqList *_F, INCI *CEq_I, INCI *F_I,
	       PolyPointList *_P, int *_IP){   /* return 0 :: no bad eq. */
  while(_C->ne--)  {	
    int j, M=_C->ne; 	/* INCI_LmR INCI_lex_GT */
    for(j=0;j<_C->ne;j++) if(INCI_lex_GT(&CEq_I[j],&CEq_I[M])) M=j;
    for(j=0;j<_P->np;j++)			
      if(Eval_Eq_on_V(&(_C->e[M]),_P->x[j],_P->n) < 0) {
	INCI AI=CEq_I[M]; Equation AE=_C->e[M]; 
	CEq_I[M]=CEq_I[_C->ne]; _C->e[M]=_C->e[_C->ne];
	CEq_I[_C->ne]=AI; _C->e[_C->ne]=AE; return ++_C->ne;}
    if(_C->e[M].c < 1) *_IP=0;
    assert(_F->ne<EQUA_Nmax);
    /* printf("#Feq=%d  #Ceq=%d\n",_F->ne,_C->ne); fflush(stdout); */
    _F->e[_F->ne]=_C->e[M]; F_I[_F->ne++]=CEq_I[M];
    if(M<_C->ne) {_C->e[M]=_C->e[_C->ne]; CEq_I[M]=CEq_I[_C->ne];}
    }
  return 0;
}

int  Finish_Find_Equations(PolyPointList *_P, VertexNumList *_V, 
		     EqList *_F, CEqList *_CEq, INCI *F_I, INCI *CEq_I){
  int IP=1;
  while(0<=_CEq->ne) if (FE_Search_Bad_Eq(_CEq,_F,CEq_I,F_I,_P,&IP)){
    assert(_V->nv<VERT_Nmax);
    _V->v[_V->nv++]=Search_New_Vertex(&(_CEq->e[_CEq->ne-1]),_P);
    Make_New_CEqs(_P,_V,_CEq,_F,CEq_I,F_I); }
  return IP;					   
}

int  Find_Equations(PolyPointList *_P, VertexNumList *_V, EqList *_F){
  /* returns IP and finds vertices and equations for _P */
  int i; 
  
  CEqList *CEq = (CEqList *) malloc(sizeof(CEqList)); 
  INCI *CEq_I = (INCI *) malloc(sizeof(INCI)*CEQ_Nmax);
  INCI *F_I = (INCI *) malloc(sizeof(INCI)*EQUA_Nmax); 
  
  CEq->ne=0;
  
  if (GLZ_Start_Simplex(_P, _V, CEq)) {
    _F->ne=CEq->ne; 
    for(i=0;i<_F->ne;i++) _F->e[i]=CEq->e[i]; 
    free(CEq); free(CEq_I); free(F_I); 
    return 0;}
    
  _F->ne=0;
  
  for (i=0;i<CEq->ne;i++) 
    if(INCI_abs(CEq_I[i]=Eq_To_INCI(&(CEq->e[i]),_P,_V))<_P->n)
      {exit(0);}
      
  i=Finish_Find_Equations(_P, _V, _F, CEq, F_I, CEq_I);
  
  free(CEq); free(CEq_I); free(F_I);
  
  return i;
}


/*  ======================================================================  */
/*  ==========		     			  	                     	==========  */
/*  ==========	  D U A L   P O L Y   &   C O M P L E T I O N 	==========  */
/*  ==========						                        	==========  */ 
/*  ======================================================================  */

void add_for_completion(Long *yDen, Long Den, EqList *_E, PolyPointList *_CP, int *old_np){
  int i,n=_CP->n;
  Long yold[POLY_Dmax];

  if(Den>1) for(i=0;i<n;i++) {
    if(yDen[i]%Den) return;
    yold[i]=yDen[i]/Den;}
  else for(i=0;i<n;i++) yold[i]=yDen[i];
  for (i=0;i<_E->ne;i++) if (Eval_Eq_on_V(&(_E->e[i]), yold, n) < 0) return;
  for (i=0;i<*old_np;i++) if (Vec_Equal(_CP->x[i],yold,n)) return;
  assert(_CP->np<POINT_Nmax);
  for(i=0;i<n;i++) _CP->x[_CP->np][i]=yold[i];
  _CP->np++;
}

void Compute_InvMat(int n, EqList *_E, int OrdFac[VERT_Nmax],
		    int BasFac[POLY_Dmax], Long *Den,
		    Long InvMat[POLY_Dmax][POLY_Dmax]){
  /* Find first POLY_Dmax linearly independent facets + Inverse Matrix */

  LRat ind[POLY_Dmax][POLY_Dmax], x[POLY_Dmax], y[POLY_Dmax], f, 
    PInvMat[POLY_Dmax][POLY_Dmax];
  int i, j, k, l, rank=0, one[POLY_Dmax];
  
  for (i=0;i<n;i++) for (j=0;j<n;j++) PInvMat[i][j]=LrI(0);
  for (i=0;i<n;i++) PInvMat[i][i]=LrI(1);
  i=0;
  while (rank<n){
    for (j=0;j<n;j++) x[j]=LrI(_E->e[OrdFac[i]].a[j]);
    for (j=0;j<n;j++) y[j]=LrI(0);
    y[rank]=LrI(1);
    for (j=0;j<rank;j++) {
      f=x[one[j]];
      for (k=0;k<n;k++) {
        x[k]=LrD(x[k],LrP(f,ind[j][k])); 
        y[k]=LrD(y[k],LrP(f,PInvMat[j][k]));  } }
    one[rank]=-1;
    for (l=0;(l<n)&&(one[rank]==-1);l++) if (x[l].N) one[rank]=l;
    if(one[rank]>-1){
      for (k=0;k<n;k++) {
        ind[rank][k]=LrQ(x[k],x[one[rank]]);
        PInvMat[rank][k]=LrQ(y[k],x[one[rank]]); }
      for (j=0;j<rank;j++) {
        f=ind[j][one[rank]];
        for (k=0;k<n;k++)         {
          ind[j][k]=LrD(ind[j][k],LrP(ind[rank][k],f));   
          PInvMat[j][k]=LrD(PInvMat[j][k],LrP(PInvMat[rank][k],f));  }     }
      BasFac[rank]=OrdFac[i];
      rank++; }  
    i++; }
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    *Den=(*Den/LFgcd(*Den,PInvMat[i][j].D))*PInvMat[i][j].D;
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    InvMat[one[i]][j]=(*Den/PInvMat[i][j].D)*PInvMat[i][j].N;

  for (i=0;i<n;i++){
    for (j=0;j<n;j++) {
      long long s=0;
      for(k=0;k<n;k++) s+=((long long) (InvMat[k][i]))*
			 ((long long) (_E->e[BasFac[j]].a[k]));
      if (s!=*Den*(i==j)) {
	puts("something wrong in Make_Dual_Poly");
	exit(0);}}} 
}

void Complete_Poly(PairMat VPM, EqList *_E, int nv, PolyPointList *_CP){
  int i, j, k, InsPoint, n=_CP->n, old_np=_CP->np;
  Long MaxDist[EQUA_Nmax], InvMat[POLY_Dmax][POLY_Dmax], Den=1;
  Long yDen[POLY_Dmax];
  int OrdFac[VERT_Nmax], BasFac[POLY_Dmax], position[POLY_Dmax];

  /* Calculate maximal distances from facets of Delta^* (Vertices of Delta) */

  for (i=0;i<_E->ne;i++) {  
    MaxDist[i]=0;
    for (j=0;j<nv;j++) 
    if (MaxDist[i]<VPM[i][j]) MaxDist[i]=VPM[i][j];}

  /* Order facets of Delta^* (Vertices of Delta) w.r.t. MaxDist   */

  OrdFac[0]=0;
  for (i=1;i<_E->ne;i++){
    InsPoint=i; 
    while (InsPoint&&(MaxDist[i]<MaxDist[OrdFac[InsPoint-1]])) InsPoint--;
    for (j=i;j>InsPoint;j--) OrdFac[j]=OrdFac[j-1];
    OrdFac[InsPoint]=i; }

  Compute_InvMat(n, _E, OrdFac, BasFac, &Den, InvMat);
  
  /* Examine all integer points of parallelogram:                         */
  /* The basic structure of the algorithm is:
  for (k=0;k<n-1;k++) position[k]=-1;      / * sets k=n-1; important!      *
  position[n-1]=-2;  / * starting point just outside the parallelogram     *
  while(k>=0){
    position[k]++;
    DO AT position;
    for(k=n-1;((position[k]==MaxDist[BasFac[k]]-1)&&(k>=0));k--) 
       position[k]=-1;  }
         / * sets k to the highest value where pos.[k] wasn't the max value; 
            resets the following max values to min values                 */
  /* Quantities linear in position can be changed with every change of
     position (here: yDen)                                             */
  
  for(i=0;i<n;i++) yDen[i]=0;
  for (k=0;k<n-1;k++) {   /* sets k=n-1; important!   */
    position[k]=-_E->e[BasFac[k]].c;   
    for(i=0;i<n;i++) yDen[i]-=_E->e[BasFac[k]].c*InvMat[i][k]; }
  position[n-1]=-_E->e[BasFac[n-1]].c-1;
  for(i=0;i<n;i++) yDen[i]-=(_E->e[BasFac[k]].c+1)*InvMat[i][n-1];
  while(k>=0){
    position[k]++;
    for(i=0;i<n;i++) yDen[i]+=InvMat[i][k];
    add_for_completion(yDen, Den, _E, _CP, &old_np);
    for(k=n-1;(k>=0);k--){
      if (position[k]!=MaxDist[BasFac[k]]-_E->e[BasFac[k]].c) break;
      position[k]=-_E->e[BasFac[k]].c;
      for (i=0;i<n;i++) yDen[i]-=MaxDist[BasFac[k]]*InvMat[i][k]; }}
      
}
