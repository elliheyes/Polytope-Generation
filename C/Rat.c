#include "Global.h"
#include "Rat.h"

Rat  rI(Long a) 			       		      /*  a -> a/1  */
{    Rat c; c.N=a; c.D=1; return c; 
}
Rat  rR(Long a, Long b)					    /* (a,b) -> a/b */
{    Rat c; register Long g=Fgcd(a,b); 
     if( (c.D=(b/g)) < 0 ) {c.D=-c.D; c.N=-a/g;} else c.N=a/g; 
#ifdef	TEST
     if(c.D<=0) {Rpr(c); puts(" *** c.D<=0 in rR! ***"); exit(0);}
#endif
     return c;
} 
Rat  irP(Long i, Rat c)					   /* (int) * (Rat) */
{    register Long g=Fgcd(i,c.D); if(g < 0) g=-g; 
     c.N *= (i/g); c.D/=g; 
#ifdef	TEST
     if(c.D<=0) {Rpr(c); puts(" *** c.D<=0 in irP! ***"); exit(0);}
#endif
     return c;
} 
Rat  rS(Rat a, Rat b)                                            /*  a + b  */
{    Rat c; register Long g=Fgcd(a.D,b.D); 
     g = Fgcd(c.N=a.N*(b.D/g)+b.N*(a.D/g), c.D=a.D*(b.D/g));
     if(g < 0) g=-g;
     c.N/=g; c.D/=g; /* LONG in line above ... */
#ifdef	TEST
     if(c.D<=0) {Rpr(c); puts(" *** c.D<=0 in rS! ***"); exit(0);}
#endif
     return c;
}
Rat  rD(Rat a, Rat b)                                            /*  a - b  */
{    Rat c; register Long g=Fgcd(a.D,b.D); 
     g = Fgcd(c.N=a.N*(b.D/g)-b.N*(a.D/g), c.D=a.D*(b.D/g));
     if(g < 0) g=-g;
     c.N/=g; c.D/=g; /* LONG in line above ... */
#ifdef	TEST
     if(c.D<=0) {Rpr(c); puts(" *** c.D<=0 in rD! ***"); exit(0);}
#endif
     return c;
}
Rat  rP(Rat a, Rat b)                                            /*  a * b  */
{    register Long g=Fgcd(a.N,b.D); register Long h=Fgcd(b.N,a.D);
     Rat c; c.N=(a.N/g)*(b.N/h); 
     if((c.D=(a.D/h)*(b.D/g)) < 0) {c.D=-c.D; c.N=-c.N;} return c;
}
Rat  rQ(Rat a, Rat b)                                            /*  a / b  */
{    register Long g=NNgcd(a.N,b.N); register Long h=Fgcd(b.D,a.D);
     Rat c; c.N=(a.N/g)*(b.D/h);
     if((c.D=(a.D/h)*(b.N/g)) < 0) {c.N=-c.N; c.D=-c.D;} return c;
}
int  rC(Rat a, Rat b)          /* Compare = [1 / 0 / -1] iff a [gt/eq/lt] b */
{    register Long g=Fgcd(a.D,b.D);
     if((g=a.N*(b.D/g)-b.N*(a.D/g))) {if(g>0) return 1; else return -1;} 
     else return 0;
}
void Rpr(Rat c)					/* write "c.N/c.D" -> outFN */
{    printf("%d/%d",(int) c.N, (int) c.D);
} 
#ifdef	TEST
void TEST_Rat()
{    Rat a, b; int i, j, k, l;
     puts(" type 4 numbers: "); scanf("%d%d%d%d",&i,&j,&k,&l);
     Rpr(rS(rR(i,j),rR(k,l))); puts(" i/j + k/l");
     Rpr(rD(rR(i,j),rR(k,l))); puts(" i/j - k/l");
}
#endif
/*  ==========  	  END of Rational Operations		==========  */



/*  ==========		(Extended) Greatest Common Divisor	 =========  */
/*
 *	Fgcd()  computes the GCD for two positive integers (F -> fast)
 *	NNgcd() computes the non-negative GCD for two arbitrary integers
 *
 *	Egcd()  computes the EGCD for two integers
 *	REgcd() recursively computes the EGCD for a number of integers
 */
Long Fgcd(register Long a, register Long b)     /* Fast greatest common div */
{
     while( a %= b ) if ( !(b %= a) ) return a;
     return  b;
}
Long NNgcd(register Long a, register Long b)  /* NonNegative gcd handling 0 */
{
     a = (a<0) ? -a:a; b = (b<0) ? -b:b; if (!b) return a;
     while( a %= b ) if ( !(b %= a) ) return a;
     return  b;
}
/*  ==========	Egcd(a,b,*x,*y)=a *x+b *y;  REgcd(Vin, *dim, Vout) = Vin.Vout; 
 *   
 *   a_i+1=a_i-1 % a_i, a_I+1==0  =>  egcd(a_0,a_1)=a_I;    Vout={x_I,y_I} 
 *   a_i==x_ia_0+y_ia_1a_i+1 and a_i+1=a_i-1 -(a_i-1 /a_i)a_i => 
 *   x_i+1=x_i-1 -(a_i-1 /a_i)x_i with x_0=1; x_1=0;  y_I=(a_I-a_0 x_I)/a_1
 *									    */
Long Egcd(register Long A0, register Long A1, Long *Vout0, Long *Vout1)  
{    register Long V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
     while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
     *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); return A1;
}
/*   
 *   g = v0 x0 + ... + v(n-1) x(n-1),	G= A g+ B vn	=>
 *   G = v0 x0 + ... + vn xn		with xn=B, xi*=A, 
 *									    */
Long REgcd(Long *Vin, int *_n, Long *Vout)  /*  recursive Egcd(a_1,...,a_n)  */
{    Long Ain[2], Aout[2], gcd; int N=*_n-1, i;
     if (*_n==2) return Egcd(*Vin,(Vin[1]),Vout,&(Vout[1]));
     *Ain=REgcd(Vin,&N,Vout); Ain[1]=Vin[N]; 	gcd=
	Egcd(*Ain,(Ain[1]),Aout,&(Aout[1]));	Vout[N]=Aout[1]; 
     for(i=0; i<N; i++) Vout[i] *= (*Aout);
     return gcd;
}
/*  ==========	   END of (Extended) Greatest Common Divisor	 =========  */


/*  ==========	 GLZ-matrix=(Egcd, B_1, B_2, ...) E.W=g, B.W=0	 =========  */
/*
 *   assumes W[i]!=0, lines 1 thru d-1 are lower triangular and		    *
 *   their diagonal entries are positive if all W[i] are positive.	    *
 *   GLZ has to be an initialized list of pointers. 	 returns gcd(W_i)   *
 *									    */
/*  F <= N/D < F+1  */
Long FloorQ(Long N,Long D)				
{    
	Long F; 
	if(D<0) {D=-D; N=-N;} 
	F=N/D; 
	return (F*D>N) ? F-1 : F; 
}

Long CeilQ(Long N,Long D)
{ 
	return -FloorQ(-N,D); 
}

Long RoundQ(Long N,Long D)
{    
	Long F; 
	if(D<0) {D=-D; N=-N;} 
	F=N/D; 
	return F+(2*(N-F*D))/D; 
}

Long W_to_GLZ(Long *W, int *d, Long **GLZ)		
{    
	int i, j; 
	Long G, *E=*GLZ, *B=GLZ[1]; 
	for(i=0;i<*d;i++) assert(W[i]!=0);
    for(i=1;i<*d;i++)for(j=0;j<*d;j++)GLZ[i][j]=0;
    G=Egcd(W[0],W[1],&E[0],&E[1]); 
    B[0]=-W[1]/G; 
    B[1]=W[0]/G;
    for(i=2;i<*d;i++){  
    	Long a, b, g=Egcd(G,W[i],&a,&b); 
    	B=GLZ[i];
        B[i]= G/g; 
        G=W[i]/g; 
        for(j=0;j<i;j++) B[j]=-E[j]*G;  /* B=Base-Line */
        for(j=0;j<i;j++) E[j]*=a;
		E[j]=b;                     /* next Egcd */
        for(j=i-1;0<j;j--){   
        	int n; 
        	Long *Y=GLZ[j], rB=RoundQ(B[j],Y[j]), rE=RoundQ(E[j],Y[j]);
            for(n=0;n<=j;n++){ 
            	B[n] -= rB*Y[n]; 
            	E[n] -= rE*Y[n]; 
            } 
		}   
		G=g;
     } 
     return G;
}

/*   Map Permutations: Do "ArgFun" for all permutations pi of *d elements */
#ifdef	__cplusplus
#define ARG_FUN		void (*ArgFun)(int *d,int *pi,int *pinv,void *info)
#else
#define	ARG_FUN 	void ( ArgFun() ) 
#endif
void Map_Permut(int *d,int *pi,int *pinv,ARG_FUN,void *AuxPtr)
{    int i, j, n_rem_perm, n_perm=1, a, b, perm_j;
     for (i=1;i<=*d;i++) n_perm*=i;
     for (i=0;i<n_perm;i++)
     {	b=i; n_rem_perm=n_perm; for (j=0;j<*d;j++) pinv[j]=-1;
	for (j=*d;j>0;j--)
	{   n_rem_perm/=j; a=b/n_rem_perm; b=b%n_rem_perm; perm_j=*d;
	    while (a>=0) { perm_j--; if (pinv[perm_j]<0) a--;}
	    pi[j-1]=perm_j; pinv[perm_j]=j-1;
	}   (ArgFun)(d,pi,pinv,AuxPtr);
     }
}
/*  ==========	   END of GLZ-matrix=(Egcd, B_1, B_2, ...)	 =========  */




LRat  LrI(LLong a) 			       		      /*  a -> a/1  */
{    LRat c; c.N=a; c.D=1; return c; 
}
LRat  LrR(LLong a, LLong b)					    /* (a,b) -> a/b */
{    LRat c; register LLong g=LFgcd(a,b); 
     if( (c.D=(b/g)) < 0 ) {c.D=-c.D; c.N=-a/g;} else c.N=a/g; 
#ifdef	TEST
     if(c.D<=0) {LRpr(c); puts(" *** c.D<=0 in rR! ***"); exit(0);}
#endif
     return c;
} 
LRat  LirP(LLong i, LRat c)					   /* (int) * (LRat) */
{    register LLong g=LFgcd(i,c.D); if(g < 0) g=-g; 
     c.N *= (i/g); c.D/=g; 
#ifdef	TEST
     if(c.D<=0) {LRpr(c); puts(" *** c.D<=0 in irP! ***"); exit(0);}
#endif
     return c;
} 
LRat  LrS(LRat a, LRat b)                                            /*  a + b  */
{    LRat c; register LLong g=LFgcd(a.D,b.D); 
     g = LFgcd(c.N=a.N*(b.D/g)+b.N*(a.D/g), c.D=a.D*(b.D/g));
     if(g < 0) g=-g;
     c.N/=g; c.D/=g; /* LLong in line above ... */
#ifdef	TEST
     if(c.D<=0) {LRpr(c); puts(" *** c.D<=0 in rS! ***"); exit(0);}
#endif
     return c;
}
LRat  LrD(LRat a, LRat b)                                            /*  a - b  */
{    LRat c; register LLong g=LFgcd(a.D,b.D); 
     g = LFgcd(c.N=a.N*(b.D/g)-b.N*(a.D/g), c.D=a.D*(b.D/g));
     if(g < 0) g=-g;
     c.N/=g; c.D/=g; /* LLong in line above ... */
#ifdef	TEST
     if(c.D<=0) {LRpr(c); puts(" *** c.D<=0 in rD! ***"); exit(0);}
#endif
     return c;
}
LRat  LrP(LRat a, LRat b)                                            /*  a * b  */
{    register LLong g=LFgcd(a.N,b.D); register LLong h=LFgcd(b.N,a.D);
     LRat c; c.N=(a.N/g)*(b.N/h); 
     if((c.D=(a.D/h)*(b.D/g)) < 0) {c.D=-c.D; c.N=-c.N;} return c;
}
LRat  LrQ(LRat a, LRat b)                                            /*  a / b  */
{    register LLong g=LNNgcd(a.N,b.N); register LLong h=LFgcd(b.D,a.D);
     LRat c; c.N=(a.N/g)*(b.D/h);
     if((c.D=(a.D/h)*(b.N/g)) < 0) {c.N=-c.N; c.D=-c.D;} return c;
}
int  LrC(LRat a, LRat b)          /* Compare = [1 / 0 / -1] iff a [gt/eq/lt] b */
{    register LLong g=LFgcd(a.D,b.D);
     if((g=a.N*(b.D/g)-b.N*(a.D/g))) {if(g>0) return 1; else return -1;} 
     else return 0;
}
void LRpr(LRat c)					/* write "c.N/c.D" -> outFN */
{    printf("%lld/%lld",(LLong) c.N, (LLong) c.D);
} 
#ifdef	TEST
void TEST_LRat()
{    LRat a, b; int i, j, k, l;
     puts(" type 4 numbers: "); scanf("%d%d%d%d",&i,&j,&k,&l);
     LRpr(LrS(rR(i,j),LrR(k,l))); puts(" i/j + k/l");
     LRpr(LrD(LrR(i,j),LrR(k,l))); puts(" i/j - k/l");
}
#endif
/*  ==========  	  END of LRational OpeLRations		==========  */



/*  ==========		(Extended) Greatest Common Divisor	 =========  */
/*
 *	LFgcd()  computes the GCD for two positive integers (F -> fast)
 *	LNNgcd() computes the non-negative GCD for two arbitrary integers
 *
 *	LEgcd()  computes the LEgcd for two integers
 *	LREgcd() recursively computes the LEgcd for a number of integers
 */
LLong LFgcd(register LLong a, register LLong b)     /* Fast greatest common div */
{
     while( a %= b ) if ( !(b %= a) ) return a;
     return  b;
}
LLong LNNgcd(register LLong a, register LLong b)  /* NonNegative gcd handling 0 */
{
     a = (a<0) ? -a:a; b = (b<0) ? -b:b; if (!b) return a;
     while( a %= b ) if ( !(b %= a) ) return a;
     return  b;
}
/*  ==========	LEgcd(a,b,*x,*y)=a *x+b *y; LREgcd(Vin, *dim, Vout) = Vin.Vout;
 *   
 *   a_i+1=a_i-1 % a_i, a_I+1==0  =>  LEgcd(a_0,a_1)=a_I;    Vout={x_I,y_I} 
 *   a_i==x_ia_0+y_ia_1a_i+1 and a_i+1=a_i-1 -(a_i-1 /a_i)a_i => 
 *   x_i+1=x_i-1 -(a_i-1 /a_i)x_i with x_0=1; x_1=0;  y_I=(a_I-a_0 x_I)/a_1
 *									    */
LLong LEgcd(register LLong A0, register LLong A1, LLong *Vout0, LLong *Vout1)  
{    register LLong V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
     while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
     *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); return A1;
}
/*   
 *   g = v0 x0 + ... + v(n-1) x(n-1),	G= A g+ B vn	=>
 *   G = v0 x0 + ... + vn xn		with xn=B, xi*=A, 
 *									    */
LLong LREgcd(LLong *Vin, int *_n, LLong *Vout)  /*  recursive LEgcd(a_1,...,a_n)  */
{    LLong Ain[2], Aout[2], gcd; int N=*_n-1, i;
     if (*_n==2) return LEgcd(*Vin,(Vin[1]),Vout,&(Vout[1]));
     *Ain=LREgcd(Vin,&N,Vout); Ain[1]=Vin[N]; 	gcd=
	LEgcd(*Ain,(Ain[1]),Aout,&(Aout[1]));	Vout[N]=Aout[1]; 
     for(i=0; i<N; i++) Vout[i] *= (*Aout);
     return gcd;
}
/*  ==========	   END of (Extended) Greatest Common Divisor	 =========  */
