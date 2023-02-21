#include "Global.h"
#include "Rat.h"

#define SL_Long		LLong		   

#define NFX_Limit       903		/* 138b->255  153e->279  165c->327 */
#define X_Limit         9999		   /* 178c->375  218->399,462,483 */
#define VPM_Limit	9999

#define TEST_GLZ_VS_SL	(0)		/* exit on difference: GL vs. SL */

#undef  WARN_BIG_NS			/* [1152] ...  1152 = sym(24-cell)  */

typedef struct {int nv, nf, ns;}  				vNF;

typedef struct {int C[VERT_Nmax], L[VERT_Nmax], s;}             PERM;

/*      =============================================================       */

int  GLZ_Make_Trian_NF(Long X[][VERT_Nmax], int *n, int *nv,
		       GL_Long G[POLYDIM][POLYDIM]);   
		       
int  SL2Z_Make_Poly_NF(Long X[][VERT_Nmax], int *n, int *nv,
		       SL_Long S[POLYDIM][POLYDIM]);   

int  Init_rVM_VPM(PolyPointList *P, VertexNumList *_V,EqList *_F,/* in */
	    	int *d,int *v,int *f, Long VM[POLYDIM][VERT_Nmax], /* out */
	    	Long VPM[VERT_Nmax][VERT_Nmax]);	/* return reflexive */
	    	
void Make_VPM_NF(int *v, int *f, Long VPM[VERT_Nmax][VERT_Nmax],      /* in */
		PERM *CL,int *ns,Long VPM_NF[VERT_Nmax][VERT_Nmax]); /* out */
		
void Aux_pNF_from_vNF(PERM *CL,int *ns,int *v,int *d,
		Long VM[POLYDIM][VERT_Nmax],			      /* in */
		Long pNF[POLYDIM][VERT_Nmax],int *t);		     /* out */
		
void New_pNF_Order(int *v,int *f,PERM *CL,int *ns,Long VPM_NF[][VERT_Nmax]);
	    
	
/*      =============================================================       */

	
void swap(int* i,int* j) {register int k; k=*i; *i=*j; *j=k;}

int  Aux_Make_Poly_NF(Long X[][VERT_Nmax], int *n, int *nv)
{    
	GL_Long G[POLYDIM][POLYDIM];
	#if	(TEST_GLZ_VS_SL)
    	int i,j,x;
     	SL_Long S[POLYDIM][POLYDIM];
     	Long XS[POLYDIM][VERT_Nmax];
     	for(i=0;i<*n;i++)for(j=0;j<*nv;j++)XS[i][j]=X[i][j];
     	x=GLZ_Make_Trian_NF(X,n,nv,G); 
     	SL2Z_Make_Poly_NF(XS,n,nv,S);
     	for(i=0;i<*n;i++)for(j=0;j<*n;j++)  assert( S[i][j]==G[i][j]);}
     	for(i=0;i<*n;i++)for(j=0;j<*nv;j++) assert(XS[i][j]==X[i][j]);
		return x;
	#else
		return GLZ_Make_Trian_NF(X,n,nv,G);
	#endif
}


GL_Long GL_Egcd(GL_Long A0, GL_Long A1, GL_Long *Vout0, GL_Long *Vout1)  
{    
	GL_Long V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
    while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
    *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); return A1;
}


GL_Long GL_RoundQ(GL_Long N,GL_Long D)
{    
	GL_Long F; if(D<0) {D=-D; N=-N;} F=N/D; return F+(2*(N-F*D))/D; 
}


GL_Long GL_W_to_GLZ(GL_Long *W, int d, GL_Long **GLZ)		
{    
	int i, j; GL_Long G, *E=*GLZ, *B=GLZ[1];for(i=0;i<d;i++)assert(W[i]!=0);
    for(i=1;i<d;i++)for(j=0;j<d;j++)GLZ[i][j]=0;
    G=GL_Egcd(W[0],W[1],&E[0],&E[1]); B[0]=-W[1]/G; B[1]=W[0]/G;
    for(i=2;i<d;i++){  
    	GL_Long a, b, g=GL_Egcd(G,W[i],&a,&b); 
    	B=GLZ[i];
        B[i]= G/g; 
        G=W[i]/g; 
        for(j=0;j<i;j++) B[j]=-E[j]*G;  /* B=Base-Line */
        for(j=0;j<i;j++) E[j]*=a;
		E[j]=b;                     /* next Egcd */
        for(j=i-1;0<j;j--){   
        	GL_Long *Y=GLZ[j],rB=GL_RoundQ(B[j],Y[j]),rE=GL_RoundQ(E[j],Y[j]);
            int n; 
            for(n=0;n<=j;n++) { B[n] -= rB*Y[n]; E[n] -= rE*Y[n]; } 
		}   
		G=g;
    } 
    return G;
}

int  GLZ_Make_Trian_NF(Long X[][VERT_Nmax], int *n, int *nv, GL_Long G[POLYDIM][POLYDIM])
{    
	int i, j, C=-1, L; 
	GL_Long g, W[POLYDIM], NF[POLYDIM][VERT_Nmax], 
	*_G[POLYDIM], NG[POLYDIM][POLYDIM];
	
	for(i=0;i<*n;i++)_G[i]=NG[i];
    
    for(i=0;i<*n;i++)for(j=0;j<*n;j++)G[i][j]=(i==j);		  /* init G */
    
    for(i=0;i<*n;i++)for(j=0;j<*nv;j++)NF[i][j]=0; 
    
    for(L=0;L<*n;L++){  
    	int N=0, p[POLYDIM];
    	/* find column C with N>0 non-zero entries */
        while(N==0){   
        	++C; 
        	for(i=0;i<*n;i++){	
        		for(j=0;j<*n;j++) NF[i][C]+=G[i][j]*X[j][C];
	   		}
	    	for(i=L;i<*n;i++) if(NF[i][C]) {W[N]=NF[i][C]; p[N++]=i;}
		}
        assert(N); 
        if(N==1) {g=W[0]; _G[0][0]=1;} 
        else g=GL_W_to_GLZ(W,N,_G);
		if(g<0) { g *= -1; for(i=0;i<N;i++)_G[0][i] *= -1; }
		NF[L][C]=g; 
		for(i=L+1;i<*n;i++) NF[i][C]=0;
        for(i=0;i<*n;i++){   
        	GL_Long Cp[POLYDIM]; 
        	for(j=0;j<N;j++) Cp[j]=G[p[j]][i];
	    	for(j=0;j<N;j++){	
	    		int k; 
	    		G[p[j]][i]=0; 
				for(k=0;k<N;k++) G[p[j]][i] += _G[j][k]*Cp[k];
	    	}
		}
        if(L!=p[0])
        	 /* swap lines G[L] <-> G[p[0]] */
        	for(i=0;i<*n;i++){
        		GL_Long A=G[L][i]; G[L][i]=G[p[0]][i]; G[p[0]][i]=A; 
			}
			/* make upper diag minimal nonneg. */
			for(i=0;i<L;i++){   
				GL_Long R=NF[i][C]/NF[L][C];
	    		if((NF[i][C]-R*NF[L][C])<0) R-=1;
	   			NF[i][C]-=R*NF[L][C]; 
	   			for(j=0;j<*n;j++)G[i][j]-=R*G[L][j];
			}
     }
     
     while(++C<*nv)for(i=0;i<*n;i++)for(j=0;j<*n;j++)NF[i][C]+=G[i][j]*X[j][C];
     
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++){ 
	 	#ifdef	SHOW_NFX_LIMIT
			g=NF[i][j]; 
			if(g<0) g=-g; 
			if(g>NFX_Limit) 
				{ fprintf(stderr,"NFX_Limit in GL -> %lld !!\n",(long long) g); return 0; } 
			else 
		#endif
		X[i][j]=NF[i][j]; 
	}
	
    return 1;
}


void TEST_rVM_VPM(int *d,int *v,int *f, Long X[POLYDIM][VERT_Nmax],
	    	Long x[VERT_Nmax][VERT_Nmax])
{    
	int i,j,err=0; 

	for(i=0;i<*v;i++){	
		for(j=0;j<*d;j++) if(abs(X[j][i])>X_Limit) err=X[j][i];
		for(j=0;j<*f;j++) if(abs(x[j][i])>VPM_Limit) err=x[j][i]; 
    }	
     
    if(err){	
    	printf("TEST_VM_VPM: limits exceeded %d\n",err);
		printf("%d %d VM[%d][%d]:\n",*v,*d,*d,*v);
		for(j=0;j<*d;j++){   
			for(i=0;i<*v;i++)printf("%3d ",(int) X[j][i]);puts("");
		}   
		puts("");
		printf("VPM[%d][%d]:\n",*f,*v);
		for(j=0;j<*f;j++){   
			for(i=0;i<*v;i++)printf("%3d ",(int) x[j][i]);puts("");
		}   
		puts("");
		exit(0);
     }
}


int  Init_rVM_VPM(PolyPointList *_P,VertexNumList *_V,EqList *_F,/* in */
	    	int *d,int *v,int *f, Long X[POLYDIM][VERT_Nmax],  /* out */
	    	Long x[VERT_Nmax][VERT_Nmax])		/* return reflexive */
{    
	int i,j, ref=1; 
	
	/* define the number of vertices, equations and points */
    *v=_V->nv; 
    *f=_F->ne; 
    *d=_P->n;
     
    /* compute VPM */
    for(j=0;j<_F->ne;j++) {
		if(_F->e[j].c!=1) ref=0; 
		for(i=0;i<_V->nv;i++)
        	x[j][i]=Eval_Eq_on_V(&_F->e[j],_P->x[_V->v[i]],_P->n);
    }
    
    for(i=0;i<_V->nv;i++){  
    	Long *pv=_P->x[_V->v[i]];
		for(j=0;j<_P->n;j++) X[j][i]=pv[j];
    }
    
    TEST_rVM_VPM(d,v,f,X,x);
    
    return ref;
}


void New_pNF_Order(int *v,int *f,PERM *CL,int *ns,Long VPM_NF[][VERT_Nmax])
{    
	int i, j;
	int pi[VERT_Nmax], c[VERT_Nmax]; 
    Long maxP[VERT_Nmax], sumP[VERT_Nmax];
    
    for(i=0;i<*v;i++){	
    	pi[i]=i; 
    	maxP[i]=sumP[i]=0; 
    	for(j=0;j<*f;j++){   
    		sumP[i]+=VPM_NF[j][i];
	    	if(VPM_NF[j][i]>maxP[i]) maxP[i]=VPM_NF[j][i];
		}	
    }
    
    for(i=0;i<*v-1;i++){	
    	int n=i; 
    	for(j=i+1;j<*v;j++){   
    		if(maxP[j]<maxP[n]) n=j; 
    		else if(maxP[j]==maxP[n]) if(sumP[j]<sumP[n]) n=j;
		}
        if(n!=i){  
        	Long aP=maxP[i]; 
        	int a=pi[i]; 
	    	maxP[i]=maxP[n]; 
	    	maxP[n]=aP; 
	    	pi[i]=pi[n]; 
	    	pi[n]=a;
	    	aP=sumP[i]; 
	    	sumP[i]=sumP[n]; 
	    	sumP[n]=aP; 
		}
    }
    
    for(i=0;i<*ns;i++){	
    	int *C=CL[i].C; 
    	for(j=0;j<*v;j++) c[j]=C[pi[j]]; 
		for(j=0;j<*v;j++) C[j]=c[j];
    }
}


void Aux_vNF_Line(int l,vNF *_X,Long x[][VERT_Nmax], PERM *CL,int *S,int *_ns)
{    
	int n=(*_ns), cf=0;	/*  cf=CompareFlag (ref. exists & o.k.      */
    Long *y, r[VERT_Nmax];	/*  r=ReferenceLine; y->X[line]		    */
    /*  go over CL (n from  *_ns-1  to  0), ns_* changes!  */
    while(n--){	
    	PERM nP[VERT_Nmax];	        
		int c=0, L=l-1, np=0, *C, ccf=cf;	/*  ccf=column compare flag */
		*nP=CL[n];
		/*  init nP (from 1st col.) */
		while(++L<_X->nf){   
			int j=0; C=nP[np].C; 
			y=x[nP[np].L[L]];
	    	while(++j<=*S) if(y[C[c]]<y[C[j]]) swap(&C[c],&C[j]);
	    	if(ccf){   
	    		Long d=y[*C]-*r;
	    		/* BAD line */
				if(d<0) ;			
				/* BETTER line */		
				else if(d){   
					*r=y[*C]; 
					cf=0; 
					*nP=nP[np]; 
					nP[np=1]=CL[n]; 
					*_ns=n+1;
		    		swap(&(nP->L[l]),&(nP->L[L]));
				}
				 /* EQUAL line */
				else{   
					swap(&(nP[np].L[l]),&(nP[np].L[L])); 
					nP[++np]=CL[n];
				}
	   		 }
	   		 	/* NEW line */
	   		else{	
	   		 	*r=y[*C]; 
	   		 	swap(&(nP[np].L[l]),&(nP[np].L[L]));
				nP[++np]=CL[n]; ccf=1;
	    	}
		}
		  
		/* check/complete nP */
		while(++c<_X->nv){   
			int s=S[c]; 
			L=np; 
			ccf=cf;
	    	if(s<c) s=S[s];	
	    	while(L--){	
	    		int j=c; 
	    		C=nP[L].C; 
	    		y=x[nP[L].L[l]];
				while(++j<=s) if(y[C[c]]<y[C[j]]) swap(&C[c],&C[j]);
				if(ccf){   
					Long d=y[C[c]]-r[c];
					/* BAD line */
		    		if(d<0) {if(--np>L) nP[L]=nP[np];}		
		    		/* BETTER line */
		    		else if(d){	
		    			r[c]=y[C[c]]; 
		    			cf=0; 
		    			np=L+1; 
		    			*_ns=n+1;
		    		}	
		    		/* else	; */			      /* EQUAL line */
				}
				else { r[c]=y[C[c]]; ccf=1; }
	    	}
		}
	
		cf=1;
		if(--(*_ns) > n) CL[n]=CL[*_ns]; 		/*  write nP to CL  */
		if(SYM_Nmax < (cf=(*_ns+np))){   
			printf("Need SYM_Nmax > %d !!\n",cf);exit(0);
		}
		for(L=0;L<np;L++) CL[(*_ns)++]=nP[L];
     }
     
     y=x[CL->L[l]];					       /* compute S */
     
     {	
     	int c=0, *C=CL->C;
		while(c<_X->nv){   
			int s=S[c]+1; 
			S[c]=c; 
			while(++c<s){   
				if(y[C[c]]==y[C[c-1]]) ++S[ S[c]=S[c-1] ]; 
	        	else S[c]=c;
	    	}
		}
     }
}

void Aux_vNF_Init(vNF *_X, Long x[][VERT_Nmax], PERM *CL, int *S, int *_ns)
{    
	int i, j, nn; 
	Long *b, *y;
    PERM P, *q, *p;		       /* b=x[nb] -> best;  y=x[nn] -> next */
    
    for(i=0;i<_X->nf;i++) P.L[i]=i;
    for(j=0;j<_X->nv;j++) P.C[j]=j; /* init P */
    
    q=CL; *q=P; b=*x;          /* P=CL[ns-1] StartPerm; maximize 0-th line */
    
    for(j=1;j<_X->nv;j++)if(b[q->C[0]]<b[q->C[j]])swap(&q->C[0],&q->C[j]);
    
    for(i=1;i<_X->nv;i++){  
    	for(j=i+1;j<_X->nv;j++) if(b[q->C[i]]<b[q->C[j]]) 
        swap(&q->C[i],&q->C[j]);
    }

	/* maximize nn-th line */
    for(nn=1;nn<_X->nf;nn++){  
    	Long d; 
    	p=&CL[*_ns]; 
    	*p=P; y=x[nn]; /* nb=*q=*b=best, nn=*p=*y=next */
        {
        	int m=0; 
        	for(j=1;j<_X->nv;j++) if(y[p->C[m]]<y[p->C[j]]) m=j;
            if(m) swap(&p->C[0],&p->C[m]);
        }
        if((d=y[p->C[0]]-b[q->C[0]]) < 0) continue;   /* d<0 => forget this */
        for(i=1;i<_X->nv;i++){   
        	int m=i; 
        	for(j=i+1;j<_X->nv;j++) if(y[p->C[m]]<y[p->C[j]]) m=j; 
            if(m>i) swap(&p->C[i],&p->C[m]);
            if(d==0) if((d=y[p->C[i]]-b[q->C[i]]) <0) break;
    	}
        if(d<0) continue;
        swap(&p->L[0],&p->L[nn]);		 /* p->L[nn]=0; p->L[0]=nn; */ 
        if(d==0) (*_ns)++; 
        else {*q=*p; *_ns=1; b=y;}                    /* d>0 => forget prev */
    }
    
    y=x[CL->L[0]]; 
    S[0]=0; 
    /* compute S */
    for(i=1;i<_X->nv;i++){
    	if(y[CL->C[i]]==y[CL->C[i-1]]) ++S[ S[i]=S[i-1] ]; 
    	else S[i]=i;
    
    }    
}

/* return "X < Y" */
int  Aux_XltY_Poly_NF(Long X[][VERT_Nmax],Long Y[][VERT_Nmax], int *n,int *nv)
{    
	int i, j; Long d;					  
    for(i=0;i<*n;i++)for(j=0;j<*nv;j++) if((d=X[i][j]-Y[i][j])){	
    	if(d<0) return 1; 
    	else return 0;
    }
    return 0;
}


void Aux_Make_Triang(PERM *CL,int ns,Long V[][VERT_Nmax],int*n,int*nv,int *t)
{    
	/* x :: make X :: if X>Y */
	int i, j, s, x=0, g=0, ps=1;		   
    Long X[POLYDIM][VERT_Nmax], Y[POLYDIM][VERT_Nmax];
    
    for(i=0;i<*n;i++) for(j=0;j<*nv;j++) X[i][j]=V[i][CL->C[j]];
    if(!Aux_Make_Poly_NF(X,n,nv)) exit(0); 		  /* t>0: print NFs */
	
	/*  -1: calc CL.s */
    if(*t){
    	if(*t<=0){ 
    		CL->s=1; 
    		if(*t+1){puts("t<-1 in Aux_Make_Triang");exit(0);} 
    	}
    } 
    
    for(s=1;s<ns;s++)CL[s].s=0;

    for(s=1;s<ns;s++)
    if(x){	
    	for(i=0;i<*n;i++) for(j=0;j<*nv;j++) X[i][j]=V[i][CL[s].C[j]];
     	if(!Aux_Make_Poly_NF(X,n,nv)) exit(0); 
		if(Aux_XltY_Poly_NF(X,Y,n,nv)) x=0;

		if(*t){   
	    	if(x==0){
	    		if(*t<0){ 
	    			int k; 
	    			for(k=g;k<s;k++) CL[k].s=0; 
	    			CL[s].s=1;*t=-1;
	    		}
				g=s; ps=1; 
	   		 } 
	    	else if(!Aux_XltY_Poly_NF(Y,X,n,nv)){ 
	    		if(*t<0) {CL[s].s=1;(*t)--;} 
	    		ps++;
	    	}
        }
    }
    else{	
    	for(i=0;i<*n;i++) for(j=0;j<*nv;j++) Y[i][j]=V[i][CL[s].C[j]];
     	if(!Aux_Make_Poly_NF(Y,n,nv)) exit(0); 
		if(Aux_XltY_Poly_NF(Y,X,n,nv)) x=1;

		if(*t){   
	    	if(x==1){	
	    		if(*t<0){ 
	    			int k; 
	    			for(k=g;k<s;k++) CL[k].s=0; 
	    			CL[s].s=1;*t=-1;
	    		}
				g=s; ps=1;
	    } 
	    else if(!Aux_XltY_Poly_NF(X,Y,n,nv)){ 
	    	if(*t<0) {CL[s].s=1;(*t)--;} 
	    	ps++;
	    }
	}
    }
    
    if(*t>0)
    printf("\nPoly NF:  NormalForm=try[%d]  #Sym(VPM)=%d  #Sym(Poly)=%d\n",g,ns,ps);
    if(x) for(i=0;i<*n;i++)for(j=0;j<*nv;j++) V[i][j]=Y[i][j];
    else  for(i=0;i<*n;i++)for(j=0;j<*nv;j++) V[i][j]=X[i][j];
}


void Make_VPM_NF(int *v, int *f, Long x[VERT_Nmax][VERT_Nmax],      /* in */
		PERM *CL,int *ns,Long VPM_NF[VERT_Nmax][VERT_Nmax])  /* out */
{    
	int i, j;
	int S[VERT_Nmax]; 
	int nsF=0, nsM=0; 

	/* X=VPM */
    volatile vNF auX; 
    vNF *_X= (vNF*) &auX; 
    
    /* define the number of vertices and number of equations */ 
    _X->nv=*v;
    _X->nf=*f;
     
    *ns=1; 
    Aux_vNF_Init(_X, x, CL, S, ns);             /* init = 1st line */
     
    for(i=1;i<_X->nf-1;i++){
    	Aux_vNF_Line(i,_X,x,CL,S,ns);  /* lines of NF */
		#ifdef	WARN_BIG_NS
			if((WARN_BIG_NS<=(*ns))||nsF) nsF=1;
		#endif
		if(*ns>nsM) nsM=*ns; 
	}
	
    _X->ns=*ns; 
    
     /* write VPM-NF to _X */
    for(i=0;i<_X->nv;i++){  
    	for(j=0;j<_X->nf;j++) /* _X->x */ VPM_NF[j][i]=x[CL->L[j]][CL->C[i]];
    }
    
    if(nsF)printf("WARNing: ns_max=%d -> ns=%d\n",nsM,*ns);
}


void Aux_pNF_from_vNF(PERM *CL,int *ns,int *v,int *d,
		Long VM[POLYDIM][VERT_Nmax],			      /* in */
		Long pNF[POLYDIM][VERT_Nmax],int *t)		     /* out */
{    
	int i,j;
    for(i=0;i<*d;i++) for(j=0;j<*v;j++) pNF[i][j]=VM[i][j];
    Aux_Make_Triang(CL,*ns,pNF,d,v,t);
}


int  Make_Poly_Sym_NF(PolyPointList *_P, VertexNumList *_V, EqList *_F, 
		      int *SymNum, int V_perm[][VERT_Nmax], 
		      Long NF[POLYDIM][VERT_Nmax])
{    
	int i, j;
	int ns;
	int t=-1;
	
	/* define the number of points, vertices and equations */
	int *d=&_P->n;
	int *v=&_V->nv;
	int *f=&_F->ne;
	
	int *C; 
	
    PERM *CL = (PERM *) malloc ( sizeof(PERM) *(SYM_Nmax+1));
    
    /* define the vertex matrix, vertex pairing matrix and the normal form vertex pairing matrix */
    Long VM[POLYDIM][VERT_Nmax];
    Long VPM[VERT_Nmax][VERT_Nmax];
    Long VPM_NF[VERT_Nmax][VERT_Nmax];

	/* initialise vertex pairing matrix */
    Init_rVM_VPM(_P,_V,_F,d,v,f,VM,VPM);

	/* compute the normal form of the vertex pairing matrix */
    Make_VPM_NF(v,f,VPM,CL,&ns,VPM_NF);
    
    /* order the normal form */
    New_pNF_Order(v,f,CL,&ns,VPM_NF);
    
    Aux_pNF_from_vNF(CL,&ns,v,d,VM,NF,&t);
    
    free(CL); 
    
    return ns;
}


void SL_swap(SL_Long *X, SL_Long *Y)
{    
	SL_Long A=*Y; *Y=*X; *X=A; 
}


SL_Long SL_Egcd(SL_Long A0, SL_Long A1, SL_Long *Vout0, SL_Long *Vout1)  
{    
	register SL_Long V0=A0, V1=A1, A2, X0=1, X1=0, X2=0;
    while((A2 = A0 % A1)) { X2=X0-X1*(A0/A1); A0=A1; A1=A2; X0=X1; X1=X2; }
    *Vout0=X1, *Vout1=(A1-(V0) * X1)/ (V1); 
    return A1;
}


int  SL2Z_Make_Poly_NF(Long X[][VERT_Nmax], int *n, int *nv, SL_Long S[POLYDIM][POLYDIM])
{
	SL_Long NF[POLYDIM][VERT_Nmax], a, b, g;
    int i, j, C=-1, L; 
    
    for(i=0;i<*n;i++)for(j=0;j<*n;j++)S[i][j]=(i==j);		  /* init S */
    
    for(i=0;i<*n;i++)for(j=0;j<*nv;j++)NF[i][j]=0;		
    
    for(L=0;L<*n-1;L++)
     {  int N=0;
        while(!N)                /* find column C with N>0 non-zero entries */
        {   ++C; for(i=0;i<*n;i++) 
	    {	for(j=0;j<*n;j++) NF[i][C]+=S[i][j]*X[j][C];
	    }
            if(NF[i=L][C]) N++; 
	    while(++i<*n) if(NF[i][C])
	    {	N++;					 /* make NF[i][C]=0 */
		if(NF[L][C])
		{   SL_Long A; 
		    g=SL_Egcd(NF[L][C],NF[i][C],&a,&b);
		    for(j=0;j<*n;j++)			/* c=-N_iC/g */
		    {	A=a*S[L][j]+b*S[i][j];		/* d= N_LC/g */
			S[i][j]=(NF[L][C]/g)*S[i][j]-(NF[i][C]/g)*S[L][j];
			S[L][j]=A;
		    }
		    NF[L][C]=g; NF[i][C]=0;
		}
		else
		{   SL_swap(&NF[L][C],&NF[i][C]);
		    for(j=0;j<*n;j++)SL_swap(&S[L][j],&S[i][j]);
		}
	    }
	    if(NF[L][C]<0) 
	    {	NF[L][C]*=-1;for(j=0;j<*n;j++)S[L][j]*=-1;	    /* sign */
	    }
	    if(N) for(i=0;i<L;i++)	 /* make upper diag minimal nonneg. */
	    {	SL_Long R=NF[i][C]/NF[L][C];
		if((NF[i][C]-R*NF[L][C])<0) R-=1;
		NF[i][C]-=R*NF[L][C]; for(j=0;j<*n;j++)S[i][j]-=R*S[L][j];
	    }
	}
     }
     while(++C<*nv)
     {	for(i=0;i<*n;i++) for(j=0;j<*n;j++) NF[i][C]+=S[i][j]*X[j][C];
	if(L)if(NF[L][C])
	{   if(NF[L][C]<0)
	    {	NF[*n-1][C]*=-1; for(j=0;j<*n;j++) S[L][j]*=-1;
	    }
	    for(i=0;i<L;i++)		 /* make upper diag minimal nonneg. */
	    {	SL_Long R=NF[i][C]/NF[L][C];
		if((NF[i][C]-R*NF[L][C])<0) R-=1;
		NF[i][C]-=R*NF[L][C]; for(j=0;j<*n;j++)S[i][j]-=R*S[L][j];
	    }
	    L=0;
	}
     }
     for(i=0;i<*n;i++)for(j=0;j<*nv;j++)
     if( labs(X[i][j]=NF[i][j]) > NFX_Limit)
     {	fprintf(stderr,"NFX_Limit in SL: I need %ld !!\n",
		labs(X[i][j]));return 0;
     }
     return 1;
}

