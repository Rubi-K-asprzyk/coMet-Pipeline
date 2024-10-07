#include <math.h>
#include <stdio.h>
// #include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
//#include "xatools.h"
//#include "SimCom.h"
//#include "xatools.h"

#define NR_END 1
#define FREE_ARG char* 
#ifndef ERRORFILE
	#define ERRORFILE "error.txt"
#endif
#ifndef MISSING
	#define MISSING -999.99F
#endif

#ifndef MISSVAL
	#define MISSVAL HUGE_VAL
#endif


/*********************************************************************/

int *ivector(long int nl, long int nh)
/*allocate a int vector with subscript range[nl..nh]*/
{
	int *v;
	
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in ivector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_ivector(int *v, long int nl, long int nh)
/*free a float vector allocated with ivector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/*********************************************************************/

long *lvector(long int nl, long int nh)
/*allocate a int vector with subscript range[nl..nh]*/
{
	long *v;
	
	v=(long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in lvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_lvector(long *v, long int nl, long int nh)
/*free a float vector allocated with ivector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/*********************************************************************/

float *vector(long int nl, long int nh)
/*allocate a float vector with subscript range[nl..nh]*/
{
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in vector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_vector(float *v, long int nl, long int nh)
/*free a float vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/********************************************************************/

double *dvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_dvector(double *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}
/********************************************************************/

char *cvector(long int nl, long int nh)
/*allocate a double vector with subscript range[nl..nh]*/
{
	char *v;

	v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char))); 
	if (!v) {
		write(ERRORFILE,"\nallocation failure in dvector()");
		exit(1);
	}
	return v-nl+NR_END;
}

void free_cvector(char *v, long int nl, long int nh)
/*free a double vector allocated with vector()*/
{
	free((FREE_ARG) (v+nl-NR_END));
}

/****************************************************************************/
int **imatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/*allocate pointers to rows*/
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in imatrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in imatrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_imatrix(int **m, long nrl,long nrh,long ncl, long nch)
/*free an int matrix allocated by imatrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

float **matrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/*allocate pointers to rows*/
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float *)));

	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_matrix(float **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

double **dmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/*allocate pointers to rows*/
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double *)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_dmatrix(double **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

char **cmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;
	
	/*allocate pointers to rows*/
	m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char *)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_cmatrix(char **m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	
/****************************************************************************/

float ***fpmatrix(long nrl, long nrh,long ncl, long nch)
/*allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float ***m;
	
	/*allocate pointers to rows*/
	m=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float **)));
	if (!m) {
		write(ERRORFILE,"\nallocation failure 1 in matrix()");
		exit(1);
	}
	m += NR_END;
	m -= nrl;
	
	/*allocate rows and set pointers to them*/
	m[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!m[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in matrix()");
		exit(1);
	}
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/*return pointer to array of pointers to rows*/
	return m;
}
	 
void free_fpmatrix(float ***m, long nrl,long nrh,long ncl, long nch)
/*free an float matrix allocated by matrix()*/
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}	

/****************************************************************************/

float ***f3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a float 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_f3tensor(float ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an float f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	


/****************************************************************************/

double ***d3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a float 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_d3tensor(double ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an float f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	

/****************************************************************************/

int ***i3tensor(long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*allocate a int 3tensor with subscript range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j, nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	int ***t;
	
	/*allocate pointers to pointers to rows*/
	t=(int ***) malloc((size_t)((nrow+NR_END)*sizeof(int **)));
	if (!t) {
		write(ERRORFILE,"\nallocation failure 1 in f3tensor()");
		exit(1);
	}
	t += NR_END;
	t -= nrl;
	
	/*allocate pointers to rows and set pointers to them*/
	t[nrl]=(int **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int *)));
	if (!t[nrl]) {
		write(ERRORFILE,"\nallocation failure 2 in f3tensor()");
		exit(1);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/*allocate rows and set pointers to them*/
	t[nrl][ncl]=(int *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(int )));
	if (!t[nrl][ncl]) {
		write(ERRORFILE,"\nallocation failure 3 in f3tensor()");
		exit(1);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*return pointer to array of pointers to rows*/
	return t;
}

void free_i3tensor(int ***t,long nrl, long nrh,long ncl, long nch,long ndl,long ndh)
/*free an int f3tensor allocated by f3tensor()*/
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}	

/****************************************************************************/

float bico(int n, int k)
/*returns the binomial coefficient Cn(haut)k(bas) = n!/(k!(n-k)!) as a floating-
point number*/
{
	float factln(int n);
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k))); 
}

/****************************************************************************/

float factln(int n)
/*returns ln(n!)*/
{
	float gammln(float xx);
	static float a[101];
	if (n<0) {
		write(ERRORFILE,"\nNegative factorial in routine factln"); 
		exit(1);
	}
	if (n <=1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}

/*******************************************************/

float gammln(float xx)
/*numerical recipes p.214, returns the value of ln(gamma(xx)) for xx>0*/
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677, 
		24.01409824083091,-1.231739572450155, 
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y; 
	return -tmp+log(2.5066282746310005*ser/x);
}

/*******************************************************/

float factrl(int n)
/*returns the value n! as a floating-point number*/
{
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;
	
	if(n<0) {
		write(ERRORFILE,"\nNegative factorial in routine factrl");
		exit(1);
	}
	if(n>32) return exp(gammln(n+1.0));
	
	while(ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}


/******************************************************************************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
		int j;
		long k;
		static long iy=0;
		static long iv[NTAB];
		float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1; 
		else *idum = -(*idum); 
		for (j=NTAB+7;j>=0;j--) { 
			k=(*idum)/IQ; 
			*idum=IA*(*idum-k*IQ)-IR*k; 
			if (*idum < 0) *idum += IM; 
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j= (int) iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=(float) AM*iy) > RNMX) return (float) RNMX;
	else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/********************************************************************************/

/* (C) Copr. 1986-92 Numerical Recipes Software . */
long timeseed(void)
/*output a seed for random generator based on system time by  multiplying sec*min*hours or give 1 if 
sec=min=hour=0*/
{
		struct tm *local;
		time_t t;
		long seed;

		seed=1;
		t = time(NULL);
		local = localtime(&t); 
		if(local->tm_sec) seed=local->tm_sec; 
		if(local->tm_min) seed *= local->tm_min; 
		if(local->tm_hour) seed *= local->tm_hour;
		return seed;
}		 

void main(int *rNS, int *rNip, int *rNG, int *rW, int *rL, int X[], int habp[], double *rS, double rFs[], int rdisp[], int *rsamehab, int *rdiffhab, int *rHabitatP, double rMatSelEsp[], int *rNhabitat , double rSelecExport[], double *rMs, double rMatTraitsEsp[], int *rTraits, double *rCompetition, double *rMAXFunc, int *rfailedRec, int *rsuccesRec, double RfEcumR[], double *rMno, double rfEcumPar[], int RdispRec[], double *rsumProbRecAll, double *rsumProbRecSucc, double *rsumProbaRecFail, int *rNgen, double rMatProbRec[], int rMatCountRec[], double rInteractionsMat[], double *requalizing,int Xhsg[], double *rstoc, double *rtest, int *rtest2, int *rerror) // NEW VARIABLES  
{                                             
	int NS,Nip, NG, W, L, NP,*NIp,s,x,y,dx,dy,xo,yo,p,i,**nisp,**nispN,po,**sps,sp,repeat;//total number of species, plots, individuals per plot, #ind of sp s in plot p, #ind of sp s in plot p for the new generation      Objects no more used: ,Xp,Yp
	double noM,Mr,Mn,Mp,Ms,selection;	//migration rates from random plots, neighboring plots, source community, within local patch         Objet supprim\E9 : ,alphadisp,deltadisp,betadisp
	double **cFsp,*Fs;			//frequency of species s in source community, cumulated freq of species s in plot p (or source community if p=-1)

	int gen,rep,*habitatp,*disp,*dispRec;			//#generations, # time steps                   Objet supprim\E9 :    HabitatSize,   habitat,  NR, *selhab,
	float randN,naf,df,*fEcum,fcum,angle,dmax, *fEcumPar;		//random number
	
	int *NSobs,samehab=0, diffhab=0, failedRec=0, succesRec=0, spO;
	int ***nihsg, h, g; // total number of ind of species s per habitat type at generation g

	double **matselhab, **traitsesp, competition, probrec, eucli, eucliF,MAXFunc, Mno, sumProbRecSucc=0., sumProbRecAll=0., sumProbaRecFail=0.;//, *SelEsp;  // matselhab = selection coefficient for each species within each habitat
	int Nhabitat,q, traits, k, *vecCountRec, NbeGen=0, Ngen, **MatCountRec, **InteractionsMatPosi;  // # of different habitat,
	double **MatProbRec,*vecProbRec, **InteractionsMat;
	long seed;
	double ProbRecFilter, ProbRecEqu, ProbRecStab, equalizing, stoc; 
	int error=0;

	NS=(*rNS);
	Nip=(*rNip);
	NG=(*rNG); #total number of mortality/replacement events
	W=(*rW);
	L=(*rL);
	selection=(*rS);
	Ngen=(*rNgen); #total number of generations where 1 generation = mortality/replacement events equal to toal number of individuals in the community

	equalizing=(*requalizing);
	stoc=(*rstoc);

	Fs=&rFs[0];

	NP=W*L;	//# of plots (=subcommunities)
	Nhabitat=(*rNhabitat);  //# of different habitat
	traits=(*rTraits);  //# of traits

	NIp=ivector(0,NP);	// # of individuals per plot 
	for(p=1;p<=NP;p++) NIp[p]=Nip;
	NIp[0]=0; 
	for(p=1;p<=NP;p++) NIp[0]+=NIp[p];

	habitatp=ivector(0,NP);		//habitat of each plot

	for (p=1;p<=NP;p++) {     //habitat of each plot importe de R
		habitatp[p]=rHabitatP[p-1];
		if(habitatp[p]<1 || habitatp[p]>Nhabitat) error = 1;  //check range of habitat name
	}

	competition=(*rCompetition);
	MAXFunc=(*rMAXFunc);



	Mr=0.;		//inter-plot 'random' migration rate
	Mn=0.1;		//inter-plot migration rate between immediate neighbouring plots
	Mp=0.;		//inter-plot migration rate within subsets of plots
	Ms=*rMs;	//migration rate from 'source' metacommunity
	Mno=*rMno;
	noM=1.-Ms-Mn-Mr-Mp; //probability than a recruited individual is not an immigrant

	//random number generators initialisation
	seed=timeseed()*(-1);

	//initialize
	nisp=imatrix(-1,NS,-1,NP);	//# of individuals per species and per plot
	nispN=imatrix(0,NS,0,NP);	//

	// NEW VARIABLE 
	nihsg=i3tensor(0,Nhabitat,1,NS,0,Ngen);	//in habitat h, number of ind per species s and generation g

	cFsp=dmatrix(0,NS,-1,NP);	//cumulated species freq per plot
	NSobs=ivector(-1,NP);		//# of species present per plot
	sps=imatrix(-1,NP,0,NS);	//correspondance between global species number and local (per plot) species number  

	matselhab=dmatrix(0,NS,0,Nhabitat);  
	traitsesp=dmatrix(0,NS,0,traits);                 

	InteractionsMat=dmatrix(0,NS*NS,0,Nhabitat);   
	InteractionsMatPosi=imatrix(0,NS,0,NS);

	//Fill matrix of affinity of each species for each habitat:  MatSelHab
	q=0;                                 
	for (i=1;i<=Nhabitat;i++) {
		for (p=1;p<=NS;p++) {
			matselhab[p][i]=rMatSelEsp[q];                           
			q++;
		} 
	}
	//Fill matrix of traits of each species   traitsesp
	q=0;
	for (i=1;i<=traits;i++) {
		for (p=1;p<=NS;p++) { 
			traitsesp[p][i]=rMatTraitsEsp[q];                           
			q++;
		} 
	}
	//Fill matrix of species-species interactions for each habitat: InteractionsMat
	q=0;
	for (i=1; i<=Nhabitat;i++) {
		for (p=1;p<=NS*NS;p++) {
			InteractionsMat[p][i]=rInteractionsMat[q];
			q++;        
		}
	}
	
	// Remplissage de la matrice InteractionsMatPosi qui indique la position des valeurs dans la matrice InteractionsMat
	q=1;
	for (i=1;i<=NS;i++){
		for (p=1;p<=NS;p++){
			InteractionsMatPosi[p][i]=q;
			q++;
		}
	}

	//cumulated species abundances in source metacommunity
	cFsp[0][-1]=0.;
	sp=0;
	for(s=1;s<=NS;s++)if(Fs[s-1]){
		sp++;
		cFsp[sp][-1]=cFsp[sp-1][-1]+Fs[s-1];
		sps[-1][sp]=s;
	}
	NSobs[-1]=sp;//total number of species (with non null freq)

	//COMPUTE PROBA DISTRIB OF POWER-EXPON DISPERSAL FUNCTION
 /*
	alphadisp=deltadisp*exp(gammln(2./betadisp))/exp(gammln(3./betadisp));
	dmax=(int)(sqrt(L*L+W*W)/2.);
	fEcum=vector(0,(10*dmax)+2);
	fcum=0.;
	for(df=0.1;df<=dmax;df+=0.1){
		fcum+=df*exp(-1.*pow(df/alphadisp,betadisp));
		fEcum[(int)(10*df)]=fcum;
	}
	for(df=0.1;df<=dmax;df+=0.1)fEcum[(int)(10*df)]/=fcum;
*/ 
	dmax=(int)(sqrt(L*L+W*W)/2.);
	fEcum=vector(0,(10*dmax)+2);
	i=0;
	for (df=0.1;df<=dmax;df+=0.1){
		fEcum[(int)(10*df)]=RfEcumR[i];
		i++; 
	}
	
	fEcumPar=vector(0,(10*dmax-5)+2);
	i=0;
	for (df=0.1;df<=dmax-1;df+=0.1){
		fEcumPar[(int)(10*df)]=rfEcumPar[i];
		i++; 
	}

	//Initialisation de disp
	
	disp=ivector(0,dmax);
	for(p=0;p<=dmax; p++) disp[p]=0;

	dispRec=ivector(0,dmax);
	for (p=0;p<=dmax; p++) dispRec[p]=0;

	// Initialisation de vecCountRec  et   MatProbRec
	vecCountRec=ivector(0,9);
	for (p=0;p<=9;p++) vecCountRec[p]=0;  
	vecProbRec=dvector(0,9);
	for (p=0;p<=9;p++) vecProbRec[p]=0.;  
	MatProbRec=dmatrix(0,Ngen,0,9);
	//for (p=0;p<=Ngen;p++)for(s=0;s<=10;s++) MatProbRec[p][s]=0.;
	MatCountRec=imatrix(0,Ngen,0,9);


	//CREATE INITIAL GENERATION
	for(p=0;p<=NP;p++)for(s=0;s<=NS;s++) nisp[s][p]=0;
	for(p=1;p<=NP;p++){
		for(i=1;i<=NIp[p];i++){
			randN=ran1(&seed);
			for(sp=1;sp<=NSobs[-1];sp++) if(cFsp[sp][-1]>randN) break;
			nisp[sps[-1][sp]][p]++;
		}
		sp=0;
		for(s=1;s<=NS;s++)if(nisp[s][p]){
			sp++;
			sps[p][sp]=s;
		}
		NSobs[p]=sp;
	}

	//sum over plots per species
	for(s=1;s<=NS;s++) nisp[s][0]=0;
	for(s=1;s<=NS;s++)for(p=1;p<=NP;p++) nisp[s][0]+=nisp[s][p];
	sp=0;
	for(s=1;s<=NS;s++)if(nisp[s][0]){
		sp++;
		sps[0][sp]=s;
	}
	NSobs[0]=sp; 

	//sum over species per plot
	for(p=0;p<=NP;p++)for(s=1;s<=NS;s++) nisp[0][p]+=nisp[s][p];

	for(p=0;p<=NP;p++)for(s=0;s<=NS;s++) nispN[s][p]=nisp[s][p];    // nispN will be the matrix where one individual will be removed  // Sampling One individual
	
// NEW VARIABLE		
	for(h=0;h<=Nhabitat;h++)for(s=1;s<=NS;s++)for(g=0;g<=Ngen;g++) nihsg[h][s][g] = 0; // initialize tensor for recording species abundance accross habitats and generations
	for(p=1;p<=NP;p++) for(s=1;s<=NS;s++) nihsg[habitatp[p]][s][0] += nisp[s][p];   // species abundances per habitat for generation 0 (initial community)
// NEW VARIABLE

	rep=0; 
	//LOOP OVER GENERATIONS
	for(gen=1;gen<=NG;gen++){   

		randN=ran1(&seed);         			                                           // Sampling One individual
		p=1+(int)(randN*NP);        //  choosing randomly a plot                    // Sampling One individual
		
		randN=ran1(&seed);                                                          // Sampling One individual
		spO=(int)(randN*NSobs[p]);                                                 // Sampling One individual
		while(cFsp[spO][p]>randN && spO>0) spO--;                                 // Sampling One individual
		while(cFsp[spO][p]<randN && spO<NSobs[p]) spO++;                           // Sampling One individual
//    for(spO=1;spO<=NSobs[p];spO++) if(cFsp[spO][p]>randN) break;        //  choosing randomly a species        // Sampling One individual
	 // spO=1+(int) (randN*NSobs[p]);                                                     // Sampling One individual
		nispN[sps[p][spO]][p]--;          //  remove one individual of that species in the chosen plot    // Sampling One individual

		//CREATE NEW GENERATION
//	for(p=0;p<=NP;p++)for(s=0;s<=NS;s++) nispN[s][p]=0;      // removing All individuals
		repeat=1;                // Sampling One individual
		while(repeat==1){       // Sampling One individual
	//	for(p=1;p<=NP;p++){            // removing All individuals
	//	for(i=1;i<=NIp[p];i++){        // removing All individuals
			randN=ran1(&seed);
			
			if(randN<Ms) po=-1; //migrant from source community
			else if (ran1(&seed)<Mno) {         // 
				po=p;  // no migration, recruitment within the same plot // NEW          	
				disp[0]++;
			}
			
			else{	//if there is a migration OR if the rate of within-plot recruitment is null, migrant according to power expon function
				naf=ran1(&seed);                                        
				if(Mno==0.) {
					df=0.0;
					do{
						df+=0.1;
					}while(fEcum[(int)(10*df)]<naf);
				}
				else{
					df=1.0;
					do{
						df+=0.1;
					}while(fEcumPar[(int)(10*df)-5]<naf);
				}
				//choose dispersal angle
				angle=ran1(&seed)*6.283185;
				//compute dx and dy and round to the closest integer
				if(cos(angle)>0) dx=(int)(df*cos(angle)+0.5);
				else dx=(int)(df*cos(angle)-0.5);
				if(sin(angle)>0) dy=(int)(df*sin(angle)+0.5);
				else dy=(int)(df*sin(angle)-0.5);

				x=1+(p-1)%L;		
				y=1+(int)((p-1)/L);
				xo=x+dx;
				if(xo<1) xo=xo+L;
				if(xo>L) xo=xo-L;
				yo=y+dy;
				if(yo<1) yo=yo+W;
				if(yo>W) yo=yo-W;
				po=L*(yo-1)+xo;
				//dx = min(abs(x-xo),L-abs(x-xo));
				//dy = min(abs(y-yo),W-abs(y-yo));
				if((int)(sqrt(dx*dx+dy*dy)+0.5) > dmax) error = 2;
				disp[(int)(sqrt(dx*dx+dy*dy)+0.5)]++;
				if(habitatp[p]==habitatp[po]) samehab++;
				else diffhab++; 
			}

			//define species
			randN=ran1(&seed);
	//	  for(sp=1;sp<=NSobs[po];sp++) if(cFsp[sp][po]>randN) break;
			sp=(int)(randN*NSobs[po]);
			while(cFsp[sp][po]>randN && sp>0) sp--;
			while(cFsp[sp][po]<randN && sp<NSobs[po]) sp++;       

	//  Apply selection and competition
			repeat=0;
			if(selection==0.) ProbRecFilter=1.;
			else {
				ProbRecFilter=matselhab[sps[po][sp]][habitatp[p]]; // recruitment proba is fixed before simulation and stored in matselhab matrix
			}

			if (equalizing==0.) ProbRecEqu=1.;  // --> Equalizing mechanisms
			else{
				ProbRecEqu=0.;
				for(s=1;s<=NSobs[p];s++) ProbRecEqu+=InteractionsMat[InteractionsMatPosi[sps[po][sp]][sps[p][s]]][habitatp[p]]*(double)(nispN[sps[p][s]][p])/(double)(NIp[p]-1); //  For all species within the plot -->  Recruitment proba determined by fitness comparison between recruited species [sps[po][sp]] and all other present species [sps[p][s]] 
			}

			if(competition==0.) ProbRecStab=1.;
			else{   // if competition, species composition of the sink plot (p) reduce the recruitment proba --> Stabilizing mechanisms
				eucliF=0.;
				ProbRecStab=0.;
				for (s=1;s<=NSobs[p];s++) {    // For all species within the plot
					eucli=0.;
					for (k=1;k<=traits;k++) eucli+=(traitsesp[sps[p][s]][k]-traitsesp[sps[po][sp]][k])*(traitsesp[sps[p][s]][k]-traitsesp[sps[po][sp]][k]);    // square of sum of difference for all traits between the immigrant and the species present in plot ###
					eucliF+=sqrt(eucli)*nispN[sps[p][s]][p]; // sum of all distances weighted by the abundance of species s     //  eucliF is the sum of the functional distances between the immigrant and the species composition of the sink plot (p)                                                             
				} // end s
				eucliF/=NIp[p];  // divided by the nbe of individuals in the plot --> mean functionnal distance of the immigrant species with all other species (weighted by the abundance of species) present in the plot sink 
				ProbRecStab=eucliF/MAXFunc;    // Recruitment probability is standardized by the maximal functional distance possible                   
		
			} // end competition
			
			probrec=stoc + competition*ProbRecStab + equalizing*ProbRecEqu + selection*ProbRecFilter; // 
			(*rtest2)=habitatp[p];
			
			
			if (probrec==0.) repeat=1;    //  si la probabilit\E9 de recrutement est \E9gal \E0 0
			else {
				if(ran1(&seed)>probrec) repeat=1;     
				else{
					repeat=0;
				}
			}
			sumProbRecAll+=probrec;                        // the sum of recruitment probability   
			if (repeat==1 && po==-1) {
				vecCountRec[0]++;    // If it is a migrant from metacommunity that failed to install
				vecProbRec[0]+=probrec;
															}
			if (repeat!=1 && po==-1) {
				vecCountRec[1]++;    // If it is a migrant from metacommunity that succesfully installed
				vecProbRec[1]+=probrec;
																}
			if (repeat==1 && p==po) {
				vecCountRec[2]++;    // If recruitment failed and no migration
				vecProbRec[2]+=probrec;
															}
			if (repeat!=1 && p==po) {
				vecCountRec[3]++;    // If recruitment succesful and no migration
				vecProbRec[3]+=probrec;
															} 
			if (repeat==1 && p!=po && po!=-1) {
				vecCountRec[4]++;    // If recruitment failed and immigration from the community
				vecProbRec[4]+=probrec;
															} 
			if (repeat!=1 && p!=po && po!=-1) {
				vecCountRec[5]++;    // If recruitment succesful and immigration from the community
				vecProbRec[5]+=probrec;
															} 
			if (repeat==1 && p!=po && habitatp[p]==habitatp[po] && po!=-1) {
				vecCountRec[6]++;  //  If recruitment failed and immigration from the community from the same habitat
				vecProbRec[6]+=probrec;
															}         
			if (repeat!=1 && p!=po && habitatp[p]==habitatp[po] && po!=-1) {
				vecCountRec[7]++;  //  If recruitment succesful and immigration from the community from the same habitat
				vecProbRec[7]+=probrec;
															}         
			if (repeat==1 && p!=po && habitatp[p]!=habitatp[po] && po!=-1) {
				vecCountRec[8]++;  //  If recruitment failed and immigration from the community from a different habitat
				vecProbRec[8]+=probrec;
															}         
			if (repeat!=1 && p!=po && habitatp[p]!=habitatp[po] && po!=-1) {
				vecCountRec[9]++;  //  If recruitment succesful and immigration from the community from a different habitat                   
				vecProbRec[9]+=probrec;
															}         
			if (repeat==1) sumProbaRecFail+=probrec;       // the sum of recruitment probability when recruitment failed
			if (repeat!=1) sumProbRecSucc+=probrec;        // the sum of recruitment probability when recruitment is successful
			if (repeat==1) failedRec++; 
			if (repeat!=1) succesRec++;					
			if (repeat!=1 && p==po) dispRec[0]++;
			if (repeat!=1 && p!=po) dispRec[(int)(sqrt(dx*dx+dy*dy)+0.5)]++;
			if (repeat!=1) nispN[sps[po][sp]][p]++;	// replace the individual removed by the individual of the chosen species 'sp' coming from plot 'po'	 // Sampling One individual		 
		}//end loop while to recruit a new individual
					//}//end loop p
	

		if (sps[po][sp]!=sps[p][spO])  {     // if the removed individual is a different species of the arriving individual
		//sum over species per plot
			nispN[0][p]=0;
			for(s=1;s<=NS;s++) nispN[0][p]+=nispN[s][p];
		}

		//compute cumulated frequencies per plot
		if (sps[po][sp]!=sps[p][spO])  {
			cFsp[0][p]=0.;
			sp=0;
			for(s=1;s<=NS;s++)if(nispN[s][p]){
				sp++;
				cFsp[sp][p]=cFsp[sp-1][p]+(double)nispN[s][p]/(double)NIp[p];
				sps[p][sp]=s;
			}
			NSobs[p]=sp;
		}
		
		rep++;
		if (rep==L*W*Nip) {      // if rep reach the total number of individual in the community
			rep=0;
			if(NbeGen>Ngen) error = 3;
			for (p=0;p<=9;p++) MatProbRec[NbeGen][p]=vecProbRec[p];        // vecProbRec[p]/(double)vecCountRec[p]; 
			for (p=0;p<=9;p++) MatCountRec[NbeGen][p]=vecCountRec[p];
			for (p=0;p<=9;p++) vecCountRec[p]=0;  
			for (p=0;p<=9;p++) vecProbRec[p]=0.;  
			NbeGen++;
		
		// NEW VARIABLES
		for(p=1;p<=NP;p++) for(s=1;s<=NS;s++) nihsg[habitatp[p]][s][NbeGen] += nispN[s][p];  // species abundances per habitat for generation NbeGen	
		}
	
	}//END loop gen

		//realized dispersal distrib
		/*sum=0.;
		for(dx=0;dx<=dmax;dx++) sum+=disp[dx];
		for(dx=0;dx<=dmax;dx++) disp[dx]/=sum;
	  */
		
		//write vectors for R output
    i=0;     // Exportation du nbe de d'individu par esp\E8ce a la fin des simulations
		for(s=1;s<=NS;s++) {
			for(p=1;p<=NP;p++) { 
			  X[i]=nispN[s][p];
			  i++;
			}
		}

//NEW VARIABLE	
		i=0;     // Exportation du nbe de d'individus par habitat, espece et generation 
		for(h=1;h<=Nhabitat;h++) {
			for(s=1;s<=NS;s++) {
				for(g=0;g<=Ngen;g++) { 
					Xhsg[i]=nihsg[h][s][g];
					i++;
				}
			}
		}

		i=0;          // Exportation du type d'habitat pour chaque sous-communaut\E9s
		for(p=1;p<=NP;p++) {
			habp[i]=habitatp[p];
			i++;
		}
		i=0;         //  realized dispersal distribution
		for(p=0;p<=dmax;p++) {
			rdisp[i]=disp[p];
			i++;
		}
    i=0;         //  realized dispersal distribution after selection/competition
		for(p=0;p<=dmax;p++) {
			RdispRec[i]=dispRec[p];
			i++;
		}
		/* 
   
    i=0;
    for(p=0;p<=9;p++) {
    rvecCountRec[i]=vecCountRec[p];
    i++;
    }
    
    i=0;
    for(p=0;p<=9;p++) {
    rvecProbRec[i]=vecProbRec[p];
    i++;
    }
   */ 
    
		i=0;     // Exportation des sommes de probabilit\E9s de recrutement \E0 chaque g\E9n\E9ration
		for(s=0;s<=(Ngen-1);s++) {
			for(p=0;p<=9;p++) {
			  rMatProbRec[i]=MatProbRec[s][p];
			  i++;
			}
		}
    
		i=0;     // Exportation des nombres de recrutement r\E9ussis et rat\E9s \E0 chaque g\E9n\E9ration
		for(s=0;s<=(Ngen-1);s++) {
			for(p=0;p<=9;p++) {
				rMatCountRec[i]=MatCountRec[s][p];
				i++;
			}
		}
    
    
    (*rtest)=ProbRecEqu;

     
		(*rsamehab)=samehab;    //  Number of dispersal item which were made between plot of same habitat (relation between habitat aggregation and dispersal abilities) 
		(*rdiffhab)=diffhab;    //  Number of dispersal item which were made between plot of different habitat (relation between habitat aggregation and dispersal abilities) 
    (*rfailedRec)=failedRec;
    (*rsuccesRec)=succesRec;
    (*rsumProbRecAll)=sumProbRecAll;
    (*rsumProbRecSucc)=sumProbRecSucc;
    (*rsumProbaRecFail)=sumProbaRecFail;
		(*rerror)=error;
    return; 
}

        /*Rprintf("nisp=%i",nisp[0][6]); 
	     Rprintf("randN=%f",randN);
       R_FlushConsole();
       R_ProcessEvents();
       return;    */
