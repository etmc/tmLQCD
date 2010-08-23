/*************************************************************
 *Written by Andreas Nube (C) 2007
 *
 *This program computes and outputs the roots of a chebycheff 
 *
 *polynomial approximating 1/x One Over X -> oox 
 *************************************************************/




#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <math.h>

#include <time.h>


#include <string.h>
#include <ctype.h>

#include <unistd.h>

#ifdef WITHGALIB
#include "oox_gawrapper.h"
#endif


#ifndef NULL
#define NULL 0
#endif

#define COMPLEX double complex

#ifndef M_PI
#define M_PI 3.1415926535897931
#endif

int wait=20;


/* will be explained later */
double norm_estimate=3.55;

/* enum for kind of opitmization criterion */

enum optikind {OPTINO,OPTIDEVIATION,OPTIMAX,OPTIRELERROR,OPTIGA};

inline double flat(){  return (double)rand()/(double)RAND_MAX;}

/**
 * this computes the routs of the chebycheff polynomial approximating 1/x
 */


static void RootsOfOneOverX(COMPLEX* roots,const int degree,
			    double epsilon){

  int i;
  double arg;

  for(i=1;i<=degree;i++){
    arg=2*M_PI*(double)i/((double)degree+1.0) ;
    roots[i-1]=0.5*(1.0+epsilon)*(1.0-cos(arg))
	       -I*sqrt(epsilon)*sin(arg);
  }


}

/**
 * this is MY bitreversal routine im proud of
 */

void ReverseBits(unsigned int *bits,int numbits){
  if(sizeof(int)==4){
  *bits= (((*bits) & 0x0000ffff) <<16) + (((*bits) & 0xffff0000)>>16);
  *bits= (((*bits) & 0x00ff00ff) <<8) + (((*bits) & 0xff00ff00)>>8);
  *bits= (((*bits) & 0x0f0f0f0f) <<4) + (((*bits) & 0xf0f0f0f0)>>4);
  *bits= (((*bits) & 0x33333333) <<2) + (((*bits) & 0xCCCCCCCC)>>2);
  *bits= (((*bits) & 0x55555555) <<1) + (((*bits) & 0xAAAAAAAA)>>1);
  *bits >>=(32-numbits);
  } else {fprintf(stderr,"error: int has not 32 bits\n");}
}

/****
 * the most efficient way (from my point of view) to compute integer powers of 2 two
 */
int pow2(int n){
  return 1<<n;
}



/**
 * reorder the roots in the bit reversal fashion
 * this is only working for powers of two ?!
 */
void BitReversalOrder(COMPLEX * roots,const int degree) {

  unsigned  int i,j,z;
  int bits=1;
  /* create a buffer for the reordered roots */
  COMPLEX *roots_buffer=NULL;

  /* allocate space for roots buffer */
  roots_buffer=(COMPLEX*)malloc(degree*sizeof(COMPLEX));
  if(roots_buffer==NULL)
    {fprintf(stderr," error allocating memory\n");
    exit(-1);
  }

  /**
   * here we calculate how many bits we need for the given degree
   * this is a kind of integer valued logarithm to the power of two
   */
  while(pow2(bits)<degree) ++bits;

  for(i=0,z=0;i<degree;){
    j=z;
    /* calculate the bitreversed index */
    ReverseBits(&j,bits);
    /**
     * bitreversion(br) is self inverse this means if the br of a
     * is b then the br of b is a if we would do the exchange of the 
     * roots for all indices and their br counterparts we would do 
     * the exchange twice and would have no br ordering
     * so the next condition ensures that we exchange only once for
     * each br pair
     */
    if(j<degree){
      fprintf(stderr,"%d\n",j);
      roots_buffer[i]=roots[j];
      ++i;
    }
    ++z;
  }
  /* write the br reordered roots back to the roots array */
  for(i=0;i<degree;i++){
    roots[i]=roots_buffer[i];
    fprintf(stderr,"%d %e %e\n",i,creal(roots[i]),cimag(roots[i]));
  }

  /* free the memory */
  free(roots_buffer);
}


/**
 * this is in principle the same as above but it supports also 
 * non power of two degrees of the polynomial
 */

void MyBitReversalOrder(COMPLEX * roots,const int degree) {

  unsigned  int i,j,z;
  int bits=1;
  COMPLEX *roots_buffer=NULL;
  roots_buffer=(COMPLEX*)malloc(degree*sizeof(COMPLEX));
  if(roots_buffer==NULL)
    {fprintf(stderr," error allocating memory\n");
    exit(-1);
  }

  while(pow2(bits)<degree) ++bits;

  /**
   * here we do the bitreversed ordering only for the FIRST
   * half of the roots, the second half is constructed by taking 
   * complex conjugate roots of the first half
   */
  for(i=0,z=0;i<degree/2;){
    j=z;
    ReverseBits(&j,bits);
    if(j<degree){
      fprintf(stderr,"%d\n",j);
      roots_buffer[i]=roots[j];
      roots_buffer[degree-i-1]=roots[degree-j-1];
      ++i;
    }
    ++z;
  }

  for(i=0;i<degree;i++){
    roots[i]=roots_buffer[i];
    fprintf(stderr,"%d %e %e\n",i,creal(roots[i]),cimag(roots[i]));
  }

  free(roots_buffer);
}

/**
 * this should be self explaining
 */

void WriteRootsToFile(char *filename,COMPLEX const *roots,int degree){

  FILE *file;
  int i;
  if((file=fopen(filename,"w"))!=NULL){

    fprintf(file,"Nr    Re   Im\n");

    for(i=0;i<degree;i++){
      fprintf(file,"%d %+1.16e %+1.16e\n",i,creal(roots[i]),cimag(roots[i]));
    }
    fclose(file);
  }

}

/**
 * evaluate the polynomial by roots
 */

static COMPLEX EvalPoly(COMPLEX const *roots,const int degree,COMPLEX value,double normierunglocal,FILE* ofile){
  COMPLEX prod=1.0;
  int i;
  for(i=0;i<degree;i++){
    prod*=value - roots[i];
    prod*=normierunglocal;
    if(ofile!=NULL)
      fprintf(ofile,"%d %e %e\n",i,creal(prod),cimag(prod) );
  }
  return prod;

}

/* */
void reorderRoots(const complex *rootsin,int degree,int *perm,complex* rootsout){
  int i;
  for(i=0;i<degree;i++){
    rootsout[i]=rootsin[perm[i]];
  }
}

/* */
inline void  permPerm(int *perm,int a,int b){
  static int c;
  c=perm[a];
  perm[a]=perm[b];
  perm[b]=c;
}


/* used for optimizing the ordering of the polynomial roots */
double relErrorTestpoints(const complex *roots,int degree,double norm,complex* testPoints,int n_points,complex *prod){
  complex buffer;
 double  overallerror=0;
 int i;
    /* calculate error for test points */
    for(i=0;i<n_points;i++){
      buffer=EvalPoly(roots,degree,testPoints[i],norm,NULL);
      buffer*=testPoints[i];
      buffer-=1;
      buffer=cabs(buffer);
      overallerror+=creal(buffer);
    }
    
    return (overallerror/=(double)n_points);
}


/* used for optimizing the ordering of the polynomial roots */
double optiDeviation(const complex *roots,int degree,double norm,complex* testPoints,int n_points,complex* prod){
  int i,j;
  double globMax=0;
  double max=0,min=0;
  double logProd=0;
  for(i=0;i<n_points;i++){
    prod[i]=testPoints[i];
  }

  for(j=0;j<degree;j++){

    for(i=0;i<n_points;i++){
      prod[i]*=norm*(testPoints[i]-roots[j]);
      logProd=log(cabs(prod[i]));
      if(i==0) min=max=logProd;
      else if(logProd >max) max=logProd;
      else if(logProd <min) min=logProd;
    }
/*     if(j==0) globMax=max-min; */
/*     else if(max-min>globMax) globMax=max-min; */
    globMax+=max-min;
  }

  return globMax/(double)degree;
}


/* used for optimizing the ordering of the polynomial roots */
double optiMaxVal(const complex *roots,int degree,double norm,complex* testPoints,int n_points,complex* prod){
  int i,j;
  double pointsMax=0;
  double globMax=0;
  double cabsprod;
  for(i=0;i<n_points;i++){
    prod[i]=testPoints[i];
  }

  for(j=0;j<degree;j++){

    for(i=0;i<n_points;i++){
      prod[i]*=norm*(testPoints[i]-roots[j]);
      cabsprod=cabs(prod[i]);
      if(i==0) pointsMax=cabsprod;
      else if(cabsprod >pointsMax) pointsMax=cabsprod;
    }
    if(j==0) globMax=pointsMax;
    else if(pointsMax>globMax) globMax=pointsMax;
  }

  return globMax;
}


/*********************************************************** 
 * a small brute-force stupid genetic algorithm to optimize 
 * root ordering
 ***********************************************************/

void OptimizeOrderMC(COMPLEX *roots,int degree,double norm,double epsilon,int n_points,double (*fn)(const complex*,int,double,complex*,int,complex*) ){


  complex *testPoints=malloc(sizeof(complex)*n_points);
  complex *error=malloc(sizeof(complex)*n_points);
  double overallerror,newoverallerror,initialoverallerror;
  complex *newRoots=malloc(sizeof(complex)*degree);
  complex* prod=(complex*)malloc(sizeof(complex)*n_points);

  int *perm=malloc(sizeof(int)*degree);
  int *permsave=malloc(sizeof(int)*degree);
  int *permtmp=malloc(sizeof(int)*degree);
  double dx=(1.-epsilon)/((double)n_points-1.0);
  int i;
  double min=1.0;

  int index1=1;
  int index2=2;
  time_t lastMinTime,actualTime;

  FILE *fileIndex;


  /* create test points */
  for(i=0;i<n_points;i++){
    testPoints[i]=epsilon+dx*(double)(i);
/*     testPoints[i]=epsilon+(1.0-epsilon)*flat(); */
  }

  /* init permutation */
  for(i=0;i<degree;i++){
    perm[i]=i;
    permsave[i]=i;
  }

  /*init newRoots */
  reorderRoots(roots,degree,perm,newRoots);





  /*   initialoverallerror=relErrorTestpoints(newRoots,degree,norm,testPoints,n_points); */
  initialoverallerror=fn(newRoots,degree,norm,testPoints,n_points,prod);

  overallerror=initialoverallerror;
	
  printf("Overall error is %e\n ",overallerror);

  /* initialize time */
  actualTime=lastMinTime=time(NULL);


  i=0;
  min=1.0;

  while(1){

    if(i%1000){
      actualTime=time(NULL);
    }
    /* if 30 seconds passed since the last minimum found -> break */
    if(actualTime-lastMinTime> wait){
      printf("We stop here since we found no new minimum within the last %d seconds\n",wait);
      break;
      
    }
    

    /* save state before shuffeling */
    for(i=0;i<degree;i++)
      permtmp[i]=perm[i];

    /* shuffle indices .... */
    for(i=0;i < 1 + (int)flat() * 1.5  ;  i++){ 
      index1=(int)(flat()*(double)(degree/2));
      index2=(int)(flat()*(double)(degree/2));
      if(flat()<0.3) {
	permPerm(perm,index1,degree-1-index1);
      } else  if(index1!=index2){
	permPerm(perm,index1,index2);
	permPerm(perm,degree-1-index1,degree-1-index2);
      }

    } 

    /* ... and apply to the roots */
    reorderRoots(roots,degree,perm,newRoots);


    /* calculate new overall error */
    /*     newoverallerror=relErrorTestpoints(newRoots,degree,norm,testPoints,n_points*j); */
    newoverallerror=fn(newRoots,degree,norm,testPoints,n_points,prod);

    /* accept/reject step */
    if(newoverallerror<overallerror){
      overallerror=newoverallerror;
      if(overallerror/initialoverallerror<min) {
	printf("Rel Overall error (%e) is %2.6f \% : rel improvement is %e \n ", overallerror,100*min,(min-overallerror/initialoverallerror)/min);

	for(i=0;i<degree;i++)
	  permsave[i]=perm[i];

	lastMinTime=time(NULL);
	min=overallerror/initialoverallerror;
      }
    } else {
      /* here we can add rules to accept even if the error becomes larger */
      if( flat()<0.00) /* accept although */{
	overallerror=newoverallerror;
      } else /* restore old state */{
	for(i=0;i<degree;i++)
	  perm[i]=permtmp[i];
	reorderRoots(roots,degree,perm,newRoots);
      }
    }

    ++i;
  }    


  /* restore last good ordering */
  for(i=0;i<degree;i++)
    perm[i]=permsave[i];


  /* take the saved state ..  */
  reorderRoots(roots,degree,permsave,newRoots);

  if((fileIndex=fopen("indices.csv","w"))!=NULL){
    printf("printing optimized order of roots:\n");
    for(i=0;i< degree;i++){
      fprintf(fileIndex," %d \n", permsave[i]);
    }
    fclose(fileIndex);
  }
  /* and write new roots */
  for(i=0;i<degree;i++)
    roots[i]=newRoots[i];


  /* free up arrays */
  free(prod);
  free(testPoints);
  free(error);
  free(perm);
  free(permsave);
  free(newRoots);



}


/**
 * here we calculate the normierung,
 * for this we need to calculate the polynomial from the roots
 * and have to start with a arbitrary normierung which is norm_estimate
 * this should be close to the expected value of the local normierung
 * because we would get numerical instabilities otherwise
 * e.g. it would be a very bad idea to start with a value of 1.0
 * for norm_estimate
 */
static COMPLEX Normierung(COMPLEX const *roots,int degree,double epsilon){
  COMPLEX x=0.5*(1+epsilon),norm;
  norm=EvalPoly(roots,degree,x,norm_estimate,NULL);
  norm*=x;
  norm=1/norm;
  return norm;
}

int main(int argc,char **argv){

  int degree=16;
  double invdegreepo;
  double epsilon=0.01/2.6;
  COMPLEX *roots;
  double *roots_double_wrapped;
  unsigned int i,j;
  COMPLEX p,err;
  double norm,normierunglocal,normierunglocal_check;
  double normierunglocal_delta,normierunglocal_olddelta;
  double error;
  double divisions=50.;
  double randomPoint;
  COMPLEX dx=(1.0)/divisions,x=0;

  FILE *file;
  FILE *prodFile;
  char filename[256];

  char optchar='\0';

  char *prefix="1overX_poly";


  enum optikind do_optimize=OPTINO;
  int output_prod_hist=0;
  /* init random number generator */

  int num_points_for_opti=-1;

  srand(time(NULL));


  i=0xabcd;
  j=i;
  ReverseBits(&j,16);
  fprintf(stderr," for checking: %x bit reversed is %x (should be %x )\n" , i ,j,0xb3d5);


  /**
   * parse arguments
   */
  while( (optchar=getopt(argc,argv,"o:d:e:pn:w:"))!=-1){
    switch (optchar) {
    case 'o':
      if(strcmp(optarg,"max")==0)
	do_optimize=OPTIMAX;
      else if (strcmp(optarg,"relerror")==0)
	do_optimize=OPTIRELERROR;
      else if (strcmp(optarg,"deviation")==0)
	do_optimize=OPTIDEVIATION;
      else if (strcmp(optarg,"ga")==0) {
#ifdef WITHGALIB
	do_optimize=OPTIGA;
#else
	do_optimize=OPTINO;
	fprintf(stderr,"oox was not compiled with galib support -> Doing NO optimization!\n");
#endif
      } else {
	fprintf(stderr," Argument -o requires an argument: max | relerror | deviation  \n");
	abort();
      }
      break;
    case 'd':
      degree=atol(optarg);
      break;
    case 'e':
      epsilon=atof(optarg);
      break;
    case 'n':
      num_points_for_opti=atoi(optarg);
      break;
    case 'w':
      wait=atoi(optarg);
      break;
    case 'p': 
      output_prod_hist = 1 ;
      break;
    case '?':
      if(optopt == 'd' || optopt == 'e' ) {
	fprintf(stderr," Argument %c requires an argument \n",optopt);
	abort();
      } else if(optopt == 'o' ){
	fprintf(stderr," Argument %c requires an argument: max | relerror | deviation  \n",optopt);
	abort();
      } else if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character `\\x%x'.\n",
		 optopt);
    default: 
      abort();
    }
  }

  invdegreepo=1.0/((double)degree);

  if(num_points_for_opti==-1)
    num_points_for_opti=(int) degree;



  dx=(1.0-epsilon)/divisions;x=epsilon;

  /*allocate memory for the roots*/
  roots=(COMPLEX*)malloc(sizeof(COMPLEX)*degree);
  if(roots==NULL)
    {fprintf(stderr," error allocating memory\n");
    exit(-1);
    }

  /*calculate the roots of the chebycheff polynomial approximating 1/x*/
  RootsOfOneOverX(roots,degree,epsilon);

  fprintf(stderr,"here come the roots\n");

  /* reorder roots in "my manner" of bitreversal order */
  MyBitReversalOrder(roots,degree);


  /*calculate the normierung */
  fprintf(stderr,"norm_estimate local %lf\n",norm_estimate);

  normierunglocal=norm_estimate;
  norm=Normierung(roots,degree,epsilon);
  norm_estimate=pow(norm,invdegreepo)*norm_estimate;

  normierunglocal_olddelta=norm_estimate-normierunglocal;

  fprintf(stderr,"First normierung local %lf (delta) -> %e\n",normierunglocal,normierunglocal_olddelta);

  normierunglocal=norm_estimate;

  j=0;
  do {
    norm=Normierung(roots,degree,epsilon);
    norm_estimate=pow(norm,invdegreepo)*norm_estimate;

    normierunglocal_delta=norm_estimate-normierunglocal;
    if(fabs(normierunglocal_delta)<fabs(normierunglocal_olddelta)){
      normierunglocal=norm_estimate;
    } else break;
    normierunglocal_olddelta=normierunglocal_delta;
    fprintf(stderr,"Second normierung local %lf diff to last %e\n",normierunglocal,normierunglocal-norm_estimate);
  } while( j++<20  );



  /*calculate the roots of the chebycheff polynomial approximating 1/x*/
/*   RootsOfOneOverX(roots,degree,epsilon); */


  /* try to optimize */
  switch (do_optimize) {
  case OPTIDEVIATION:
    OptimizeOrderMC(roots,degree,normierunglocal,epsilon,num_points_for_opti,optiDeviation);
    break;
  case OPTIMAX:
    OptimizeOrderMC(roots,degree,normierunglocal,epsilon,num_points_for_opti,optiMaxVal);
    break;
  case OPTIRELERROR:
    OptimizeOrderMC(roots,degree,normierunglocal,epsilon,num_points_for_opti,relErrorTestpoints);
    break;
  case OPTINO:
    fprintf(stderr,"Doing NO optization of the root ordering !\n");
    break;
  case OPTIGA:
#ifdef WITHGALIB
    roots_double_wrapped=malloc(sizeof(double)*2*degree);
    /*convert to double array */
    for(j=0;j<degree;j++){
      roots_double_wrapped[2*j]=creal(roots[j]);
      roots_double_wrapped[2*j+1]=cimag(roots[j]);
    }

    /* perform ga */
    initGAObject(degree,roots_double_wrapped,normierunglocal,epsilon,num_points_for_opti);

    /*convert from double array */
    for(j=0;j<degree;j++){
      roots[j]=roots_double_wrapped[2*j]+I*roots_double_wrapped[2*j+1];
    }
#endif
    break;
  default:
    fprintf(stderr,"Error unknown type of optimization criterion \n");
    break;
  }

  /* put out product hist */
  if(output_prod_hist)
  for(x=1;creal(x)<=100.0;x+=5){
    randomPoint=1./x;
    sprintf(filename,"%s.%d.%1.16e.prod.at.%lf",prefix,degree,epsilon,randomPoint);
    if((prodFile=fopen(filename,"w"))!=NULL){
      EvalPoly(roots,degree,randomPoint,normierunglocal,prodFile);
      fclose(prodFile);
    }
  }


  /* write the reordered roots to a file */
  sprintf(filename,"%s_deg_%d_eps_%1.16e.roots",prefix,degree,epsilon);
  WriteRootsToFile(filename,roots,degree);

  /* write the normierung to a file */
  sprintf(filename,"%s_deg_%d_eps_%1.16e.const",prefix,degree,epsilon);
  file=fopen(filename,"w");
  if(file!=NULL){
    fprintf(file,"%1.16f\n",normierunglocal);
    fclose(file);
  }

  /* check if what we've written matches the contents of normierung local */
  file=fopen(filename,"r");
  if(file!=NULL){
    fscanf(file,"%lf\n",&normierunglocal_check);
    fclose(file);
    fprintf(stderr, "normierung local - contents of file = %e \n", normierunglocal - normierunglocal_check);
  }


  /* output some values of poly */
  for(i=0,x=0;i<=divisions;i++){
    x=i*dx;
    p=EvalPoly(roots,degree,x,normierunglocal,NULL);
    err=x*p-1.;
    printf("%d %e %e %e %e %e\n",i,creal(x),creal(p),cimag(p),creal(err),cimag(err)); 
  }
    /*  x=0.5;
  p=EvalPoly(roots,degree,x,normierunglocal,NULL);
  printf("%d %e %e %e \n",i,creal(x),creal(p),cimag(p));*/

  /* print the error estimate of the approximation */
  error=2*pow((1-sqrt(epsilon))/(1+sqrt(epsilon)),(double)(degree+1));
  fprintf(stderr,"error of poly is of order %e\n",
	  error );

  /* ... and also the normierung again for checking purpose */
  fprintf(stderr,"normierung local %e\n",
	  normierunglocal );

  /* free memory */
  free(roots);

  /* finnish ;-)*/
  exit(0);

}
