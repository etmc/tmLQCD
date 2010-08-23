
#include <ga/ga.h>
#include <ga/garandom.h>

#include <complex>
//#include "MyRootOrderGenome.h"


 #include "oox_gawrapper.h"

using namespace std;



struct poly_params {
  int degree;
  std::complex<double> *roots;
  double norm;
  double epsilon;
  int n_points;
  double *points;
  std::complex<double> *prod;
  bool *indexMap1;
  bool *indexMap2;
};


float objectiveFn(GAGenome &g);
void IdInitializer(GAGenome &g);
void RandomInitializer(GAGenome &g);

int SwapMutator(GAGenome &g,float pmut);
int SimpleSwapMutator(GAGenome &g,float pmut);
int MyUniformCrossover(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
bool check(GA1DArrayGenome<int> &g,const char *c);


double objectiveDeviation(poly_params* pp,GA1DArrayGenome<int> &g);


complex<double> testPP(poly_params* pp,complex<double> value=complex<double>(0.5,0.0)){
  complex<double> prod(1.0,0.);
  int i;
  for(i=0;i<pp->degree;i++){
    prod*=value - pp->roots[i];
    prod*=pp->norm;
  }
  return prod;
}


void initGAObject(int degree,double *roots,double norm,double epsilon,int n_points){

   complex<double> *wr_roots=new complex<double>[degree];
   unsigned int seed=0;

   GARandomSeed(seed);

   for(int i = 0;i<degree;i++){
     wr_roots[i]=complex<double>(roots[2*i],roots[2*i+1]);
   }


   poly_params *pp=new poly_params;


   pp->degree=degree;
   pp->roots=wr_roots;
   pp->norm=norm;
   pp->epsilon=epsilon;
   pp->n_points=n_points;


   pp->points=new double[pp->n_points];

   pp->prod=new complex<double>[pp->n_points];
   pp->indexMap1=new bool[pp->degree];
   pp->indexMap2=new bool[pp->degree];


  double dx=(1.-epsilon)/((double)pp->n_points-1.0);

  /* create test points */
  for(int i=0;i<pp->n_points;i++){
    pp->points[i]=epsilon+dx*(double)(i);
//     pp->points[i]=epsilon+(1.0-epsilon)*GARandomDouble();
  }


  cout << "result of testPP = " << testPP(pp) << "\n";


  GA1DArrayGenome<int> newobj(degree,objectiveFn,(void*)pp);

  if(pp->degree < 200 ){
    newobj.initializer(RandomInitializer);
  } else {
    newobj.initializer(IdInitializer);
  }


  newobj.mutator(SwapMutator);
  newobj.crossover(MyUniformCrossover);


  GASimpleGA ga(newobj);


  ga.parameters("settings.txt");

  ga.pMutation(ga.pMutation()*10.0/(double)pp->degree);

  ga.evolve(seed);

  //  ga.evolve();

  cout << ga.statistics() << "\n";

  const GA1DArrayGenome<int> &bestgenome=DYN_CAST(const GA1DArrayGenome<int> &,
				       ga.statistics().bestIndividual());

  cout << " here comes the best individual \n";
  cout << bestgenome << "\n";

   for(int i = 0;i<degree;i++){
     roots[2*i]=real(pp->roots[bestgenome.gene(i)]);
     roots[2*i+1]=imag(pp->roots[bestgenome.gene(i)]);
   }

  delete [] pp->roots;
  delete [] pp->points;
  delete [] pp->prod;
  delete [] pp->indexMap1;
  delete [] pp->indexMap2;
  delete pp;

}

float objectiveFn(GAGenome &g){
  GA1DArrayGenome<int> &child=DYN_CAST(GA1DArrayGenome<int>&,g);
  double objVal=(double)objectiveDeviation((poly_params*)child.userData(),child);
  cout << "Objval = " << objVal;
  return 1.0/ objVal ;
}


void IdInitializer(GAGenome &g){

  GA1DArrayGenome<int> &child=DYN_CAST(GA1DArrayGenome<int>&,g);

  int n=child.length();
  for( int i=0;i<n;i++){
    child.gene(i,i);
  }
  check(child,"IdInitializer");
}

void RandomInitializer(GAGenome &g){

  GA1DArrayGenome<int> &child=DYN_CAST(GA1DArrayGenome<int>&,g);

  bool *indexMap=((poly_params*)child.userData())->indexMap1;

  int n=child.length();
  int rnd;
  int j,c;

  for(int i=0;i<n;i++){
    indexMap[i]=true;
  }

  for( int i=0;i<n/2;i++){
    rnd = GARandomInt(0,n-1);


    c=0;j=0;
    while(1){
      j=j%n;
      if(indexMap[j] && indexMap[n-1-j]) {
	++c;
	if(c>=rnd) break;
      }
      ++j;
    }

    child.gene(i,j);
    child.gene(n-1-i,n-1-j);

    indexMap[j]=false;
    indexMap[n-1-j]=false;

  }
  check(child,"RandomInitializer");
}




/* used for optimizing the ordering of the polynomial roots */
double objectiveDeviation(poly_params* pp,GA1DArrayGenome<int> &g){
  int i,j;
  double globMax=0;
  double max=0,min=0;
  double logProd=0,absprod;
  const  double eps=1.e-10;
  const  double large=1.e16;
  bool first;

  for(i=0;i<pp->n_points;i++){
    pp->prod[i]=complex<double>(pp->points[i],0.0);
  }

  for(j=0;j<pp->degree;j++){

    first=true;
    for(i=0;i<pp->n_points;i++){
      pp->prod[i]*=pp->norm*(pp->points[i]-pp->roots[g.gene(j)]);
      absprod=abs(pp->prod[i]);
      
//       if(absprod>large){ absprod=large;}
//       if(absprod<eps){ absprod=eps;}

      logProd=log(absprod);
      if(first) {min=max=logProd; first=false;}
      else if(logProd >max) max=logProd;
      else if(logProd <min) min=logProd;
    }
/*     if(j==0) globMax=max-min; */
/*     else if(max-min>globMax) globMax=max-min; */
    globMax+=max-min;
  }

  return globMax/pp->degree;
}

bool check(GA1DArrayGenome<int> &g,const char *msg=" "){
  return true;
  bool *indexMap=((poly_params*)g.userData())->indexMap1;
  int n=g.length();

  for(int i=0;i<n;i++)
    indexMap[i]=true;

  for(int i=0;i<n/2;i++){
    if(g.gene(i)!=n-1-g.gene(n-1-i) || !(indexMap[g.gene(i)] && indexMap[n-1-g.gene(i)]) ){
      fprintf(stdout,"Error validity check failed : %s \n",msg);
      cout << "Object:\n" << g << "\n";
      return false;
    }

    indexMap[g.gene(i)]=false;
  }

  return true;

}

int SwapMutator(GAGenome &g,float pmut){

  GA1DArrayGenome<int> &child=DYN_CAST(GA1DArrayGenome<int>&,g);

  register int n, i;
  int index1, index2;
  if(pmut <= 0.0) return(0);

  float nMut = pmut * STA_CAST(float,child.length());
  int length = child.length();
  if(nMut < 1.0){		// we have to do a flip test on each bit
    nMut = 0;
    for(i=0; i<length/2; i++){
      if(GAFlipCoin(pmut)){
	index1=i;
	index2=GARandomInt(0, length/2-1);
	if(index1!=index2){
	  child.swap(index1, index2);
	  child.swap(length-1-index1, length-1-index2);
	  nMut+=2;
	} else {
	  child.swap(index1, length-1-index1);
	  nMut+=1;
	}
      }
    }
  }
  else{				// only flip the number of bits we need to flip
    for(n=0; n<nMut; n++){
      index1=GARandomInt(0, length/2-1);
      index2=GARandomInt(0, length/2-1);
      if(index1!=index2){
	child.swap(index1, index2);
	child.swap(length-1-index1, length-1-index2);
	n+=2;
      } else {
	child.swap(index1, length-1-index1);
	n+=1;
      }
    }
  }
  check(child,"SwapMutator");
  return(STA_CAST(int,nMut));

}

int SimpleSwapMutator(GAGenome &g,float pmut){

  GA1DArrayGenome<int> &child=DYN_CAST(GA1DArrayGenome<int>&,g);

  int index1,index2;

  int length = child.length();
  index1=GARandomInt(0, length/2-1);
  index2=GARandomInt(0, length/2-1);

  if(index1!=index2){
    child.swap(index1, index2);
    child.swap(length-1-index1, length-1-index2);

    check(child,"SimpleSwapMutator");
    return 2;
  } else {
    child.swap(index1, length-1-index1);
    check(child,"SimpleSwapMutator");
    return 1;
  }

}


// #define RECURSIVE

 int MyUniformCrossover(const GAGenome& p1, const GAGenome& p2,
		 GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<int> &mom=DYN_CAST(const GA1DArrayGenome<int> &, p1);
  const GA1DArrayGenome<int> &dad=DYN_CAST(const GA1DArrayGenome<int> &, p2);

  int n=0;
  int i;
  int mgi,dgi;

  if(c1 && c2){
    GA1DArrayGenome<int> &sis=DYN_CAST(GA1DArrayGenome<int> &, *c1);
    GA1DArrayGenome<int> &bro=DYN_CAST(GA1DArrayGenome<int> &, *c2);

    bool *indexMapSis=((poly_params*)sis.userData())->indexMap1;
    bool *indexMapBro=((poly_params*)sis.userData())->indexMap2;

    if(sis.length() == bro.length() &&
       mom.length() == dad.length() &&
       sis.length() == mom.length()){
      int length=sis.length();

      /* initialize sisters and brothers indexmap */
      for(i=0;i<length;i++)
	indexMapSis[i]=true;
      for(i=0;i<length;i++)
	indexMapBro[i]=true;

      for(i=0; i<length/2; i++){

	mgi=mom.gene(i);
	dgi=dad.gene(i);

	if(GARandomBit()){

	  if(indexMapSis[(int)mgi] && indexMapSis[length-1-(int)mgi] ){
	    sis.gene(i, mgi);
	    indexMapSis[(int)mgi]=false;
	  } else if(indexMapSis[(int)dgi] && indexMapSis[length-1-(int)dgi] ){
	    sis.gene(i, dgi);
	    indexMapSis[(int)dgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=mom;
	    bro=dad;
	    return 0;
#endif
	  }
	  
	  if(indexMapBro[(int)dgi] && indexMapBro[length-1-(int)dgi] ){
	    bro.gene(i, dgi);
	    indexMapBro[(int)dgi]=false;
	  } else if(indexMapBro[(int)mgi] && indexMapBro[length-1-(int)mgi] ){
	    bro.gene(i, mgi);
	    indexMapBro[(int)mgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=mom;
	    bro=dad;
	    return 0;
#endif
	  }
	}
	else{


	  if(indexMapSis[(int)dgi] && indexMapSis[length-1-(int)dgi] ){
	    sis.gene(i, dgi);
	    indexMapSis[(int)dgi]=false;
	  } else  if(indexMapSis[(int)mgi] && indexMapSis[length-1-(int)mgi] ) {
	    sis.gene(i, mgi);
	    indexMapSis[(int)mgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=dad;
	    bro=mom;
	    return 0;
#endif
	  }
	  
	  if(indexMapBro[(int)mgi] && indexMapBro[length-1-(int)mgi] ){
	    bro.gene(i, mgi);
	    indexMapBro[(int)mgi]=false;
	  } else if(indexMapBro[(int)dgi] && indexMapBro[length-1-(int)dgi] ){
	    bro.gene(i, dgi);
	    indexMapBro[(int)dgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=dad;
	    bro=mom;
	    return 0;
#endif
	  }

	}
      }

      for(i=length/2; i<length; i++){
	sis.gene(i,length-1-sis.gene(length-1-i));
	bro.gene(i,length-1-bro.gene(length-1-i));
      }

      if( !check(sis,"Sister MyUniformCrossover") || !check(bro,"Brother MyUniformCrossover")){
	cout << "dad\n" << dad << "\nmom\n" << mom;
      }

    }
    n = 2;
  }
  else if(c1 || c2){
    GA1DArrayGenome<int> &sis = (c1 ? 
			       DYN_CAST(GA1DArrayGenome<int> &, *c1) : 
			       DYN_CAST(GA1DArrayGenome<int> &, *c2));

    if(mom.length() == dad.length() && sis.length() == mom.length()){

      int length=sis.length();
      bool *indexMapSis=((poly_params*)sis.userData())->indexMap1;

      /* initialize sisters and brothers indexmap */
      for(i=0;i<length;i++)
	indexMapSis[i]=true;

      for(i=0; i<length/2; i++){

	mgi=mom.gene(i);
	dgi=dad.gene(i);


	if(GARandomBit()){
	  if(indexMapSis[(int)mgi] && indexMapSis[length-1-(int)mgi] ){
	    sis.gene(i, mgi);
	    indexMapSis[(int)mgi]=false;
	  } else if(indexMapSis[(int)dgi] && indexMapSis[length-1-(int)dgi] ){
	    sis.gene(i, dgi);
	    indexMapSis[(int)dgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=mom;
	    return 0;
#endif
	  }
	} else {

	  if(indexMapSis[(int)dgi] && indexMapSis[length-1-(int)dgi] ){
	    sis.gene(i, dgi);
	    indexMapSis[(int)dgi]=false;
	  } else if(indexMapSis[(int)mgi] && indexMapSis[length-1-(int)mgi] ){
	    sis.gene(i, mgi);
	    indexMapSis[(int)mgi]=false;
	  } else {
#ifdef RECURSIVE
	    return MyUniformCrossover(p1,p2,c1,c2);
#else
	    sis=dad;
	    return 0;
#endif
	  }


	}
      }

      for(i=length/2; i<length; i++){
	sis.gene(i,length-1-sis.gene(length-1-i));
      }

      if( !check(sis,"Sister MyUniformCrossover") ){
	cout << "dad\n" << dad << "\nmom\n" << mom;
      }

    }
    n = 1;
  }

  return n;
}
