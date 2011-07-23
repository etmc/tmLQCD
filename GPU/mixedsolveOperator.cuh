/***********************************************************************
 *
 * Copyright (C)
 * original code from Florian Burger
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  
 * File: mixedsolveOperator.cuh
 *
 * operators for (non-EO) mixed_solver in mixedsolve.cu,
 * derived from template interface class template<class RealT>class mixedsolveOperator;
 * in mixedsolve.cu
 * 
 *
 **************************************************************************/



template<class RealT>class MixedsolveOperatorDirac:public MixedsolveOperator<RealT>
{
protected:
  template<class>friend class MixedsolveOperatorDiracDaggerDirac;
  dim3 tm_dirac_kappaBlockdim;
  dim3 tm_dirac_kappaGriddim;

public:
  int rescalekappa; 

  MixedsolveOperatorDirac(int rescalekappaD=0)
  : tm_dirac_kappaBlockdim(BLOCK,1,1),tm_dirac_kappaGriddim(( VOLUME>=BLOCK ? int(VOLUME/BLOCK)+1 : 1 ),1,1), // this is the partitioning for the Dirac-Kernel
    rescalekappa(rescalekappaD)
  { }


  virtual void gpuInit(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {
    #ifdef USETEXTURE
      //Bind texture gf
      bind_texture_gf(gf);
      //Bind texture spinor to spin4 (D_tm is always applied to spin4)
      bind_texture_spin(spinTmp,1);
    #endif
  }

  virtual void gpu(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {
    // D Ddagger    --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // mu -> -mu for twisted term
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!

    #ifdef USETEXTURE
      unbind_texture_spin(1);
    #endif
      // GAMMA5, mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>> (spinin,spinTmp);
    dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spinTmp,1);
    #endif
      //D_tm 
      dev_tm_dirac_kappa<RealT> <<<tm_dirac_kappaGriddim, tm_dirac_kappaBlockdim>>> (gf, spinTmp, spinout, dev_nn);
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    //GAMMA5 mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>>(spinout,spinTmp);
    dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
      bind_texture_spin(spinTmp,1);
    #endif
    //D_tm
    dev_tm_dirac_kappa<RealT> <<<tm_dirac_kappaGriddim, tm_dirac_kappaBlockdim>>> (gf, spinTmp, spinout, dev_nn);
  }

  virtual void gpuDeinit(dev_spinorM(RealT)* spininout,dev_spinorM(RealT)* spinTmp,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim,const RealT scaleparam)
  {//we have to invert D^+D instead of D because first operator is hermitian , which is a requirement for conjugated gradient algorithm
    if(rescalekappa == 1)
    { //want D^-1 rescaled by 2*kappa
      /// maybe move this block into mixedsolveFunction::gpuDeinit(...) ? - "rescalekappa" can be a public member, which has to set before the dev_cg call
   
      //multiply with D^dagger
    
      #ifdef USETEXTURE
        unbind_texture_spin(1);
      #endif
        dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>>(spininout,spinTmp);
        dev_swapmu <<<1,1>>> ();
      #ifdef USETEXTURE
        bind_texture_spin(spinTmp,1);
      #endif
        dev_tm_dirac_kappa<RealT> <<<tm_dirac_kappaGriddim, tm_dirac_kappaBlockdim >>> (gf, spinTmp, spininout, dev_nn);
        dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>>(spininout,spinTmp);
        dev_swapmu <<<1,1>>> ();
      #ifdef USETEXTURE
        unbind_texture_spin(1);
      #endif


      //go over to non-kappa, Ddagger = g5 D g5
      dev_skalarmult_spinor_field<RealT> <<<linAlgGriddim, linAlgBlockdim >>>(spinTmp,RealT(1.0/(scaleparam*scaleparam)), spininout);  
      //dev_tm_dirac_kappa<<<tm_dirac_kappaGriddim, tm_dirac_kappaBlockdim >>>(gf, spin3, spinout, nn_grid);
    }
  }


  virtual void check(spinor* const conjungateBasisPSpininTmp,spinor* const spinout,const int volume)
  {
    printf("Applying double precision Dirac-Op...\n");
   
    Q_pm_psi_gpu(spinout, conjungateBasisPSpininTmp);
    //diff(residueRSpininout, residueRSpininout, spinTmp ,volume);
  }

  virtual void checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int volume)
  {//=^ multiplication D^+ of inverted (D^+D)^-1 => D^-1
    Q_minus_psi_gpu(spinTmp, spinin);
    assign(spinout, spinTmp, volume);
  }
};


/*template<class RealT>class MixedsolveOperatorDiracDaggerDirac:public MixedsolveOperator<RealT>
{
protected:
  MixedsolveOperatorDirac<RealT> operatorDirac;

public:

  MixedsolveOperatorDiracDaggerDirac()
  : operatorDirac(0)
  { }


  virtual void gpuInit(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {  operatorDirac.gpuInit(spinin,spinTmp,spinout,gf,dev_nn,linAlgGriddim,linAlgBlockdim);  }


  virtual void gpu(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {
    // D Ddagger    --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // mu -> -mu for twisted term
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!

    #ifdef USETEXTURE
     unbind_texture_spin(1);//because it is bind to spin2==spinTmp
    #endif
    #ifdef USETEXTURE
      bind_texture_spin(spinin,1);
    #endif
      //D_tm 
      dev_tm_dirac_kappa<RealT> <<<operatorDirac.tm_dirac_kappaGriddim, operatorDirac.tm_dirac_kappaBlockdim>>> (gf, spinin, spinTmp, dev_nn);
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
      // GAMMA5, mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>> (spinTmp,spinout);
    dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spinout,1);
    #endif
      //D_tm 
      dev_tm_dirac_kappa<RealT> <<<operatorDirac.tm_dirac_kappaGriddim, operatorDirac.tm_dirac_kappaBlockdim>>> (gf, spinout, spinTmp, dev_nn);
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    //GAMMA5 mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>>(spinTmp,spinout);
    dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spinin,1);
    #endif

   //dev_skalarmult_add_assign_spinor_field<<<linAlgGriddim,linAlgBlockdim>>>(spinout,mStarSquare,spinin,spinout);  
  }

  
  virtual void check(spinor* const conjungateBasisPSpininTmp,spinor* const spinout,const int volume)
  {
    printf("Applying double precision Dirac-Op...\n");
    Q_pm_psi(spinout, conjungateBasisPSpininTmp);//D^+D statt DD^+ wie ~OperatorDirac

    //assign_add_mul_r(spinout, conjungateBasisPSpininTmp, mStarSquare, volume);
    //operatorDirac.check(conjungateBasisPSpininTmp,spinout,volume);
  }

  virtual void checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int volume)
  {  
    assign(spinout, spinin, volume);
  }
};*/


template<class RealT>class MixedsolveOperatorDiracDaggerDirac:public MixedsolveOperator<RealT>
{
protected:
  MixedsolveOperatorDirac<RealT> operatorDirac;

public:

  MixedsolveOperatorDiracDaggerDirac()
  : operatorDirac(0)
  { }


  #ifdef USETEXTURE
    virtual void gpuInit(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
    {  bind_texture_gf(gf);  }//Bind texture gf
  #endif

  virtual void gpu(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {
    // D Ddagger    --   Ddagger = gamma5 D gamma5  for Wilson Dirac Operator
    // mu -> -mu for twisted term
    // DO NOT USE tm_dirac_dagger_kappa here, otherwise spin2 will be overwritten!!!

    #ifdef USETEXTURE
      bind_texture_spin(spinin,1);//correct?
    #endif
      //D_tm 
      dev_tm_dirac_kappa<RealT> <<<operatorDirac.tm_dirac_kappaGriddim, operatorDirac.tm_dirac_kappaBlockdim>>> (gf, spinin, spinTmp, dev_nn);
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
      // GAMMA5, mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>> (spinTmp,spinout);
    dev_swapmu <<<1,1>>> ();
    #ifdef USETEXTURE
     bind_texture_spin(spinout,1);
    #endif
      //D_tm 
      dev_tm_dirac_kappa<RealT> <<<operatorDirac.tm_dirac_kappaGriddim, operatorDirac.tm_dirac_kappaBlockdim>>> (gf, spinout, spinTmp, dev_nn);
    #ifdef USETEXTURE
     unbind_texture_spin(1);
    #endif
    //GAMMA5 mu -> -mu
    dev_gamma5<RealT> <<<linAlgGriddim, linAlgBlockdim>>>(spinTmp,spinout);
    dev_swapmu <<<1,1>>> ();

   //dev_skalarmult_add_assign_spinor_field<<<linAlgGriddim,linAlgBlockdim>>>(spinout,mStarSquare,spinin,spinout);  
  }

  
  virtual void check(spinor* const conjungateBasisPSpininTmp,spinor* const spinout,const int volume)
  {
    printf("Applying double precision Dirac-Op...\n");
    Q_pm_psi(spinout, conjungateBasisPSpininTmp);//D^+D statt DD^+ wie ~OperatorDirac

    //assign_add_mul_r(spinout, conjungateBasisPSpininTmp, mStarSquare, volume);
    //operatorDirac.check(conjungateBasisPSpininTmp,spinout,volume);
  }

  virtual void checkDeinit(spinor* const spinin,spinor* const spinTmp,spinor* const spinout,int volume)
  {  assign(spinout, spinin, volume);  }
};


template<class RealT>class MixedsolveOperatorDiracDaggerDiracDiracDaggerDirac:public MixedsolveOperatorDiracDaggerDirac<RealT>
{
public:

  virtual void gpu(dev_spinorM(RealT)* spinin,dev_spinorM(RealT)* spinTmp,dev_spinorM(RealT)* spinout,dev_su3_2vM(RealT)* gf,int* dev_nn,const dim3& linAlgGriddim,const dim3& linAlgBlockdim)
  {
    MixedsolveOperatorDiracDaggerDirac<RealT>::gpu(spinin ,spinTmp,spinout,gf,dev_nn,linAlgGriddim,linAlgBlockdim);
    MixedsolveOperatorDiracDaggerDirac<RealT>::gpu(spinout,spinTmp,spinout,gf,dev_nn,linAlgGriddim,linAlgBlockdim);
  }

  
  virtual void check(spinor* const conjungateBasisPSpininTmp,spinor* const spinout,const int volume)
  {
    printf("Applying double precision Dirac-Op 2x...\n");
    Q_pm_psi(spinout, conjungateBasisPSpininTmp);//D^+D statt DD^+ wie ~OperatorDirac
    assign(conjungateBasisPSpininTmp, spinout, volume);
    Q_pm_psi(spinout, conjungateBasisPSpininTmp);//D^+D statt DD^+ wie ~OperatorDirac
  }
};

