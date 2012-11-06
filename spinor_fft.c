
#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#ifdef MPI
# include <mpi.h>
#endif
#include "global.h"
#include <io/params.h>
#include <io/gauge.h>
#include <io/spinor.h>
#include <io/utils.h>
#include "start.h"
#include "linalg_eo.h"

#include "spinor_fft.h"
#include "mpi_init.h"
#include "init/init.h"

#ifdef HAVE_FFTW
  #include <fftw3.h>
#endif

#include <string.h>

void  spinor_fft_print_reduct_dims(int *remaining_dims,FILE *logFile);

#ifdef MPI
void check_mpi_comm_membership(MPI_Comm commself,MPI_Comm commcheck,const char *name_a,const char *name_b,FILE *logFile);
#endif

#ifdef HAVE_FFTW
fftw_plan spinor_fftw_plan2d(spinor *spinor_in,spinor *spinor_out,int dim0,int dim1,int howmany,unsigned int forward,int fftw_flags);
#endif

void spinor_fft_transpose_xp_t(spinor *fieldout,spinor* fieldin,int dim0,int dim1,int forward,double mulp);



/**
 * accumulates pieces of the spinor field on nodes with index 0 in the dimensions given in which
 * the collected data is returned
 */
void spinor_fft_reduce_2d(spinor *localSpinorField,int *collectionRank,spinor*** field_collection,spinor **membuff){
  /* this implementation is intended for four dimensional parallelisation */
#if (defined  PARALLELXYZT  && defined MPI && defined HAVE_FFTW)

  int sendRecvCoord[4];
  int i;
  int dims[]={g_nproc_t,g_nproc_x,g_nproc_y,g_nproc_z};


  /* logfile variables */
  char *logFilePrefix="Process";
  char logFileName[512];
  FILE *logFile;
  const int MSG_LOCALDATA = 457;
  MPI_Status ierr;
  MPI_Datatype mpi_local_spinor;
  const int which[]={0,1};


  (*field_collection)=NULL;
  (*membuff)=NULL;

/*   int result; */
  sprintf(logFileName,"./%s_%02d.log",logFilePrefix,g_cart_id);
  logFile=fopen(logFileName,"a");


  MPI_Type_contiguous(VOLUME, field_point, &mpi_local_spinor);
  MPI_Type_commit(&mpi_local_spinor);


  for(i=0;i<4;i++)
    sendRecvCoord[i]=g_proc_coords[i];

  if( g_proc_coords[which[0]] == 0 && g_proc_coords[which[1]] == 0 ){

      /* i am one of the nodes where data is accumulated */
      spinor **accu_field;
      spinor **fft_field;
      spinor *memory_buffer_accu_field;
      spinor *memory_buffer_fft_field;
      int REDUCTIONVOLUME=1;
      int recvRank;
      MPI_Request *requests;
      MPI_Status *status;
      int request_count=0;
      int num_requests;
      fftw_plan local_2d_fft_forward;

      *collectionRank=TRUE;

      /* calculate the number of reduced 2d volume accumulated in this node */
      
      /* number of spinor fields in local units */
      REDUCTIONVOLUME*=dims[which[0]]*dims[which[1]];

      /* number of receive messages */
      num_requests=REDUCTIONVOLUME-1;

      /* reserve space for receive messages */
      requests=(MPI_Request*)malloc(sizeof(MPI_Request)*num_requests);
      status=(MPI_Status*)malloc(sizeof(MPI_Status)*num_requests);

      fprintf(logFile,"reduction volume = %d\n",REDUCTIONVOLUME);

      /* allocate space for spinor field collection */
      allocate_spinor_field_array(&accu_field,&memory_buffer_accu_field,VOLUME,REDUCTIONVOLUME);
      allocate_spinor_field_array(&fft_field,&memory_buffer_fft_field,VOLUME,REDUCTIONVOLUME);


      /* receive from certain nodes pieces of the spinor field */
      for(sendRecvCoord[which[0]] = 0 ; sendRecvCoord[which[0]]< dims[which[0]] ; sendRecvCoord[which[0]]++){
	for(sendRecvCoord[which[1]] = 0 ; sendRecvCoord[which[1]]< dims[which[1]] ; sendRecvCoord[which[1]]++){
	  if( sendRecvCoord[which[0]] != 0 || sendRecvCoord[which[1]]  != 0){

	    MPI_Cart_rank(g_cart_grid,sendRecvCoord,&recvRank);

	    MPI_Irecv(accu_field[sendRecvCoord[which[0]]*dims[which[1]]+sendRecvCoord[which[1]] ] /* buffer */,
		     1, /* how may */
		     mpi_local_spinor, /* mpi data type */
		     recvRank, /* from whom i get it */
		     MSG_LOCALDATA, /* msg id */
		     g_cart_grid, /* communicator , status */
		     requests+request_count);
	    ++request_count;

	  }
	}
      }


      /* wait until all request finished */
      MPI_Waitall(num_requests, requests, status);

      assign(accu_field[0],localSpinorField,VOLUME);

      /* transpose in xp-t space */
      spinor_fft_transpose_xp_t(fft_field[0],accu_field[0],dims[0],dims[1],TRUE,1.);

      /* create fftw plan */
      local_2d_fft_forward=spinor_fftw_plan2d(fft_field[0],accu_field[0],T*dims[0],LX*dims[1],LY*LZ,1,FFTW_ESTIMATE);
      fftw_execute(local_2d_fft_forward);
      fftw_destroy_plan(local_2d_fft_forward);

/*       assign(accu_field[0],fft_field[0],VOLUME*REDUCTIONVOLUME); */


      free_spinor_field_array(&memory_buffer_fft_field); memory_buffer_fft_field=NULL;

/*       free_spinor_field_array(&memory_buffer_accu_field); memory_buffer_accu_field=NULL; */
      (*field_collection)=accu_field;
      (*membuff)=memory_buffer_accu_field;
      free(requests); requests = NULL;
      free(status); status=NULL;

    } else {
      int sendRank;
      MPI_Request request;
      MPI_Status status;

      *collectionRank=FALSE;

      /* coordinates of the "root" */
      sendRecvCoord[which[0]]=0;
      sendRecvCoord[which[1]]=0;

      MPI_Cart_rank(g_cart_grid,sendRecvCoord,&sendRank); 

      MPI_Isend(localSpinorField,1,mpi_local_spinor,sendRank,MSG_LOCALDATA,g_cart_grid,&request);

      MPI_Wait(&request,&status);

    }


    MPI_Type_free(&mpi_local_spinor);

    fclose(logFile);

#else
    if(g_proc_id==0)
      fprintf(stderr,"Error: Please choose FOUR dimensional parallelization!!!\n");

#endif
}


/**
 * accumulates pieces of the spinor field on nodes with index 0 in the dimensions given in which
 * the collected data is returned
 */
void spinor_fft_redist_2d(spinor *localSpinorField,int collectionRank,spinor** field_collection,spinor *membuff){
  /* this implementation is intended for four dimensional parallelisation */
#if ( defined PARALLELXYZT && defined MPI && defined HAVE_FFTW)

  int sendRecvCoord[4];
  int i;
  int dims[]={g_nproc_t,g_nproc_x,g_nproc_y,g_nproc_z};


  /* logfile variables */
  char *logFilePrefix="Process";
  char logFileName[512];
  FILE *logFile;
  const int MSG_LOCALDATA = 5687;
  MPI_Status ierr;
  MPI_Datatype mpi_local_spinor;
  const int which[]={0,1};



/*   int result; */
  sprintf(logFileName,"./%s_%02d.log",logFilePrefix,g_cart_id);
  logFile=fopen(logFileName,"a");


  /* new mpi type */
  MPI_Type_contiguous(VOLUME, field_point, &mpi_local_spinor);
  MPI_Type_commit(&mpi_local_spinor);


  for(i=0;i<4;i++)
    sendRecvCoord[i]=g_proc_coords[i];

    if( collectionRank == TRUE ){

      /* i am one of the nodes where data is accumulated */
      spinor **accu_field=field_collection;
      spinor **fft_field;
      spinor *memory_buffer_accu_field=membuff;
      spinor *memory_buffer_fft_field;
      int REDUCTIONVOLUME=1;
      int sendRank;
      MPI_Request *requests;
      MPI_Status *status;
      int request_count=0;
      int num_requests;
      fftw_plan local_2d_fft_backward;


      /* calculate the number of reduced 2d volume accumulated in this node */
      
      /* number of spinor fields in local units */
      REDUCTIONVOLUME*=dims[which[0]]*dims[which[1]];

      /* number of receive messages */
      num_requests=REDUCTIONVOLUME-1;

      /* reserve space for receive messages */
      requests=(MPI_Request*)malloc(sizeof(MPI_Request)*num_requests);
      status=(MPI_Status*)malloc(sizeof(MPI_Status)*num_requests);

      fprintf(logFile,"reduction volume = %d\n",REDUCTIONVOLUME);

      /* allocate space for spinor field collection */
      allocate_spinor_field_array(&fft_field,&memory_buffer_fft_field,VOLUME,REDUCTIONVOLUME);



      /* create fftw plan */
      local_2d_fft_backward=spinor_fftw_plan2d(accu_field[0],fft_field[0],T*dims[0],LX*dims[1],LY*LZ,0,FFTW_ESTIMATE);
      fftw_execute(local_2d_fft_backward);
      fftw_destroy_plan(local_2d_fft_backward);


/*       assign(fft_field[0],accu_field[0],VOLUME*REDUCTIONVOLUME); */

      /* transpose in xp-t space */
      spinor_fft_transpose_xp_t(accu_field[0],fft_field[0],dims[0],dims[1],FALSE,1./(double)(T*dims[0] * LX*dims[1]));



      /* receive from certain nodes pieces of the spinor field */
      for(sendRecvCoord[which[0]] = 0 ; sendRecvCoord[which[0]]< dims[which[0]] ; sendRecvCoord[which[0]]++){
	for(sendRecvCoord[which[1]] = 0 ; sendRecvCoord[which[1]]< dims[which[1]] ; sendRecvCoord[which[1]]++){
	  if( sendRecvCoord[which[0]] != 0 || sendRecvCoord[which[1]]  != 0){

	    MPI_Cart_rank(g_cart_grid,sendRecvCoord,&sendRank);

	    MPI_Isend(accu_field[sendRecvCoord[which[0]]*dims[which[1]]+sendRecvCoord[which[1]] ] /* buffer */,
		     1, /* how may */
		     mpi_local_spinor, /* mpi data type */
		     sendRank, /* from whom i get it */
		     MSG_LOCALDATA, /* msg id */
		     g_cart_grid, /* communicator , status */
		     requests+request_count);
	    ++request_count;

	  }
	}
      }

      assign(localSpinorField,accu_field[0],VOLUME);



      /* wait until all request finished */
      MPI_Waitall(num_requests, requests, status);


      free_spinor_field_array(&memory_buffer_fft_field); memory_buffer_fft_field=NULL; fft_field=NULL;
      free_spinor_field_array(&memory_buffer_accu_field); memory_buffer_accu_field=NULL; accu_field=NULL;

      free(requests); requests = NULL;
      free(status); status=NULL;

    } else {
      int recvRank;
      MPI_Request request;
      MPI_Status status;


      /* coordinates of the "root" */
      sendRecvCoord[which[0]]=0;
      sendRecvCoord[which[1]]=0;

      MPI_Cart_rank(g_cart_grid,sendRecvCoord,&recvRank); 

      MPI_Irecv(localSpinorField,1,mpi_local_spinor,recvRank,MSG_LOCALDATA,g_cart_grid,&request);

      MPI_Wait(&request,&status);

    }

    MPI_Type_free(&mpi_local_spinor);

    fclose(logFile);

#else
    if(g_proc_id==0)
      fprintf(stderr,"Error: Please choose FOUR dimensional parallelization!!!\n");

#endif
}


#ifdef HAVE_FFTW
fftw_plan spinor_fftw_plan2d(spinor *spinor_in,spinor *spinor_out,int dim0,int dim1,int howmany_wospin,unsigned int forward,int fftw_flags){

/*    int index_s = gsi(get_index(it, ix, iy, iz, T, L)); */
/*    double *xi_ = xi + index_s; */

  int Dim1[2];
/*    cerr << "Trying to create a plan for T=" << T << " L=" << L ; */
/*    cerr.flush(); */

  int rank=2;

  int stride=12*howmany_wospin;
  int dist=1;
  int howmany=12*howmany_wospin;
  fftw_plan plan;


  Dim1[0]=dim0;
  Dim1[1]=dim1;


  if(fftw_flags==-1){fftw_flags=FFTW_ESTIMATE;}
  if(forward){
    plan=fftw_plan_many_dft(rank, Dim1, howmany, (fftw_complex*)spinor_in, NULL, stride, dist, 
				      (fftw_complex*)spinor_out,NULL,stride,dist,
				      FFTW_FORWARD,fftw_flags);
  } else {
    plan=fftw_plan_many_dft(rank, Dim1, howmany, (fftw_complex*)spinor_in, NULL, stride, dist, 
				      (fftw_complex*)spinor_out,NULL,stride,dist,
				      FFTW_BACKWARD,fftw_flags);
  }
/*    if(plan!=NULL) cerr << "  [OK]"<< endl; */
/*    else cerr << "  [FAIL]"<< endl; */
/*    cerr.flush(); */

 return plan;

}
#endif

void spinor_fft_transpose_xp_t(spinor *fieldout,spinor* fieldin,int dim0,int dim1,int forward,double mulp){
  int LXYZ=LX*LY*LZ;
  int xyz,tp,xp,t;
  spinor *spin1,*spin2;
  if(forward == TRUE){
    for(tp=0;tp<dim0;tp++){
      for(xp=0;xp<dim1;xp++){
	for(t=0;t<T;t++){
	  spin1=fieldin  + LXYZ * ((tp * dim1 + xp ) *   T  + t  );
	  /**                             ^      ^       ^    ^
	   *                               \      \     /     |
	   *                                \      \   /     /
	   *                                 \      \ /     /
	   *                                  \      X     /
	   *                                   \    / \   /
	   *                                    \  /   \ /
	   *                                     \/     X
	   *                                     /\    / \
	   *                                    /  \  /   \
	   *                                   /    \/     \
	   *                                  /     /\      \
	   *                                 /     /  \      \
	   *                                /     /    \      \
	   *                               /     /      \      |
	   *                              V     V        V     V        */
	  spin2=fieldout + LXYZ * ((tp *  T   + t  ) * dim1 + xp );
	  for(xyz=0;xyz<LXYZ;xyz++){
	    memcpy(spin2+xyz,spin1+xyz,sizeof(spinor));
	    _vector_mul((spin2+xyz)->s0,mulp,(spin2+xyz)->s0);
	    _vector_mul((spin2+xyz)->s1,mulp,(spin2+xyz)->s1);
	    _vector_mul((spin2+xyz)->s2,mulp,(spin2+xyz)->s2);
	    _vector_mul((spin2+xyz)->s3,mulp,(spin2+xyz)->s3);
	  }
	  /* optionally multiply with mulp */
	}
      }
    }
  } else {
    for(tp=0;tp<dim0;tp++){
      for(xp=0;xp<dim1;xp++){
	for(t=0;t<T;t++){
	  spin1=fieldin  + LXYZ * ((tp *  T   +  t ) * dim1 + xp );
	  /**                             ^      ^       ^    ^
	   *                               \      \     /     |
	   *                                \      \   /     /
	   *                                 \      \ /     /
	   *                                  \      X     /
	   *                                   \    / \   /
	   *                                    \  /   \ /
	   *                                     \/     X
	   *                                     /\    / \
	   *                                    /  \  /   \
	   *                                   /    \/     \
	   *                                  /     /\      \
	   *                                 /     /  \      \
	   *                                /     /    \      \
	   *                               /     /      \      |
	   *                              V     V        V     V        */
	  spin2=fieldout + LXYZ * ((tp * dim1 + xp ) *   T  +  t  );
	  for(xyz=0;xyz<LXYZ;xyz++){
	    memcpy(spin2+xyz,spin1+xyz,sizeof(spinor));
	    _vector_mul((spin2+xyz)->s0,mulp,(spin2+xyz)->s0);
	    _vector_mul((spin2+xyz)->s1,mulp,(spin2+xyz)->s1);
	    _vector_mul((spin2+xyz)->s2,mulp,(spin2+xyz)->s2);
	    _vector_mul((spin2+xyz)->s3,mulp,(spin2+xyz)->s3);
	  }
	  /* optionally multiply with mulp */
	}
      }
    }
  }
}


#ifdef MPI
void check_mpi_comm_membership(MPI_Comm commself,MPI_Comm commcheck,const char *name_a,const char *name_b,FILE *logFile){
  int result;
  fprintf(logFile,"checking %s against %s : \n" , name_a,name_b);
    MPI_Comm_compare(MPI_COMM_SELF,commcheck,&result);
    switch(result){
    case MPI_CONGRUENT: fprintf(logFile,"CONGRUENT\n"); break;
    case MPI_IDENT: fprintf(logFile,"IDENTICAL\n"); break;
    case MPI_SIMILAR: fprintf(logFile,"SIMILAR\n"); break;
    case MPI_UNEQUAL: fprintf(logFile,"UNEQUAL\n"); break;
    default : fprintf(logFile,"unknown relation ??\n");break;
    }
}
#endif

void  spinor_fft_print_reduct_dims(int *remaining_dims,FILE *logFile){
  int i;

  fprintf(logFile,"Reducing spinor_field to dims : ");

    for(i=0;i<4;i++){
      if( remaining_dims[i]==TRUE){
	switch(i){
	case 0: fprintf(logFile," T"); break;
	case 1: fprintf(logFile," X"); break;
	case 2: fprintf(logFile," Y"); break;
	case 3: fprintf(logFile," Z"); break;
	default: fprintf(logFile," sorry we are in QCD, unknown dimension -> extra dimensions ??"); break;
	}
      }
    }
  
    fprintf(logFile,"\n");
  
}
