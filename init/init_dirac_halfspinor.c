/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2012 Carsten Urbach
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#if (defined SPI)
# include "DirectPut.h"
#endif
#include "global.h"
#include "su3.h"
#include "init_dirac_halfspinor.h"

#ifdef BGQ
#  define SPI_ALIGN_BASE 0x7f
#else
#  define SPI_ALIGN_BASE ALIGN_BASE
#endif

halfspinor ** NBPointer_;
halfspinor * HalfSpinor_;
halfspinor * HalfSpinor ALIGN;
halfspinor *** NBPointer;
halfspinor * sendBuffer, * recvBuffer;
halfspinor * sendBuffer_, * recvBuffer_;

/* The single precision versions */
halfspinor32 ** NBPointer32_;
halfspinor32 * HalfSpinor32_;
halfspinor32 * HalfSpinor32 ALIGN;
halfspinor32 *** NBPointer32;
halfspinor32 * sendBuffer32, * recvBuffer32;
halfspinor32 * sendBuffer32_, * recvBuffer32_;


int init_dirac_halfspinor() {
  int j=0, k;
  int x, y, z, t;

  NBPointer = (halfspinor***) calloc(4,sizeof(halfspinor**));
  NBPointer_ = (halfspinor**) calloc(16,(VOLUME+RAND)*sizeof(halfspinor*));
  NBPointer[0] = NBPointer_;
  NBPointer[1] = NBPointer_ + (8*(VOLUME+RAND)/2);
  NBPointer[2] = NBPointer_ + (16*(VOLUME+RAND)/2);
  NBPointer[3] = NBPointer_ + (24*(VOLUME+RAND)/2);

  if((void*)(HalfSpinor_ = (halfspinor*)calloc(4*(VOLUME)+1, sizeof(halfspinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }

  HalfSpinor = (halfspinor*)(((unsigned long int)(HalfSpinor_)+ALIGN_BASE+1)&~ALIGN_BASE);

#ifdef MPI
  if((void*)(sendBuffer_ = (halfspinor*)calloc(RAND/2+8, sizeof(halfspinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  sendBuffer = (halfspinor*)(((unsigned long int)(sendBuffer_)+SPI_ALIGN_BASE+1)&~SPI_ALIGN_BASE);
  if((void*)(recvBuffer_ = (halfspinor*)calloc(RAND/2+8, sizeof(halfspinor))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  recvBuffer = (halfspinor*)(((unsigned long int)(recvBuffer_)+SPI_ALIGN_BASE+1)&~SPI_ALIGN_BASE);
#endif

  for(int ieo = 0; ieo < 2; ieo++) {
    for(int i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+1)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(int mu = 0; mu < 4; mu++) {
	NBPointer[ieo][8*i + 2*mu + 0] = &HalfSpinor[ 8*g_lexic2eosub[ g_idn[j][mu] ] + 2*mu + 0];
	NBPointer[ieo][8*i + 2*mu + 1] = &HalfSpinor[ 8*g_lexic2eosub[ g_iup[j][mu] ] + 2*mu + 1];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	k = (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i] = &sendBuffer[ k ];
      }
      if(t == T-1) {
	k = (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
	NBPointer[ieo][8*i + 1] = &sendBuffer[ k ];
      }
#endif
#if ((defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == 0) {
	k = (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 2] = &sendBuffer[ k ];
      }
      if(x == LX-1) {
	k = (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
	NBPointer[ieo][8*i + 3] = &sendBuffer[ k ];
      }
#endif
#if ((defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == 0) {
	k = (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 4] = &sendBuffer[ k ];
      }
      if(y == LY-1) {
	k = (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
	NBPointer[ieo][8*i + 5] = &sendBuffer[ k ];
      }
#endif
#if ((defined PARALLELXYZ) || (defined PARALLELXYZT))
      if(z == 0) {
	k = (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 6] = &sendBuffer[ k ];
      }
      if(z == LZ-1) {
	k = (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
	NBPointer[ieo][8*i + 7] = &sendBuffer[ k ];
      }
#endif
    }
    for(int i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(int mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
#ifdef MPI
#endif
  }
  for(int ieo = 2; ieo < 4; ieo++) {
    for(int i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+0)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(int mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = &HalfSpinor[8*i + mu];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == T-1) {
	NBPointer[ieo][8*i]     = &recvBuffer[ (g_lexic2eosub[ g_iup[j][0] ] - VOLUME/2)];
      }
      if(t == 0) {
	NBPointer[ieo][8*i + 1] = &recvBuffer[ (g_lexic2eosub[ g_idn[j][0] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == LX-1) { 
	NBPointer[ieo][8*i + 2] = &recvBuffer[ (g_lexic2eosub[ g_iup[j][1] ] - VOLUME/2)];
      }
      if(x == 0) {
	NBPointer[ieo][8*i + 3] = &recvBuffer[ (g_lexic2eosub[ g_idn[j][1] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == LY-1) {
	NBPointer[ieo][8*i + 4] = &recvBuffer[ (g_lexic2eosub[ g_iup[j][2] ] - VOLUME/2)];
      }
      if(y == 0) {
	NBPointer[ieo][8*i + 5] = &recvBuffer[ (g_lexic2eosub[ g_idn[j][2] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXYZ) || (defined PARALLELXYZT))
      if(z == LZ-1) {
	NBPointer[ieo][8*i + 6] = &recvBuffer[ (g_lexic2eosub[ g_iup[j][3] ] - VOLUME/2)];
      }
      if(z == 0) {
	NBPointer[ieo][8*i + 7] = &recvBuffer[ (g_lexic2eosub[ g_idn[j][3] ] - VOLUME/2)];
      }
#endif
    }
    for(int i = VOLUME/2; i < (VOLUME+RAND)/2; i++) {
      for(int mu = 0; mu < 8; mu++) {
	NBPointer[ieo][8*i + mu] = NBPointer[ieo][0];
      }
    }
  }
#if (defined SPI && defined MPI)
  // here comes the SPI initialisation
  uint64_t messageSizes[NUM_DIRS];
  uint64_t roffsets[NUM_DIRS], soffsets[NUM_DIRS];

#if (defined PARALLELT)
  spi_num_dirs = 2;
#endif
#if (defined PARALLELXT)
  spi_num_dirs = 4;
#endif
#if (defined PARALLELXYT)
  spi_num_dirs = 6;
#endif
#if (defined PARALLELXYZT)
  spi_num_dirs = 8;
#endif

  totalMessageSize = 0;
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    // message sizes in Bytes
    if(i == 0 || i == 1) messageSizes[i] = LX*LY*LZ*6*sizeof(double);
    else if(i == 2 || i == 3) messageSizes[i] = T*LY*LZ*6*sizeof(double);
    else if(i == 4 || i == 5) messageSizes[i] = T*LX*LZ*6*sizeof(double);
    else if(i == 6 || i == 7) messageSizes[i] = T*LX*LY*6*sizeof(double);

    soffsets[i] = totalMessageSize;
    totalMessageSize += messageSizes[i];
  }
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    // forward here is backward on the right neighbour
    // and the other way around...
    if(i%2 == 0) {
      roffsets[i] = soffsets[i] + messageSizes[i];
    }
    else {
      roffsets[i] = soffsets[i] - messageSizes[i-1];
    }
  }

  Personality_t pers;
  int rc = 0;
  // get the CNK personality
  Kernel_GetPersonality(&pers, sizeof(pers));
  int mypers[6];
  mypers[0] = pers.Network_Config.Acoord;
  mypers[1] = pers.Network_Config.Bcoord;
  mypers[2] = pers.Network_Config.Ccoord;
  mypers[3] = pers.Network_Config.Dcoord;
  mypers[4] = pers.Network_Config.Ecoord;

  get_destinations(mypers);

  // adjust the SPI pointers to the send and receive buffers
  SPIrecvBuffers = (char *)(recvBuffer);
  SPIsendBuffers = (char *)(sendBuffer);

  // Setup the FIFO handles
  rc = msg_InjFifoInit ( &injFifoHandle,
			 0,                      /* startingSubgroupId */
			 0,                      /* startingFifoId     */
			 spi_num_dirs,           /* numFifos   */
			 INJ_MEMORY_FIFO_SIZE+1, /* fifoSize */
			 NULL                    /* Use default attributes */
			 );
  if(rc != 0) {
    fprintf(stderr, "msg_InjFifoInit failed with rc=%d\n",rc);
    exit(1);
  }

  // Set up base address table for reception counter and buffer
  setup_mregions_bats_counters(totalMessageSize);

  // Create descriptors
  // Injection Direct Put Descriptor, one for each neighbour
  SPIDescriptors =
    ( MUHWI_Descriptor_t *)(((uint64_t)SPIDescriptorsMemory+64)&~(64-1));
  create_descriptors(SPIDescriptors, messageSizes, soffsets, roffsets, spi_num_dirs);  

  // test communication
  for(unsigned int i = 0; i < RAND/2; i++) {
    sendBuffer[i].s0.c0 = (double)g_cart_id;
    sendBuffer[i].s0.c1 = (double)g_cart_id;
    sendBuffer[i].s0.c2 = (double)g_cart_id;
    sendBuffer[i].s1.c0 = (double)g_cart_id;
    sendBuffer[i].s1.c1 = (double)g_cart_id;
    sendBuffer[i].s1.c2 = (double)g_cart_id;
  }

  // Initialize the barrier, resetting the hardware.
  rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);
  if(rc) {
    printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
    exit(__LINE__);
  }
  // reset the recv counter 
  recvCounter = totalMessageSize;
  global_barrier(); // make sure everybody is set recv counter
  
  //#pragma omp for nowait
  for (unsigned int j = 0; j < spi_num_dirs; j++) {
    descCount[ j ] =
      msg_InjFifoInject ( injFifoHandle,
			  j,
			  &SPIDescriptors[j]);
  }
  // wait for receive completion
  while ( recvCounter > 0 );

  _bgq_msync();

  j = 0;
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    if(i == 0) k = g_nb_t_up;
    if(i == 1) k = g_nb_t_dn;
    if(i == 2) k = g_nb_x_up;
    if(i == 3) k = g_nb_x_dn;
    if(i == 4) k = g_nb_y_up;
    if(i == 5) k = g_nb_y_dn;
    if(i == 6) k = g_nb_z_up;
    if(i == 7) k = g_nb_z_dn;
    for(int mu = 0; mu < messageSizes[i]/sizeof(halfspinor); mu++) {
      if(k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c0) ||
	 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c1) ||
	 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c2) ||
	 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c0) ||
	 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c1) ||
	 k != (int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s1.c2)) {
	if(g_cart_id == 0) {
	  printf("SPI exchange doesn't work for dir %d: %d != %d at point %d\n", 
		 i, k ,(int)creal(recvBuffer[ soffsets[i]/sizeof(halfspinor) + mu ].s0.c0), mu);
	}
	j++;
      }
    }
  }
  if(j > 0) {
    printf("hmm, SPI exchange failed on proc %d...\n!", g_cart_id);
  }
  else {
    if(g_cart_id == 0) printf("# SPI exchange successfully tested\n");
  }
  
#endif // SPI
  return(0);
}


int init_dirac_halfspinor32() {
  int j=0, k;
  int x, y, z, t, mu;
  
  NBPointer32 = (halfspinor32***) calloc(4,sizeof(halfspinor32**));
  NBPointer32_ = (halfspinor32**) calloc(16,(VOLUME+RAND)*sizeof(halfspinor32*));
  NBPointer32[0] = NBPointer32_;
  NBPointer32[1] = NBPointer32_ + (8*(VOLUME+RAND)/2);
  NBPointer32[2] = NBPointer32_ + (16*(VOLUME+RAND)/2);
  NBPointer32[3] = NBPointer32_ + (24*(VOLUME+RAND)/2);

  if((void*)(HalfSpinor32_ = (halfspinor32*)calloc(8*(VOLUME+RAND)+1, sizeof(halfspinor32))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(-1);
  }

  HalfSpinor32 = (halfspinor32*)(((unsigned long int)(HalfSpinor32_)+ALIGN_BASE)&~ALIGN_BASE);

#ifdef MPI
  //re-use memory from 64Bit version
  sendBuffer32 = (halfspinor32*)sendBuffer;
  recvBuffer32 = (halfspinor32*)recvBuffer;
#endif

  for(int ieo = 0; ieo < 2; ieo++) {
    for(int i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+1)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 4; mu++) {
	NBPointer32[ieo][8*i + 2*mu + 0] = &HalfSpinor32[ 8*g_lexic2eosub[ g_idn[j][mu] ] + 2*mu + 0];
	NBPointer32[ieo][8*i + 2*mu + 1] = &HalfSpinor32[ 8*g_lexic2eosub[ g_iup[j][mu] ] + 2*mu + 1];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == 0) {
	k = (g_lexic2eosub[g_idn[j][0]] - VOLUME/2);
	NBPointer32[ieo][8*i] = &sendBuffer32[ k ];
      }
      if(t == T-1) {
	k = (g_lexic2eosub[g_iup[j][0]] - VOLUME/2);
	NBPointer32[ieo][8*i + 1] = &sendBuffer32[ k ];
      }
#endif
#if ((defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == 0) {
	k = (g_lexic2eosub[g_idn[j][1]] - VOLUME/2);
	NBPointer32[ieo][8*i + 2] = &sendBuffer32[ k ];
      }
      if(x == LX-1) {
	k = (g_lexic2eosub[g_iup[j][1]] - VOLUME/2);
	NBPointer32[ieo][8*i + 3] = &sendBuffer32[ k ];
      }
#endif
#if ((defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == 0) {
	k = (g_lexic2eosub[g_idn[j][2]] - VOLUME/2);
	NBPointer32[ieo][8*i + 4] = &sendBuffer32[ k ];
      }
      if(y == LY-1) {
	k = (g_lexic2eosub[g_iup[j][2]] - VOLUME/2);
	NBPointer32[ieo][8*i + 5] = &sendBuffer32[ k ];
      }
#endif
#if ((defined PARALLELXYZ) || (defined PARALLELXYZT))
      if(z == 0) {
	k = (g_lexic2eosub[g_idn[j][3]] - VOLUME/2);
	NBPointer32[ieo][8*i + 6] = &sendBuffer32[ k ];
      }
      if(z == LZ-1) {
	k = (g_lexic2eosub[g_iup[j][3]] - VOLUME/2);
	NBPointer32[ieo][8*i + 7] = &sendBuffer32[ k ];
      }
#endif
    }
  }
  for(int ieo = 2; ieo < 4; ieo++) {
    for(int i = 0; i < VOLUME/2; i++) {
      j = g_eo2lexic[i + ((ieo+0)%2)*(VOLUME+RAND)/2];
      /* get (t,x,y,z) from j */
      t = j/(LX*LY*LZ);
      x = (j-t*(LX*LY*LZ))/(LY*LZ);
      y = (j-t*(LX*LY*LZ)-x*(LY*LZ))/(LZ);
      z = (j-t*(LX*LY*LZ)-x*(LY*LZ) - y*LZ);
      for(mu = 0; mu < 8; mu++) {
	NBPointer32[ieo][8*i + mu] = &HalfSpinor32[8*i + mu];
      }
#if ((defined PARALLELT) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(t == T-1) {
	NBPointer32[ieo][8*i]     = &recvBuffer32[ (g_lexic2eosub[ g_iup[j][0] ] - VOLUME/2)];
      }
      if(t == 0) {
	NBPointer32[ieo][8*i + 1] = &recvBuffer32[ (g_lexic2eosub[ g_idn[j][0] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELX) || (defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXT) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(x == LX-1) { 
	NBPointer32[ieo][8*i + 2] = &recvBuffer32[ (g_lexic2eosub[ g_iup[j][1] ] - VOLUME/2)];
      }
      if(x == 0) {
	NBPointer32[ieo][8*i + 3] = &recvBuffer32[ (g_lexic2eosub[ g_idn[j][1] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXY) || (defined PARALLELXYZ) || (defined PARALLELXYT) || (defined PARALLELXYZT))
      if(y == LY-1) {
	NBPointer32[ieo][8*i + 4] = &recvBuffer32[ (g_lexic2eosub[ g_iup[j][2] ] - VOLUME/2)];
      }
      if(y == 0) {
	NBPointer32[ieo][8*i + 5] = &recvBuffer32[ (g_lexic2eosub[ g_idn[j][2] ] - VOLUME/2)];
      }
#endif
#if ((defined PARALLELXYZ) || (defined PARALLELXYZT))
      if(z == LZ-1) {
	NBPointer32[ieo][8*i + 6] = &recvBuffer32[ (g_lexic2eosub[ g_iup[j][3] ] - VOLUME/2)];
      }
      if(z == 0) {
	NBPointer32[ieo][8*i + 7] = &recvBuffer32[ (g_lexic2eosub[ g_idn[j][3] ] - VOLUME/2)];
      }
#endif
    }
  }
#if (defined SPI && defined MPI)
  // here comes the SPI initialisation
  uint64_t messageSizes[NUM_DIRS];
  uint64_t roffsets[NUM_DIRS], soffsets[NUM_DIRS];

  uint64_t tMS = 0;
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    // message sizes in Bytes
    if(i == 0 || i == 1) messageSizes[i] = LX*LY*LZ*6*sizeof(float);
    else if(i == 2 || i == 3) messageSizes[i] = T*LY*LZ*6*sizeof(float);
    else if(i == 4 || i == 5) messageSizes[i] = T*LX*LZ*6*sizeof(float);
    else if(i == 6 || i == 7) messageSizes[i] = T*LX*LY*6*sizeof(float);

    soffsets[i] = tMS;
    tMS += messageSizes[i];
  }
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    // forward here is backward on the right neighbour
    // and the other way around...
    if(i%2 == 0) {
      roffsets[i] = soffsets[i] + messageSizes[i];
    }
    else {
      roffsets[i] = soffsets[i] - messageSizes[i-1];
    }
  }

  // Create descriptors
  // Injection Direct Put Descriptor, one for each neighbour
  SPIDescriptors32 =
    ( MUHWI_Descriptor_t *)(((uint64_t)SPIDescriptorsMemory32+64)&~(64-1));
  create_descriptors(SPIDescriptors32, messageSizes, soffsets, roffsets, spi_num_dirs);  

  // test communication
  for(unsigned int i = 0; i < RAND/2; i++) {
    sendBuffer32[i].s0.c0 = (double)g_cart_id;
    sendBuffer32[i].s0.c1 = (double)g_cart_id;
    sendBuffer32[i].s0.c2 = (double)g_cart_id;
    sendBuffer32[i].s1.c0 = (double)g_cart_id;
    sendBuffer32[i].s1.c1 = (double)g_cart_id;
    sendBuffer32[i].s1.c2 = (double)g_cart_id;
  }

  // Initialize the barrier, resetting the hardware.
  int rc = MUSPI_GIBarrierInit ( &GIBarrier, 0 /*comm world class route */);
  if(rc) {
    printf("MUSPI_GIBarrierInit returned rc = %d\n", rc);
    exit(__LINE__);
  }
  // reset the recv counter 
  recvCounter = totalMessageSize/2;
  global_barrier(); // make sure everybody is set recv counter
  
  //#pragma omp for nowait
  for (unsigned int j = 0; j < spi_num_dirs; j++) {
    descCount[ j ] =
      msg_InjFifoInject ( injFifoHandle,
			  j,
			  &SPIDescriptors32[j]);
  }
  // wait for receive completion
  while ( recvCounter > 0 );

  _bgq_msync();

  j = 0;
  for(unsigned int i = 0; i < spi_num_dirs; i++) {
    if(i == 0) k = g_nb_t_up;
    if(i == 1) k = g_nb_t_dn;
    if(i == 2) k = g_nb_x_up;
    if(i == 3) k = g_nb_x_dn;
    if(i == 4) k = g_nb_y_up;
    if(i == 5) k = g_nb_y_dn;
    if(i == 6) k = g_nb_z_up;
    if(i == 7) k = g_nb_z_dn;
    for(int mu = 0; mu < messageSizes[i]/sizeof(halfspinor32); mu++) {
      if(k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s0.c0) ||
	 k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s0.c1) ||
	 k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s0.c2) ||
	 k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s1.c0) ||
	 k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s1.c1) ||
	 k != (int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s1.c2)) {
	if(g_cart_id == 0) {
	  printf("32 Bit SPI exchange doesn't work for dir %d: %d != %d at point %d\n", 
		 i, k ,(int)creal(recvBuffer32[ soffsets[i]/sizeof(halfspinor32) + mu ].s0.c0), mu);
	}
	j++;
      }
    }
  }
  if(j > 0) {
    printf("hmm, SPI exchange failed for 32 Bit halfspinor on proc %d...\n!", g_cart_id);
  }
  else {
    if(g_cart_id == 0) printf("# 32 Bit SPI exchange successfully tested\n");
  }

#endif
  return(0);
}
