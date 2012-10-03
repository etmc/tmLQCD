/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
#  include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <stdint.h>
#ifdef MPI
#  include <mpi.h>
#endif
#ifdef OMP
#  include <omp.h>
#endif
#include "global.h"
#include "DirectPut.h"

// actual number of directions
unsigned int spi_num_dirs = NUM_DIRS;
// total Message Size
// in bytes!
uint64_t totalMessageSize;
// Allocate static memory for descriptors
char SPIDescriptorsMemory[ NUM_DIRS * sizeof(MUHWI_Descriptor_t) + 64 ];
char SPIDescriptorsMemory32[ NUM_DIRS * sizeof(MUHWI_Descriptor_t) + 64 ];
// pointer to descriptor array
MUHWI_Descriptor_t *SPIDescriptors;
MUHWI_Descriptor_t *SPIDescriptors32;

const int batsubgroupID = 0;
int do_dynamic      = 1;
// Enable different zone routing modes
uint8_t  zoneRoutingMask = 0;
unsigned zoneRoutingId   = 0;
// stay on bubble bits
uint8_t stayOnBubbleMask  = 0;
unsigned stayOnBubbleFlag = 0;

// pointers to send and receive buffers
char * SPIrecvBuffers;
char * SPIsendBuffers;

// neighbour destination cache
struct { 
  MUHWI_Destination_t dest;
  uint8_t             hintsABCD;
  uint8_t             hintsE;
} nb2dest[NUM_DIRS];

// receive counter
volatile uint64_t recvCounter;

// counter for injected messages
uint64_t descCount[NUM_DIRS];

// base addess table slot for receive buffer and counter
uint32_t recvBufBatId = 0, recvCntrBatId = 1;

// physical address of send buffers
uint64_t sendBufPAddr;

msg_InjFifoHandle_t injFifoHandle;

void setup_mregions_bats_counters(const int bufferSize) {
  const uint64_t buffersSize =  bufferSize;

  // allocate bat entries for the recive buffer and the receive counter
  
  uint32_t batIds[2] = { recvBufBatId, recvCntrBatId };
  MUSPI_BaseAddressTableSubGroup_t batSubGrp;
  
  int rc =  Kernel_AllocateBaseAddressTable( batsubgroupID/*subgrpId*/,
					     &batSubGrp,
					     2,/*nbatids*/
					     batIds,
					     0 /* "User" use */);
  
  if (rc != 0) {
    fprintf(stderr, "Kernel_AllocateBaseAddressTable failed with rc=%d\n", rc);
    exit(1);
  }
  
  // Receive buffer bat is set to the PA addr of the receive buffer
  Kernel_MemoryRegion_t memRegion;
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   SPIrecvBuffers,
				   buffersSize);
  if ( rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    exit(1);
  }
  
  uint64_t paAddr = 
    (uint64_t)SPIrecvBuffers - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  
  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvBufBatId,
			      paAddr );
  
  if(rc != 0) {
    printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
    exit(1);
  }
  
  // Receive counter bat is set to the MU style atomic PA addr of the receive counter
  if( (uint64_t)(&recvCounter) & 0x7 ) {
    printf("ERROR: recv counter is not 8 byte aligned\n");
    exit(1);
  }
  
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   (void *)&recvCounter,
				   sizeof(recvCounter));
  if(rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    exit(1);
  }
  
  paAddr = 
    (uint64_t)&recvCounter - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  
  uint64_t paAddrAtomic =  MUSPI_GetAtomicAddress(paAddr,MUHWI_ATOMIC_OPCODE_STORE_ADD);
  
  rc = MUSPI_SetBaseAddress ( &batSubGrp,
			      recvCntrBatId,
			      paAddrAtomic );
  
  if(rc != 0) {
    printf("MUSPI_SetBaseAddress failed with rc=%d\n",rc);
    exit(1);
  }
  
  // Get the send buffers physical address
  rc = Kernel_CreateMemoryRegion ( &memRegion,
				   SPIsendBuffers,
				   buffersSize);
  if(rc != 0) {
    printf("Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
    exit(1);
  }
  
  sendBufPAddr = 
    (uint64_t)SPIsendBuffers - 
    (uint64_t)memRegion.BaseVa + 
    (uint64_t)memRegion.BasePa;
  return;
}


void create_descriptors(MUHWI_Descriptor_t * descriptors, uint64_t * messageSizes, uint64_t * soffsets, 
			uint64_t * roffsets, const unsigned int num_dirs) {
  uint64_t anyFifoMap = 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_AP | 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_BP | 
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_CP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_DP |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EM |
    MUHWI_DESCRIPTOR_TORUS_FIFO_MAP_EP;
  
  uint64_t offset;
  static int  did_print =0;
 
  // loop over directions
  // CHECK offset needs to be adjusted for QCD case
  for(unsigned int i = 0; i < num_dirs; i++) {
    // Injection Direct Put Descriptor Information Structure
    MUSPI_Pt2PtDirectPutDescriptorInfo_t dinfo;
    
    memset( (void*)&dinfo, 0x00, sizeof(dinfo) );
      
    dinfo.Base.Payload_Address = sendBufPAddr + soffsets[i];
    dinfo.Base.Message_Length  = messageSizes[i];
    dinfo.Base.Torus_FIFO_Map  = anyFifoMap;
      
    dinfo.Base.Dest = nb2dest[i].dest;
      
    dinfo.Pt2Pt.Hints_ABCD = nb2dest[i].hintsABCD; 

    if(do_dynamic) {	  
      dinfo.Pt2Pt.Misc1 =
	nb2dest[i].hintsE |
	MUHWI_PACKET_USE_DYNAMIC_ROUTING |  
	MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;
      
      dinfo.Pt2Pt.Misc2 = 
	MUHWI_PACKET_VIRTUAL_CHANNEL_DYNAMIC | 
	zoneRoutingMask | 
	stayOnBubbleMask;
      if ( (g_cart_id ==0) && (did_print ==0)) 
	printf("# SPI using dynamic routing  zoneRoutingMask=%d stayOnBubbleMask=%d\n",
	       zoneRoutingMask, stayOnBubbleMask);
    }
    else {	    	    
      dinfo.Pt2Pt.Misc1 =
	nb2dest[i].hintsE |
	MUHWI_PACKET_USE_DETERMINISTIC_ROUTING |  
	MUHWI_PACKET_DO_NOT_ROUTE_TO_IO_NODE;
	
      dinfo.Pt2Pt.Misc2 = 
	MUHWI_PACKET_VIRTUAL_CHANNEL_DETERMINISTIC | 
	zoneRoutingMask | 
	stayOnBubbleMask;
      if ( (g_cart_id ==0) && (did_print ==0)) printf("# SPI using deterministic routing\n");
    }
    did_print++;
    
    dinfo.Pt2Pt.Skip  = 8; // for checksumming, skip the header 	      
    dinfo.DirectPut.Rec_Payload_Base_Address_Id = recvBufBatId;
    dinfo.DirectPut.Rec_Payload_Offset          = roffsets[i];
    dinfo.DirectPut.Rec_Counter_Base_Address_Id = recvCntrBatId;
    dinfo.DirectPut.Rec_Counter_Offset          = 0;
      
    dinfo.DirectPut.Pacing = MUHWI_PACKET_DIRECT_PUT_IS_NOT_PACED;
      
    int rc = MUSPI_CreatePt2PtDirectPutDescriptor(&descriptors[i],
						  &dinfo );
    if (rc != 0) {
      fprintf(stderr, "MUSPI_CreatePt2PtDirectPutDescriptor failed with rc=%d\n",rc);
      exit(1);
    }
  }
}


int get_destinations(int * mypers) {

  int tmp[6];
#if (defined PARALLELT || defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Status mstatus;
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_up, 0, 
	       (void*)tmp, 6, MPI_INT, g_nb_t_dn, 0,
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[1].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_t_dn, 1, 
	       (void*)tmp, 6, MPI_INT, g_nb_t_up, 1, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[0].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
#endif
#if (defined PARALLELXT || defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_up, 2, 
	       (void*)tmp, 6, MPI_INT, g_nb_x_dn, 2, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[3].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_x_dn, 3, 
	       (void*)tmp, 6, MPI_INT, g_nb_x_up, 3, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[2].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
#endif
#if (defined PARALLELXYT || defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_up, 4, 
	       (void*)tmp, 6, MPI_INT, g_nb_y_dn, 4, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[5].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_y_dn, 5, 
	       (void*)tmp, 6, MPI_INT, g_nb_y_up, 5, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[4].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
#endif
#if (defined PARALLELXYZT)
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_up, 6, 
	       (void*)tmp, 6, MPI_INT, g_nb_z_dn, 6, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[7].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
  MPI_Sendrecv((void*)mypers, 6, MPI_INT, g_nb_z_dn, 7, 
	       (void*)tmp, 6, MPI_INT, g_nb_z_up, 7, 
	       g_cart_grid, &mstatus);
  MUSPI_SetUpDestination( &nb2dest[6].dest, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4] );  
#endif
  return(0);
}

typedef struct msg_InjFifoInfo
{
  MUSPI_InjFifoSubGroup_t  subgroup[BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  uint32_t                 numFifosInSubgroup[BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  void                    *fifoMemoryPtr [BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP * 
                                          BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  void                    *fifoPtr       [BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP * 
                                          BGQ_MU_NUM_FIFO_SUBGROUPS_PER_NODE];
  uint32_t                 startingSubgroupId;
  uint32_t                 startingFifoId;
  uint32_t                 numFifos;
  uint32_t                 numSubgroups;
} msg_InjFifoInfo_t;


uint64_t msg_InjFifoInject ( msg_InjFifoHandle_t injFifoHandle,
                             uint32_t            relativeFifoId,
                             MUHWI_Descriptor_t *descPtr ) {
  msg_InjFifoInfo_t *info = (msg_InjFifoInfo_t*)injFifoHandle.pOpaqueObject;
  
  uint32_t globalFifoId = (info->startingSubgroupId * BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP) +
    info->startingFifoId + relativeFifoId;
  
  uint32_t subgroupId   = globalFifoId / BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP;  
  uint64_t rc = MUSPI_InjFifoInject (MUSPI_IdToInjFifo( globalFifoId % BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP,
							&info->subgroup[subgroupId] ),
				     descPtr);
  return rc;
}

void msg_InjFifoTerm ( msg_InjFifoHandle_t injFifoHandle ) {
  return; /*Simple library do nothing! */
}

int msg_InjFifoInit ( msg_InjFifoHandle_t *injFifoHandlePtr,
                      uint32_t             startingSubgroupId,
                      uint32_t             startingFifoId,
                      uint32_t             numFifos,
                      size_t               fifoSize,
                      Kernel_InjFifoAttributes_t  *injFifoAttrs ) {  

  void                *buffer = NULL;
  uint32_t endingFifoId; // Relative to a subgroup
  uint32_t numFifosInSubgroup;
  int rc;
  uint32_t subgroupId = startingSubgroupId;
  uint32_t fifoIds[BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP];
  Kernel_InjFifoAttributes_t attrs[BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP];
  Kernel_InjFifoAttributes_t defaultAttrs;
  uint64_t lock_cache;

  memset ( &defaultAttrs, 0x00, sizeof(defaultAttrs) );
  if(injFifoAttrs == NULL) {
    injFifoAttrs = &defaultAttrs;
  }

  // Malloc space for the info structure
  msg_InjFifoInfo_t *info;
  info = (msg_InjFifoInfo_t *) memalign(32, sizeof(msg_InjFifoInfo_t));
  if( !info ) return -1;
  
    // Initialize the info structure
  info->startingSubgroupId = startingSubgroupId;
  info->startingFifoId     = startingFifoId;
  info->numFifos           = numFifos;
  info->numSubgroups       = 0;

  // Malloc space for the injection fifos.  They are 64-byte aligned.
  for (unsigned int i = 0; i < numFifos; i++) {
    info->fifoPtr[i] = (uint64_t*)memalign(64, fifoSize);
    if ( !info->fifoPtr[i] ) return -1;
  }
  
  // Process one subgroup at a time.
  // - Allocate the fifos.
  // - Init the MU MMIO for the fifos.
  // - Activate the fifos.
  while ( numFifos > 0 ) {
    info->numSubgroups++;
    
    // startingFifoId is the starting fifo number relative to the
    // subgroup we are working on.
    // Determine endingFifoId, the ending fifo number relative to
    // the subgroup we are working on.
    endingFifoId = startingFifoId + numFifos-1;
    if ( endingFifoId > (BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP-1) ) {
      endingFifoId = BGQ_MU_NUM_INJ_FIFOS_PER_SUBGROUP-1;
    }
    numFifosInSubgroup = endingFifoId - startingFifoId + 1;
    info->numFifosInSubgroup[subgroupId] = numFifosInSubgroup;
    
    // Init structures for allocating the fifos...
    // - fifo Ids
    // - attributes
    for (unsigned int i = 0; i < numFifosInSubgroup; i++) {
      fifoIds[i] = startingFifoId + i;
      memcpy(&attrs[i], injFifoAttrs, sizeof(attrs[i]));
    }
    
    // Allocate the fifos
    rc = Kernel_AllocateInjFifos (subgroupId,
				  &info->subgroup[subgroupId], 
				  numFifosInSubgroup,
				  fifoIds,
				  attrs);
    if ( rc ) {
      printf("msg_InjFifoInit: Kernel_AllocateInjFifos failed with rc=%d\n",rc);
      return rc;
    }
    
    // Init the MU MMIO for the fifos.
    for (unsigned int i = 0; i < numFifosInSubgroup; i++) {
      Kernel_MemoryRegion_t memRegion;
      rc = Kernel_CreateMemoryRegion ( &memRegion,
				       info->fifoPtr[numFifos-i-1],
				       fifoSize );
      if ( rc ) {
	printf("msg_InjFifoInit: Kernel_CreateMemoryRegion failed with rc=%d\n",rc);
	return rc;
      }
  
      // initialise the Fifos
      rc = Kernel_InjFifoInit (&info->subgroup[subgroupId], 
			       fifoIds[i],
			       &memRegion,
			       (uint64_t)info->fifoPtr[numFifos-i-1] -
			       (uint64_t)memRegion.BaseVa,
			       fifoSize-1);    
      if ( rc ) {
	printf("msg_InjFifoInit: Kernel_InjFifoInit failed with rc=%d\n",rc);
	return rc;
      }
    }
    
    // Activate the fifos.
    rc = Kernel_InjFifoActivate (&info->subgroup[subgroupId],
				 numFifosInSubgroup,
				 fifoIds,
				 KERNEL_INJ_FIFO_ACTIVATE);
    if ( rc ) {
      printf("msg_InjFifoInit: Kernel_InjFifoActivate failed with rc=%d\n",rc);
      return rc;
    }
    
    startingFifoId = 0; // Next subgroup will start at fifo 0.
    
    subgroupId++;       // Next subgroup.
    numFifos -= numFifosInSubgroup;
  }
  
  injFifoHandlePtr->pOpaqueObject = (void *)info;
  return 0;
}


void global_barrier() {
  int rc = 0;
  uint64_t timeoutCycles = 60UL * 1600000000UL; // about 60 sec at 1.6 ghz
  rc = MUSPI_GIBarrierEnter ( &GIBarrier );
  if (rc) {
    printf("MUSPI_GIBarrierEnter failed returned rc = %d\n", rc);
    exit(1);
  }
  
  // Poll for completion of the barrier.
  rc = MUSPI_GIBarrierPollWithTimeout ( &GIBarrier, timeoutCycles);
  if( rc ) {
    printf("MUSPI_GIBarrierPollWithTimeout failed returned rc = %d\n", rc);
    DelayTimeBase (200000000000UL);
    exit(1);
  }
  return;
}
