//
// Friends of Friends halo finder
//
// This code is based on a serial FoF halo finder by
// University of Washington N-BODY SHOP
// http://www-hpcc.astro.washington.edu/tools/fof.html
//


/*
						FOF v1.1

			A Group Finder for N-body Simulations

					October 26, 1994
*/

//
// Memory allocation:
//   mem1: move_particles2_min buffer / kdtree, fifo, igrp, halo data
//   mem2: snapshot particles
//

// TODO: memory allocation assumes fof_linking_parameter= 0.2
//       possible shortage of memory for larger l >~ 1.0

// TOTO: NHaloAlloc is resolution dependent. Possble lack of memory for
//       high resolution (with lots of small haloes)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "particle.h"
#include "msg.h"
#include "comm.h"
#include "timer.h"
#include "domain.h"

static const int nBucket= 16;

#define ROOT	   1
#define LOWER(i)   (i<<1)
#define UPPER(i)   ((i<<1)+1)
#define PARENT(i)  (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)
#define SETNEXT(i) { while (i&1) i=i>>1; ++i; }

#define LEFT  2
#define RIGHT 1
#define DUAL  3
#define GLOBAL 4

typedef struct bndBound {
  float fMin[3];
  float fMax[3];
} BND;

typedef struct kdNode {
  float fSplit;
  BND bnd;
  int iDim;
  int pLower;
  int pUpper;
} KDN;

typedef struct kdContext {
  int nBucket;
  int nActive;
  float fPeriod[3];
  int nLevels;
  int nNodes;
  int nSplit;
  ParticleMinimum *p;
  int* iGroup;
  KDN *kdNodes;
  int nGroup;
} *KD;

#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
  float dx,dy,dz,dx1,dy1,dz1,fDist2,fMax2;	\
  dx = c[cp].bnd.fMin[0]-x;			\
  dx1 = x-c[cp].bnd.fMax[0];			\
  dy = c[cp].bnd.fMin[1]-y;			\
  dy1 = y-c[cp].bnd.fMax[1];			\
  dz = c[cp].bnd.fMin[2]-z;			\
  dz1 = z-c[cp].bnd.fMax[2];			\
  if (dx > 0.0) {				\
    if (dx1+lx < dx) {				\
      dx1 += lx;				\
      dx -= lx;					\
      sx = x+lx;				\
      fDist2 = dx1*dx1;				\
      fMax2 = dx*dx;				\
    }						\
    else {					\
      sx = x;					\
      fDist2 = dx*dx;				\
      fMax2 = dx1*dx1;				\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dx1 > 0.0) {					\
    if (dx+lx < dx1) {					\
      dx += lx;						\
      dx1 -= lx;					\
      sx = x-lx;					\
      fDist2 = dx*dx;					\
      fMax2 = dx1*dx1;					\
    }							\
    else {						\
      sx = x;						\
      fDist2 = dx1*dx1;					\
      fMax2 = dx*dx;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sx = x;						\
    fDist2 = 0.0;					\
    if (dx < dx1) fMax2 = dx*dx;			\
    else fMax2 = dx1*dx1;				\
  }							\
  if (dy > 0.0) {					\
    if (dy1+ly < dy) {					\
      dy1 += ly;					\
      dy -= ly;						\
      sy = y+ly;					\
      fDist2 += dy1*dy1;				\
      fMax2 += dy*dy;					\
    }							\
    else {						\
      sy = y;						\
      fDist2 += dy*dy;					\
      fMax2 += dy1*dy1;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dy1 > 0.0) {					\
    if (dy+ly < dy1) {					\
      dy += ly;						\
      dy1 -= ly;					\
      sy = y-ly;					\
      fDist2 += dy*dy;					\
      fMax2 += dy1*dy1;					\
    }							\
    else {						\
      sy = y;						\
      fDist2 += dy1*dy1;				\
      fMax2 += dy*dy;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sy = y;						\
    if (dy < dy1) fMax2 += dy*dy;			\
    else fMax2 += dy1*dy1;				\
  }							\
  if (dz > 0.0) {					\
    if (dz1+lz < dz) {					\
      dz1 += lz;					\
      dz -= lz;						\
      sz = z+lz;					\
      fDist2 += dz1*dz1;				\
      fMax2 += dz*dz;					\
    }							\
    else {						\
      sz = z;						\
      fDist2 += dz*dz;					\
      fMax2 += dz1*dz1;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dz1 > 0.0) {					\
    if (dz+lz < dz1) {					\
      dz += lz;						\
      dz1 -= lz;					\
      sz = z-lz;					\
      fDist2 += dz*dz;					\
      fMax2 += dz1*dz1;					\
    }							\
    else {						\
      sz = z;						\
      fDist2 += dz1*dz1;				\
      fMax2 += dz*dz;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sz = z;						\
    if (dz < dz1) fMax2 += dz*dz;			\
    else fMax2 += dz1*dz1;				\
  }							\
  if (fMax2 < fBall2) goto ContainedCell;		\
}

void kdTime(KD,int *,int *);
void kdBuildTree(KD);
int kdFoF(KD,float);
//int kdTooSmall(KD,int);
//void kdFinish(KD);

//
// Added part for parallel version
//

typedef struct {
  int nfof;
  int boundary;
  float x0[3], dx_sum[3];
  float v_sum[3];
  int merge_to;
} HaloInfo;

static float BoxSize, HalfBoxSize;
static float xleft, xright;

static KDN* KdNodes;
static int  NKdNodeAlloc;
static int *igrp, *Fifo, *Map;
static HaloInfo* Halo;
static int NHalo, NHaloAlloc, NHaloBuf;
static const int nfof_min= 32;

// Particles exported to next Node
static float *dx_buf, *dx_recv_buf;
static int* igrp_buf;
static int np_export, np_export_alloc;
static HaloInfo* halo_buf;
static int nhalo_export, nhalo_export_alloc;

static void* Buf; static size_t BufSize;

static int* GlobalLinking= 0;
static int NGlobalLinking, NGlobalLinking1;

static inline void export_particle(ParticleMinimum const * const p)
{
  // Add particle p to the export buffer dx_buf, igrp_buf
  if(np_export >= np_export_alloc)
    msg_abort(7200, "Error: Not enough space for fof particle export: %d\n", 
	      np_export);

  dx_buf[3*np_export    ]= p->x[0];
  dx_buf[3*np_export + 1]= p->x[1];
  dx_buf[3*np_export + 2]= p->x[2];
  igrp_buf[np_export]= nhalo_export;

  np_export++;
}


static inline void export_halo(HaloInfo const * const h)
{
  // Add halo h to export buffer halo_buf
  if(nhalo_export >= nhalo_export_alloc)
    msg_abort(7300, "Error: not enough space for halo export: %d\n", 
	      nhalo_export);

  halo_buf[nhalo_export]= *h;

  nhalo_export++;
}

static inline void export_empty_halo(HaloInfo const * const h)
{
  // Copy halo h to export buffer halo_buf (but nfof=0)/
  // Halo information remains in the current node
  if(nhalo_export >= nhalo_export_alloc)
    msg_abort(7301, "Error: not enough space for halo export: %d\n", 
	      nhalo_export);

  HaloInfo* const h0= halo_buf + nhalo_export;
  h0->nfof= 0;
  h0->boundary= h->boundary;
  for(int i=0; i<3; i++) {
    h0->x0[i]= h->x0[i];
    h0->dx_sum[i]= 0.0f;
    h0->v_sum[i]= 0.0f;
  }

  nhalo_export++;
}


static inline void halo_particle(HaloInfo* const h, ParticleMinimum const * const p)
{
  if(NHalo >= NHaloAlloc)
    msg_abort(7050, "Not enough space for halo\n");

  // for each particle belongs to halo h
  if(p->x[0] > xright) {
    export_particle(p);
    h->boundary |= RIGHT;
  }
  else if(p->x[0] < xleft) {
    h->boundary |= LEFT;
  }

  if(h->nfof == 0) {
    for(int i=0; i<3; i++)
      h->x0[i]= p->x[i];
  }
  else {
    for(int i=0; i<3; i++) {
      float dx= p->x[i] - h->x0[i];
      if(dx > HalfBoxSize)
	dx -= BoxSize;
      else if(dx < -HalfBoxSize)
	dx += BoxSize;

      h->dx_sum[i] += dx;
    }
  }

  for(int i=0; i<3; i++)
    h->v_sum[i] += p->v[i];

  

  h->nfof++;
}
  
static inline void halo_clear(HaloInfo* const h)
{
  h->nfof= 0;
  h->boundary= 0;
  h->merge_to= 0;
  for(int i=0; i<3; i++) {
    h->x0[i]= 0.0f;
    h->dx_sum[i]= 0.0f;
    h->v_sum[i]= 0.0f;
  }
}

static inline void merge_halo(int i, int j)
{
  // Merge halo i to j
  // Algorithm of finding equivalent class. See Numerical Recepie
  while(Halo[i].merge_to != i)
    i= Halo[i].merge_to;

  while(Halo[j].merge_to != j)
    j= Halo[j].merge_to;

  assert(0 <= i && i < NHaloBuf);
  assert(0 <= j && j < NHalo); // ** slow check

  //msg_printf("halo merge %d -> %d\n", i, j);
  //printf("halo merge %d -> %d\n", i, j);

  Halo[i].merge_to= j;
}


void kdSelect(KD kd,int d,int k,int l,int r)
{
  ParticleMinimum* p = kd->p;
  while(r > l) {
    double v = p[k].x[d];
    ParticleMinimum t = p[r];
    p[r] = p[k];
    p[k] = t;
    int i = l - 1;
    int j = r;
    while(1) {
      while(i < j) if (p[++i].x[d] >= v) break;
      while(i < j) if (p[--j].x[d] <= v) break;
      t = p[i];
      p[i] = p[j];
      p[j] = t;
      if (j <= i) break;
    }
    p[j] = p[i];
    p[i] = p[r];
    p[r] = t;
    if (i >= k) r = i - 1;
    if (i <= k) l = i + 1;
  }
}


void kdCombine(KDN const * const p1, KDN const * const p2, KDN * const pOut)
{
  // Combine the bounds.
  for (int j=0; j<3; j++) {
    if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
      pOut->bnd.fMin[j] = p2->bnd.fMin[j];
    else
      pOut->bnd.fMin[j] = p1->bnd.fMin[j];
    if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
      pOut->bnd.fMax[j] = p2->bnd.fMax[j];
    else
      pOut->bnd.fMax[j] = p1->bnd.fMax[j];
  }
}


void kdUpPass(KD kd, int iCell)
{
  int l, u, pj,j;

  KDN* c = kd->kdNodes;
  if (c[iCell].iDim != -1) {
    l = LOWER(iCell);
    u = UPPER(iCell);
    kdUpPass(kd, l);
    kdUpPass(kd, u);
    kdCombine(&c[l], &c[u], &c[iCell]);
  }
  else {
    l = c[iCell].pLower;
    u = c[iCell].pUpper;
    for (j=0;j<3;++j) {
      c[iCell].bnd.fMin[j] = kd->p[u].x[j];
      c[iCell].bnd.fMax[j] = kd->p[u].x[j];
    }
    for (pj=l;pj<u;++pj) {
      for (j=0;j<3;++j) {
	if (kd->p[pj].x[j] < c[iCell].bnd.fMin[j])
	  c[iCell].bnd.fMin[j] = kd->p[pj].x[j];
	if (kd->p[pj].x[j] > c[iCell].bnd.fMax[j])
	  c[iCell].bnd.fMax[j] = kd->p[pj].x[j];
      }
    }
  }
}

void kdBuildTree(KD kd)
{
  BND bnd;
  
  int n = kd->nActive;
  kd->nLevels = 1;
  int l = 1;
  while(n > kd->nBucket) {
    n = n>>1;
    l = l<<1;
    ++kd->nLevels;
  }
  kd->nSplit = l;
  kd->nNodes = l<<1;
  kd->kdNodes= KdNodes;
  assert(NKdNodeAlloc >= kd->nNodes);
  assert(kd->kdNodes != NULL);

  // Calculate Bounds.
  for (int j=0; j<3; ++j) {
    bnd.fMin[j] = kd->p[0].x[j];
    bnd.fMax[j] = kd->p[0].x[j];
  }
  for (int i=1; i<kd->nActive; ++i) {
    for (int j=0; j<3; ++j) {
      if (bnd.fMin[j] > kd->p[i].x[j]) 
	bnd.fMin[j] = kd->p[i].x[j];
      else if (bnd.fMax[j] < kd->p[i].x[j])
	bnd.fMax[j] = kd->p[i].x[j];
    }
  }

  // Set up ROOT node

  KDN* c = kd->kdNodes;
  c[ROOT].pLower = 0;
  c[ROOT].pUpper = kd->nActive-1;
  c[ROOT].bnd = bnd;
  int i = ROOT;
  while(1) {
    assert(c[i].pUpper - c[i].pLower + 1 > 0);
    if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
      int d = 0;
      for (int j=1; j<3; ++j) {
	if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
	    c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
      }
      c[i].iDim = d;
      
      int m = (c[i].pLower + c[i].pUpper)/2;
      kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);
      
      c[i].fSplit = kd->p[m].x[d];
      c[LOWER(i)].bnd = c[i].bnd;
      c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
      c[LOWER(i)].pLower = c[i].pLower;
      c[LOWER(i)].pUpper = m;
      c[UPPER(i)].bnd = c[i].bnd;
      c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
      c[UPPER(i)].pLower = m+1;
      c[UPPER(i)].pUpper = c[i].pUpper;
      int diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
      assert(diff == 0 || diff == 1);
      i = LOWER(i);
    }
    else {
      c[i].iDim = -1;
      SETNEXT(i);
      if (i == ROOT) break;
    }
  }
  kdUpPass(kd,ROOT);
}


int kdFoF(KD kd,float fEps)
{
  float dx,dy,dz,x,y,z,sx,sy,sz,fDist2;
  
  ParticleMinimum* const p = kd->p;
  KDN* const c= kd->kdNodes;
  const float lx = kd->fPeriod[0]*2; // No need for periodic wrapup in x
  const float ly = kd->fPeriod[1];
  const float lz = kd->fPeriod[2];
  const float fEps2 = fEps*fEps;

  const int nFifo = kd->nActive;


  int iHead= 0;
  int iTail= 0;
  int iGroup= 0;
  
  assert(Map);

  NHalo= 0;
  np_export= 0;
  nhalo_export= 0;

  assert(igrp);
  for (int pn=0; pn<kd->nActive; ++pn)
    igrp[pn]= 0;

  for (int pn=0; pn<kd->nActive; ++pn) {
    if (igrp[pn]) continue;

    ++iGroup; 

    halo_clear(&Halo[NHalo]);

    //
    // Mark it and add to the do-fifo.
    //
    igrp[pn]= iGroup; halo_particle(&Halo[NHalo], &p[pn]);
    Fifo[iTail++] = pn;
    if (iTail == nFifo) iTail = 0;

    while (iHead != iTail) {
      int pi = Fifo[iHead++];
      if (iHead == nFifo) iHead = 0;

      //
      // Now do an fEps-Ball Gather!
      //
      x = p[pi].x[0];
      y = p[pi].x[1];
      z = p[pi].x[2];
      int cp = ROOT;
      while (1) {
	INTERSECT(c,cp,fEps2,lx,ly,lz,x,y,z,sx,sy,sz);
	//
	// We have an intersection to test.
	//
	if (c[cp].iDim >= 0) {
	  cp = LOWER(cp);
	  continue;
	}
	else {
	  for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	    if (igrp[pj]) continue;
	    dx = sx - p[pj].x[0];
	    dy = sy - p[pj].x[1];
	    dz = sz - p[pj].x[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < fEps2) {
	      //
	      // Mark it and add to the do-fifo.
	      //
	      igrp[pj]= iGroup; halo_particle(&Halo[NHalo], &p[pj]);
	      Fifo[iTail++] = pj;
	      if (iTail == nFifo) iTail = 0;
	    }
	  }
	  SETNEXT(cp);
	  if (cp == ROOT) break;
	  continue;
	}
      ContainedCell:
	for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	  if (igrp[pj]) continue;
	  //
	  // Mark it and add to the do-fifo.
	  //
	  igrp[pj] = iGroup; halo_particle(&Halo[NHalo], &p[pj]);
	  Fifo[iTail++] = pj;
	  if (iTail == nFifo) iTail = 0;
	}
      GetNextCell:
	SETNEXT(cp);
	if (cp == ROOT) break;
      }
    }

    if(Halo[NHalo].boundary == 3) {
      export_empty_halo(&Halo[NHalo]);
      Map[iGroup]= NHalo;
      NHalo++; 
    }
    else if(Halo[NHalo].boundary == RIGHT) {
      // Halo exported to the right node
      export_halo(&Halo[NHalo]);
      Map[iGroup]= -2; // exported
    }
    else if(Halo[NHalo].nfof >= nfof_min || 
       (Halo[NHalo].boundary == LEFT)) {
      Map[iGroup]= NHalo;
      NHalo++; 
    }
    else
      Map[iGroup]= -1; // too small
  }

  // remap igrp
  for (int pn=0; pn<kd->nActive; ++pn)
    igrp[pn]= Map[igrp[pn]];

  kd->nGroup = NHalo; // *** +plus1 ? well not used anymore...anyway
  

  return NHalo;
}


/*
void kdFinish(KD kd)
{
  //free(kd->p);
  //free(kd->kdNodes);
  //free(kd);
}
*/

//
// MPI Communication
//


// send particle positions for FOF linking
static int fof_send_buffer_positions(Snapshot* const snapshot, int* igrp) 
{
  ParticleMinimum* const p= snapshot->p;
  const int np_local= snapshot->np_local;
  const int np_alloc= snapshot->np_allocated;

  msg_printf(verbose, "Exchanging FOF buffer positions.\n");

  if(comm_right_edge()) {
    for(int i=0; i<np_export; i++)
      dx_buf[3*i] -= BoxSize;
  }

  int nrecv= comm_get_nrecv(ToRight, np_export);
  if(np_local + nrecv > np_alloc)
    msg_abort(7100, "Error: Not enough space for FOF buffer particles: "
	            "%d + %d particles ", np_local, nrecv);

  //msg_printf("FOF buffer particles %d %d\n", np_export, nrecv);
  comm_sendrecv(ToRight, dx_buf, 3*np_export, dx_recv_buf, 3*nrecv, MPI_FLOAT);

  // Buffer particle positions added as particles np_local <= i 
  for(int i=0; i<nrecv; ++i) {
    p[np_local + i].id  = -1;         // buffer particles have positions only
    p[np_local + i].x[0]= dx_recv_buf[3*i  ];
    p[np_local + i].x[1]= dx_recv_buf[3*i+1];
    p[np_local + i].x[2]= dx_recv_buf[3*i+2];    
  }

  comm_sendrecv(ToRight, igrp_buf, np_export, igrp+np_local, nrecv, MPI_INT);
  for(int i=0; i<nrecv; ++i)
    igrp[i+np_local] += NHalo;
  
  int nsend_global;
  MPI_Reduce(&np_export, &nsend_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  msg_printf(info, "%d particle positions copied for FOF.\n", nsend_global);

  return np_local + nrecv;
}

static int fof_send_halo()
{
  //MPI_Barrier(MPI_COMM_WORLD); // !!!debug necessary?
  int nrecv= comm_get_nrecv(ToRight, nhalo_export);
  if(NHalo + nrecv > NHaloAlloc)
    msg_abort(7400, "Error: Not enough space for fof_send_halo(): "
	            "%d + %d halos\n", NHalo, nrecv);

  // Buffer Halo info added at the end of "halo"
  comm_sendrecv(ToRight, halo_buf, nhalo_export*sizeof(HaloInfo), 
		Halo+NHalo, nrecv*sizeof(HaloInfo), MPI_BYTE);

  msg_printf(debug, "Number of halos (local) %d + %d\n", NHalo, nrecv);

  NHaloBuf= NHalo + nrecv;

  return NHalo + nrecv;
}

//
// Link buffer particles and merge halo info across the node boundary
//
int link_buffer_particles(KD kd, float fEps, int nhalo, int np_local, int np_plus_buffer)
{
  float dx,dy,dz,x,y,z,sx,sy,sz,fDist2;
  
  ParticleMinimum* const p = kd->p;
  KDN* const c= kd->kdNodes;
  const float lx = kd->fPeriod[0]*2;
  const float ly = kd->fPeriod[1];
  const float lz = kd->fPeriod[2];
  const float fEps2 = fEps*fEps;

  for(int i=0; i<nhalo; i++)
    Halo[i].merge_to= i;

  for (int pn=np_local; pn<np_plus_buffer; ++pn) {
    // Now do an fEps-Ball Gather!
    x = p[pn].x[0];
    y = p[pn].x[1];
    z = p[pn].x[2];

    int cp = ROOT;
    while(1) {
      INTERSECT(c,cp,fEps2,lx,ly,lz,x,y,z,sx,sy,sz);
      // We have an intersection to test.
      if (c[cp].iDim >= 0) {
	cp = LOWER(cp);
	continue;
      }
      else {
	for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	  dx = sx - p[pj].x[0];
	  dy = sy - p[pj].x[1];
	  dz = sz - p[pj].x[2];
	  fDist2 = dx*dx + dy*dy + dz*dz;
	  if (fDist2 < fEps2) {
	    merge_halo(igrp[pn], igrp[pj]);
	  }
	}
	SETNEXT(cp);
	if (cp == ROOT) break;
	continue;
      }
    ContainedCell:
      for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	merge_halo(igrp[pn], igrp[pj]);
      }
    GetNextCell:
      SETNEXT(cp);
      if (cp == ROOT) break;
    }
  }
  
  return NHalo;
}

int top_merge(int ihalo)
{
  while(Halo[ihalo].merge_to != ihalo) ihalo= Halo[ihalo].merge_to;
  assert(0 <= ihalo && ihalo < NHaloBuf);
  return ihalo;
}

void merge_halo_info(HaloInfo* const h, const int nhalo)
{
  for(int i=0; i<nhalo; i++) {
    if(h[i].merge_to == i) continue;

    int top= i;
    while(h[top].merge_to != top) top= h[top].merge_to;
    assert(0 <= top && top < nhalo);

    for(int k=0; k<3; k++) {
      float dx0= h[i].x0[k] - h[top].x0[k];
      if(dx0 <= -HalfBoxSize) dx0 += BoxSize;
      if(dx0 > HalfBoxSize) dx0 -= BoxSize;

      h[top].dx_sum[k] += h[i].dx_sum[k] + dx0*h[i].nfof;
      h[top].v_sum[k]  += h[i].v_sum[k];
    }

    h[top].nfof += h[i].nfof;
    h[i].nfof= 0;
  }
}

int delete_small_halos(const int nhalo)
{
  int copyto=0;

  for(int i=0; i<nhalo; i++) {
    if(Halo[i].nfof >= nfof_min ||
       (Halo[i].boundary & GLOBAL) == GLOBAL) {
      Map[i]= copyto;
      Halo[copyto++]= Halo[i];
    }
    else
      Map[i]= -1;
  }

  NHalo= copyto;

  return copyto;
}

void setup_global_linking(const int nhalo)
{
  // Setup linking information of
  // Boundary == 3 halos (which touches both left and right boundary)
  // linking information (top of merge_to linking) is collected

  int nboundary3= 0;

  for(int i=0; i<NHalo; i++)
    if(Halo[i].boundary == 3) nboundary3++;

  NGlobalLinking1= nboundary3;

  for(int i=NHalo; i<nhalo; i++) 
    if(Halo[i].boundary == 3) nboundary3++;

  NGlobalLinking= nboundary3;

  if(GlobalLinking) free(GlobalLinking);
  GlobalLinking= malloc(sizeof(int)*nboundary3); assert(GlobalLinking);


  int n= 0;
  for(int i=0; i<nhalo; i++) {
    if((Halo[i].boundary & DUAL) == DUAL) {
      int itop= top_merge(i);
      GlobalLinking[n]= itop;
      Halo[itop].boundary |= GLOBAL;
      n++;
    }
  }

  assert(n == NGlobalLinking);

  msg_printf(normal, "Set %d global linking information.\n", NGlobalLinking);
}

void remap_global_linking()
{
  for(int i=0; i<NGlobalLinking; i++) {
    int inew= Map[GlobalLinking[i]];
    assert(0 <= inew && inew < NHalo);
    GlobalLinking[i]= Map[GlobalLinking[i]];
  }
}

void global_linking(HaloInfo* h, const int nhalo_total)
{
  // Link same boundary=3 halo, original one to exported one
  msg_printf(verbose, "global_linking. nhalo_total= %d\n", nhalo_total);

  const int this_node= comm_this_node();
  const int nnode= comm_nnode();
  int *nglobal_linking= 0, *disp= 0, *ihalo_offset= 0, *nhalo= 0;
  

  if(this_node == 0) {
    nglobal_linking= malloc(sizeof(int)*nnode*4); assert(nglobal_linking);
    disp= nglobal_linking + nnode;
    ihalo_offset= disp + nnode;
    nhalo= ihalo_offset + nnode;
  }


  // Gather GlobalLiking[]

  int ret= 
  MPI_Gather(&NGlobalLinking, 1, MPI_INT, 
	     nglobal_linking, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  int* global_linking= 0;
  int nglobal_linking_tot= 0;

  if(this_node == 0) {
    int offset= 0;

    for(int i=0; i<nnode; i++) {
      disp[i]= offset;
      offset += nglobal_linking[i];
    }

    global_linking= malloc(sizeof(int)*offset); assert(global_linking);

    nglobal_linking_tot= offset;
  }


  ret= MPI_Gatherv(GlobalLinking, NGlobalLinking, MPI_INT,
	      global_linking, nglobal_linking, disp, MPI_INT,
	      0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);
  msg_printf(normal, "Global linking information %d.\n", nglobal_linking_tot);

  ret= 
    MPI_Gather(&NHalo, 1, MPI_INT, nhalo, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  //
  if(this_node == 0) {
    int offset= 0;
    for(int i=0; i<nnode; i++) {
      ihalo_offset[i]= offset;
      offset += nhalo[i];
    }

    // initialize merge_to
    for(int i=0; i<nhalo_total; i++)
      h[i].merge_to= i;
    
    // Linking
    int nboundary3= NGlobalLinking1;
    int count= 0;
    for(int i=0; i<nnode; i++) {
      int inode= comm_node(i);
      int rnode= comm_node(i+1); // right of node i
      
      int ioffset= disp[inode];
      int roffset= disp[rnode] + nglobal_linking[rnode] - nboundary3;
      
      for(int k=0; k<nboundary3; k++) {
	// Link each boundary=3 haloes
	assert(0 <= ioffset+k && ioffset+k < nglobal_linking_tot);
	assert(0 <= roffset+k && roffset+k < nglobal_linking_tot);

	int ihalo= ihalo_offset[inode] + global_linking[ioffset+k];
	int rhalo= ihalo_offset[rnode] + global_linking[roffset+k];
	
	assert(0 <= ihalo && ihalo < nhalo_total);
	assert(0 <= rhalo && rhalo < nhalo_total);
	
	// merge ihalo to rhalo
	while(h[ihalo].merge_to != ihalo)
	  ihalo= h[ihalo].merge_to; // top of ihalo
	
	while(h[rhalo].merge_to != rhalo)
	  rhalo= h[rhalo].merge_to;
	
	assert(0 <= ihalo && ihalo < nhalo_total);
	assert(0 <= rhalo && rhalo < nhalo_total);

	h[ihalo].merge_to= rhalo;
	count++;
      }
      nboundary3= nglobal_linking[rnode] - nboundary3;
    }
    
    assert(nboundary3 == NGlobalLinking1);
    assert(2*count == nglobal_linking_tot);

    free(nglobal_linking);
    free(global_linking);

    printf("global_linking done\n");
  }

}

//
// Main interface called from main.c
//

static int n_kd_nodes(int n)
{
  int l = 1;
  while(n > nBucket) {
    n= n >> 1;
    l= l << 1;
  }
  return l<<1;
}

size_t fof_calc_memory(const int np_alloc, const int nc)
{
  int nNodes= n_kd_nodes(np_alloc);
  size_t size= sizeof(KDN)*nNodes;

  size += sizeof(int)*np_alloc*3; // igrp, Map, Fifo

  int n_halo_alloc= 0.2*np_alloc;
  int nhalo_export_alloc= nc*nc;
  size += sizeof(HaloInfo)*(n_halo_alloc + nhalo_export_alloc);

  int np_export_alloc= nc*nc;
  size += (sizeof(float)*6 + sizeof(int))*np_export_alloc;

  return size;
}

void fof_init(const int np_alloc, const int nc, void* mem, size_t mem_size)
{
  size_t bytes= 0;

  Buf= mem;
  BufSize= mem_size;

  NHaloAlloc= 0.2*np_alloc; // depends on collapsed fraction (set to 0.2)
  Halo= (HaloInfo*) mem; mem= Halo + NHaloAlloc;
  bytes += sizeof(HaloInfo)*NHaloAlloc;

  np_export_alloc= nc*nc; // mean is linking_parameter (0.2)*nc*nc;
  
  dx_buf= mem; mem= dx_buf + 3*np_export_alloc;

  dx_recv_buf= mem; mem= dx_recv_buf + 3*np_export_alloc;
  bytes += sizeof(float)*6*np_export_alloc;

  igrp_buf= mem; mem= igrp_buf + np_export_alloc;
  bytes += sizeof(int)*np_export_alloc;

  nhalo_export_alloc= nc*nc;

  halo_buf= mem; mem= halo_buf + nhalo_export_alloc;
  bytes += sizeof(HaloInfo)*nhalo_export_alloc;


  int nNodes= n_kd_nodes(np_alloc);
  KdNodes= (KDN*) mem; mem= KdNodes + nNodes;
  NKdNodeAlloc= nNodes;  
  bytes += sizeof(KDN)*nNodes;
  
  igrp= (int*) mem;        
  Fifo= igrp + np_alloc;
  Map= Fifo + np_alloc;

  mem= Map + np_alloc; bytes += sizeof(int)*3*np_alloc;


  msg_printf(info, "%d Mbytes allocated for FOF halo finding within mem1\n", bytes/(1024*1024));

  assert(bytes <= mem_size);
}

void fof_find_halos(Snapshot* snapshot, const float ll)
{
  xleft= comm_xmin() + ll;
  xright= comm_xmax() - ll;
  BoxSize= snapshot->boxsize; assert(BoxSize > 0.0f);
  HalfBoxSize= 0.5f*BoxSize;

  msg_printf(verbose, "FOF halo finding...\n");
                                                         timer_start(comm);
  domain_wrap_min(snapshot, Buf, BufSize);
  domain_decompose_min(snapshot, Buf, BufSize);
                                                         timer_stop(comm);

  struct kdContext kdcontext;
  KD kd= &kdcontext;
  kd->nBucket= nBucket;
  kd->p= snapshot->p;
  kd->iGroup= igrp;
  kd->kdNodes= 0;
  kd->nActive= snapshot->np_local;
  for (int j=0; j<3; j++) 
    kd->fPeriod[j] = BoxSize;

                                                         timer_start(kd_build);
  kdBuildTree(kd);
  msg_printf(verbose, "KD tree built.\n");
                                                         timer_stop(kd_build);
							 timer_start(kd_link);
  int nGroup= kdFoF(kd, ll);
  msg_printf(debug, "Number of initial groups:%d\n", nGroup);
                                                         timer_stop(kd_link);

                                                         timer_start(comm);
  int np_with_buffer= fof_send_buffer_positions(snapshot, igrp);
  int nhalo_with_buffer= fof_send_halo();
  

  link_buffer_particles(kd, ll, nhalo_with_buffer, 
			snapshot->np_local, np_with_buffer);

  merge_halo_info(Halo, nhalo_with_buffer);
                                                         timer_stop(comm);
                                                         timer_start(global);
  timer_start(global);
  setup_global_linking(nhalo_with_buffer);
  
  delete_small_halos(nhalo_with_buffer);
  
  remap_global_linking();
                                                         timer_stop(global);
}

void fof_write_halos(char filename[])
{
                                                         timer_start(global);
  // velocity is internal (no LPT correction and scaling)
  const int this_node= comm_this_node();
  const int nnode= comm_nnode();

  int* const nhalo_recv= malloc(sizeof(int)*nnode*2);
  int* const disp= nhalo_recv + nnode;

  int ret= 
    MPI_Gather(&NHalo, 1, MPI_INT, nhalo_recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  int nhalo=0;
  int new_allocation= 0;

  HaloInfo* buf1= 0;

  if(this_node == 0) {
    for(int i=0; i<nnode; i++) {
      disp[i] = nhalo*sizeof(HaloInfo);
      nhalo += nhalo_recv[i];
      nhalo_recv[i] *= sizeof(HaloInfo);
    }
    
    msg_printf(info, "nhalo %d, %d Mbytes gathered to node 0\n", nhalo,
	       (int)(sizeof(HaloInfo)*nhalo/(1024*1024)));


    if(BufSize >= sizeof(HaloInfo)*(NHalo + nhalo)) {
      buf1= Halo + NHalo;
    }
    else {
      msg_printf(info, "Newly allocating %d Mbytes for all halos\n",
		 (int)(sizeof(HaloInfo)*nhalo/(1024*1024)));
      buf1= malloc(sizeof(HaloInfo)*nhalo);
      new_allocation= 1;
    }
      
    if(buf1 == 0) {
      msg_abort(8100, "Error: Unable to allocate memory for halo output: %d Mbytes\n", sizeof(HaloInfo)*nhalo/(1024*1024));
    }
  }

  ret= MPI_Gatherv(Halo, NHalo*sizeof(HaloInfo), MPI_BYTE,
	      buf1, nhalo_recv, disp, MPI_BYTE,
	      0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  // "Global Linking" merge haloes extends more than 1 node
  HaloInfo* h= (HaloInfo*) buf1;
  global_linking(h, nhalo); //, nnode, nhalo_recv);
  merge_halo_info(h, nhalo);

  int nhalo_written= 0;

  if(this_node == 0) {
    FILE* fp= fopen(filename, "w");
    if(fp == 0)
      msg_abort(8110, "Error: Unable to write halo to file %s\n", filename);
    
    for(int i=0; i<nhalo; ++i) {
      if(h[i].nfof < nfof_min) continue;

      nhalo_written++;

      for(int k=0; k<3; k++) {
	float x= h[i].x0[k] + h[i].dx_sum[k]/h[i].nfof;
	if(x < 0.0f) x += BoxSize;
	if(x >= BoxSize) x -= BoxSize; 
	// not elseif because x += BoxSize could be equal to Boxsize for 
	// x<0 very close to zero due to finite precision

	// fof.txt 0 <= x < boxsize guaranteed
	h[i].x0[k]= x;
	
	h[i].v_sum[k]= h[i].v_sum[k]/h[i].nfof;
      }
      
      fprintf(fp, "%d %e %e %e %e %e %e\n", 
	      h[i].nfof, h[i].x0[0], h[i].x0[1], h[i].x0[2],
	      h[i].v_sum[0], h[i].v_sum[1], h[i].v_sum[2]);
    }
    fclose(fp);

    msg_printf(normal, "%d halos written to %s\n", nhalo_written, filename);
  }


  free(nhalo_recv);
  if(new_allocation)
    free(buf1);
                                                             timer_stop(global);
}
  
