#ifndef GMX_PME_MALLOC_FREE_H
#define GMX_PME_MALLOC_FREE_H

#include<stdio.h>
#include<stdlib.h>

#include "pme-type.h"

#define N 121687
#define GRID 72
#define GRIDP 75
#define ORDER 4


void *save_calloc(const char *name, const char *file, int line,
                  size_t nelem, size_t elsize);

#ifdef __cplusplus
#define snew(ptr, nelem) \
    gmx_snew_impl(#ptr, __FILE__, __LINE__, (ptr), (nelem))

template <typename T> static inline
void gmx_snew_impl(const char *name, const char *file, int line,
                   T * &ptr, size_t nelem)
{
    ptr = (T *)save_calloc(name, file, line, nelem, sizeof(T));
}

#else

#define snew(ptr, nelem) \
    (ptr) = save_calloc(#ptr, __FILE__, __LINE__, (nelem), sizeof(*(ptr)))

#endif

void save_free(const char *name, const char *file, int line, void *ptr);

#define sfree(ptr) save_free(#ptr, __FILE__, __LINE__, (ptr))


void init_recipbox(matrix rbox);
void init_nnxyz(int * &nnx, int * &nny, int * &nnz);
void init_pmegrid(pmegrid_t &grid);
void init_g2t(int ** &g2t);
void init_pmegrids(pmegrids_t &grids);
void init_coordinate(rvec * &x);
void init_idx(ivec * &idx);
void init_fractx(rvec * &fractx);
void init_coefficient(real * &coefficient);
void init_ind(int * &ind);
void init_splinedata(splinedata_t * &spline);
void init_atomcomm(pme_atomcomm_t &atc);
tgmx_pme_t *pme_init();
void pmegrids_destroy(pmegrids_t &grids);
void splinedata_destroy(splinedata_t * &spline);
void atomcomm_destroy(pme_atomcomm_t &atc);
int gmx_pme_destroy(tgmx_pme_t * &pmedata);

#endif
