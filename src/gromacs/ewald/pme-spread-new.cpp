/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "pme-spread.h"

#include "config.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/simd/simd.h"

#include "pme-type.h"

#include "pme-simd.h"

/*  #define WRITE_DATA */

#define SUB_GRID_X 4
#define SUB_GRID_Y 8
#define SUB_GRID_Z 8

/* TODO consider split of pme-spline from this file */

static void calc_interpolation_idx(struct gmx_pme_t *pme, pme_atomcomm_t *atc,
                                   int start, int grid_index, int end, int thread)
{
    int             i;
    int            *idxptr, tix, tiy, tiz;
    real           *xptr, *fptr, tx, ty, tz;
    real            rxx, ryx, ryy, rzx, rzy, rzz;
    int             nx, ny, nz;
    int            *thread_idx = NULL;
    int            *tpl_n      = NULL;
    int             thread_i;

    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    for (i = start; i < end; i++)
    {
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];

        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );

        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fptr[XX] = tx - tix + pme->fshx[tix];
        fptr[YY] = ty - tiy + pme->fshy[tiy];
        fptr[ZZ] = tz - tiz;

        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];
    } // end of for
}

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE(order)                     \
    {                                              \
        for (int j = 0; (j < DIM); j++)            \
        {                                          \
            real dr, div;                          \
            real data[PME_ORDER_MAX];              \
                                                   \
            dr  = xptr[j];                         \
                                               \
            /* dr is relative offset from lower cell limit */ \
            data[order-1] = 0;                     \
            data[1]       = dr;                          \
            data[0]       = 1 - dr;                      \
                                               \
            for (int k = 3; (k < order); k++)      \
            {                                      \
                div       = 1.0/(k - 1.0);               \
                data[k-1] = div*dr*data[k-2];      \
                for (int l = 1; (l < (k-1)); l++)  \
                {                                  \
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)* \
                                       data[k-l-1]);                \
                }                                  \
                data[0] = div*(1-dr)*data[0];      \
            }                                      \
            /* differentiate */                    \
            dtheta[j][i*order+0] = -data[0];       \
            for (int k = 1; (k < order); k++)      \
            {                                      \
                dtheta[j][i*order+k] = data[k-1] - data[k]; \
            }                                      \
                                               \
            div           = 1.0/(order - 1);                 \
            data[order-1] = div*dr*data[order-2];  \
            for (int l = 1; (l < (order-1)); l++)  \
            {                                      \
                data[order-l-1] = div*((dr+l)*data[order-l-2]+    \
                                       (order-l-dr)*data[order-l-1]); \
            }                                      \
            data[0] = div*(1 - dr)*data[0];        \
                                               \
            for (int k = 0; k < order; k++)        \
            {                                      \
                theta[j][i*order+k]  = data[k];    \
            }                                      \
        }                                          \
    }

static void make_bsplines(splinevec theta, splinevec dtheta, int order,
                          rvec fractx[], int nr, int ind[], real coefficient[],
                          gmx_bool bDoSplines)
{
    /* construct splines for local atoms */
    int   i, ii;
    real *xptr;

    for (i = 0; i < nr; i++)
    {
        /* With free energy we do not use the coefficient check.
         * In most cases this will be more efficient than calling make_bsplines
         * twice, since usually more than half the particles have non-zero coefficients.
         */
        ii = ind[i];
        if (bDoSplines || coefficient[ii] != 0.0)
        {
            xptr = fractx[ii];
            assert(order >= 4 && order <= PME_ORDER_MAX);
            switch (order)
            {
                case 4:  CALC_SPLINE(4);     break;
                case 5:  CALC_SPLINE(5);     break;
                default: CALC_SPLINE(order); break;
            }
        }
    }
}

/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE(order)                            \
    for (ithx = 0; (ithx < order); ithx++)                    \
    {                                                    \
        index_x = (i0+ithx)*pny*pnz;                     \
        valx    = coefficient*thx[ithx];                          \
                                                     \
        for (ithy = 0; (ithy < order); ithy++)                \
        {                                                \
            valxy    = valx*thy[ithy];                   \
            index_xy = index_x+(j0+ithy)*pnz;            \
                                                     \
            for (ithz = 0; (ithz < order); ithz++)            \
            {                                            \
                index_xyz        = index_xy+(k0+ithz);   \
                grid[index_xyz] += valxy*thz[ithz];      \
            }                                            \
        }                                                \
    }


static void spread_coefficients_bsplines_thread(pmegrid_t                         *pmegrid,
                                                pme_atomcomm_t                    *atc,
                                                splinedata_t                      *spline,
                                                struct pme_spline_work *work)
{

    /* spread coefficients from home atoms to local grid */
    real          *grid;
    int            i, nn, n, ithx, ithy, ithz, i0, j0, k0;
    int       *    idxptr;
    int            order, norder, index_x, index_xy, index_xyz;
    real           valx, valxy, coefficient;
    real          *thx, *thy, *thz;
    int            pnx, pny, pnz, ndatatot;
    int            offx, offy, offz;

#if defined PME_SIMD4_SPREAD_GATHER && !defined PME_SIMD4_UNALIGNED
    real           thz_buffer[GMX_SIMD4_WIDTH*3], *thz_aligned;

    thz_aligned = gmx_simd4_align_r(thz_buffer);
#endif

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx*pny*pnz;
    grid     = pmegrid->grid;
    for (i = 0; i < ndatatot; i++)
    {
        grid[i] = 0;
    }

    order = pmegrid->order;

    for (nn = 0; nn < spline->n; nn++)
    {
        n           = spline->ind[nn];
        coefficient = atc->coefficient[n];

        if (coefficient != 0)
        {
            idxptr = atc->idx[n];
            norder = nn*order;

            i0   = idxptr[XX] - offx;
            j0   = idxptr[YY] - offy;
            k0   = idxptr[ZZ] - offz;

            thx = spline->theta[XX] + norder;
            thy = spline->theta[YY] + norder;
            thz = spline->theta[ZZ] + norder;

            switch (order)
            {
                case 4:
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_SPREAD_SIMD4_ORDER4
#else
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme-simd4.h"
#else
                    DO_BSPLINE(4);
#endif
                    break;
                case 5:
#ifdef PME_SIMD4_SPREAD_GATHER
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 5
#include "pme-simd4.h"
#else
                    DO_BSPLINE(5);
#endif
                    break;
                default:
                    DO_BSPLINE(order);
                    break;
            }
        }
    }
}

void spread_on_grid(struct gmx_pme_t *pme,
                    pme_atomcomm_t *atc, pmegrids_t *grids,
                    gmx_bool bCalcSplines, gmx_bool bSpread,
                    real *fftgrid, gmx_bool bDoSplines, int grid_index)
{
    int nthread, thread;

    nthread = pme->nthread;
    assert(nthread > 0);

    if (bCalcSplines)
    {
#pragma omp parallel for num_threads(nthread) schedule(static)
        for (thread = 0; thread < nthread; thread++)
        {
            int start, end;

            start = atc->n* thread   /nthread;
            end   = atc->n*(thread+1)/nthread;

            /* Compute fftgrid index for all atoms,
             * with help of some extra variables.
             */
            calc_interpolation_idx(pme, atc, start, grid_index, end, thread);
        }
    }

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        splinedata_t *spline;
        pmegrid_t *grid = NULL;

        /* make local bsplines  */
        if (grids == NULL || !pme->bUseThreads)
        {
            spline = &atc->spline[0];

            spline->n = atc->n;

            if (bSpread)
            {
                grid = &grids->grid;    // 取地址使得grid成指针
            }
        }

        if (bCalcSplines)
        {
            make_bsplines(spline->theta, spline->dtheta, pme->pme_order,
                          atc->fractx, spline->n, spline->ind, atc->coefficient, bDoSplines);
        }

        if (bSpread)
        {
            /* put local atoms on grid. */
            spread_coefficients_bsplines_thread(grid, atc, spline, pme->spline_work);
        }
    }
}




/*
* [MODIFY] theta --- spline->theta
* dtheta --- spline->dtheta
* xptr --- fptr
*/ 

/* Macro to force loop unrolling by fixing order.
 * This gives a significant performance gain.
 */
#define CALC_SPLINE_NEW(order)                     \
    {                                              \
        for (int j = 0; (j < DIM); j++)            \
        {                                          \
            real dr, div;                          \
            real data[PME_ORDER_MAX];              \
                                                   \
            dr  = fptr[j];                         \
                                               \
            /* dr is relative offset from lower cell limit */ \
            data[order-1] = 0;                     \
            data[1]       = dr;                          \
            data[0]       = 1 - dr;                      \
                                               \
            for (int k = 3; (k < order); k++)      \
            {                                      \
                div       = 1.0/(k - 1.0);               \
                data[k-1] = div*dr*data[k-2];      \
                for (int l = 1; (l < (k-1)); l++)  \
                {                                  \
                    data[k-l-1] = div*((dr+l)*data[k-l-2]+(k-l-dr)* \
                                       data[k-l-1]);                \
                }                                  \
                data[0] = div*(1-dr)*data[0];      \
            }                                      \
            /* differentiate */                    \
            spline->dtheta[j][i*order+0] = -data[0];       \
            for (int k = 1; (k < order); k++)      \
            {                                      \
                spline->dtheta[j][i*order+k] = data[k-1] - data[k]; \
            }                                      \
                                               \
            div           = 1.0/(order - 1);                 \
            data[order-1] = div*dr*data[order-2];  \
            for (int l = 1; (l < (order-1)); l++)  \
            {                                      \
                data[order-l-1] = div*((dr+l)*data[order-l-2]+    \
                                       (order-l-dr)*data[order-l-1]); \
            }                                      \
            data[0] = div*(1 - dr)*data[0];        \
                                               \
            for (int k = 0; k < order; k++)        \
            {                                      \
                spline->theta[j][i*order+k]  = data[k];    \
            }                                      \
        }                                          \
    }

void spread_on_grid_new(struct gmx_pme_t *pme,
                        pme_atomcomm_t *atc, pmegrids_t *grids,
                        gmx_bool bCalcSplines, gmx_bool bSpread,
                        real *fftgrid, gmx_bool bDoSplines, int grid_index)
{
    int             i;
    int             start, end;
    int             nx, ny, nz;
    real            rxx, ryx, ryy, rzx, rzy, rzz;
    int            *idxptr, tix, tiy, tiz;
    real           *xptr, *fptr, tx, ty, tz;
    pmegrid_t      *pmegrid;
    real           *grid;
    splinedata_t   *spline;
    int            pnx, pny, pnz, ndatatot;
    int            offx, offy, offz;
    real           coefficient;
    int            order, norder;
    int            i0, j0, k0;
    real          *thx, *thy, *thz;
    int            ithx, ithy, ithz;
    int            index_x, index_xy, index_xyz;
    real           valx, valxy;

    start = 0;
    end = atc->n;

    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;    

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    spline = &atc->spline[0];

    pmegrid = &grids->grid;

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx*pny*pnz;
    grid = grids->grid.grid;
    order = pmegrid->order;

    // 内存分配时已经使用calloc将grid清0了
    // for (i = 0; i < ndatatot; i++)
    // {
    //     grid[i] = 0;
    // }


    for(i = start; i < end; i++) 
    {
        coefficient = atc->coefficient[i];
        xptr   = atc->x[i];
        idxptr = atc->idx[i];
        fptr   = atc->fractx[i];
        norder = i*order;

        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );

        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        /* Because decomposition only occurs in x and y,
         * we never have a fraction correction in z.
         */
        fptr[XX] = tx - tix + pme->fshx[tix];
        fptr[YY] = ty - tiy + pme->fshy[tiy];
        fptr[ZZ] = tz - tiz;

        idxptr[XX] = pme->nnx[tix];
        idxptr[YY] = pme->nny[tiy];
        idxptr[ZZ] = pme->nnz[tiz];

        if(coefficient != 0.0)
        {
            CALC_SPLINE_NEW(4);

            i0   = idxptr[XX] - offx;
            j0   = idxptr[YY] - offy;
            k0   = idxptr[ZZ] - offz;

            thx = spline->theta[XX] + norder;
            thy = spline->theta[YY] + norder;
            thz = spline->theta[ZZ] + norder;

#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_SPREAD_SIMD4_ORDER4
#else
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme-simd4.h"
#else
                    DO_BSPLINE(4);
#endif

        }
    
    }

}



/* This has to be a macro to enable full compiler optimization with xlC (and probably others too) */
#define DO_BSPLINE_v1(order)                            \
    for (ithx = 0; (ithx < order); ithx++)                    \
    {                                                    \
        index_x = (i0+ithx)*sgdy*sgdz;                     \
        valx    = coefficient*thx[ithx];                          \
                                                     \
        for (ithy = 0; (ithy < order); ithy++)                \
        {                                                \
            valxy    = valx*thy[ithy];                   \
            index_xy = index_x+(j0+ithy)*sgdz;            \
                                                     \
            for (ithz = 0; (ithz < order); ithz++)            \
            {                                            \
                index_xyz        = index_xy+(k0+ithz);   \
                sub_grid[index_xyz] += valxy*thz[ithz];      \
            }                                            \
        }                                                \
    }

void spread_on_grid_new_v1(struct gmx_pme_t *pme,
                        pme_atomcomm_t *atc, pmegrids_t *grids,
                        gmx_bool bCalcSplines, gmx_bool bSpread,
                        real *fftgrid, gmx_bool bDoSplines, int grid_index)
{
    int             i, ii, jj, kk;
    int             start, end;
    int             nx, ny, nz;
    real            rxx, ryx, ryy, rzx, rzy, rzz;
    int            *idxptr, tix, tiy, tiz;
    real           *xptr, *fptr, tx, ty, tz;
    pmegrid_t      *pmegrid;
    real           *grid;
    splinedata_t   *spline;
    int            pnx, pny, pnz, ndatatot;
    int            offx, offy, offz;
    real           coefficient;
    int            order, norder;
    int            i0, j0, k0;
    real          *thx, *thy, *thz;
    int            ithx, ithy, ithz;
    int            index_x, index_xy, index_xyz;
    real           valx, valxy;
    real           *sub_grid;               // 存储子网格内容
    int            sub_grid_base[DIM];      // 存储子网格相对于原网格的索引
    int            sub_grid_sub_order[DIM]; // 计算子网格边长减去order的大小
    int            sgdx, sgdy, sgdz;        // 子网格三维边长 sub_grid_dimension[DIM]
    int            g_index, sg_index;       // 使用子网格更新网格时，计算二者的索引
    // 边界处理{SUBGRID > order时发生}：防止子网格边缘超出原网格边缘(防止上界即可),然后网格更新时使用这个值
    int            sub_boundx, sub_boundy, sub_boundz; 
    int            idxptr_offset[DIM];      // 原子索引相对于子网格的偏移索引 idxptr = atc->idx[i] ; idxptr[i][DIM](pme->nnx[tix]) - sub_grid_base[DIM]

    start = 0;
    end = atc->n;

    nx  = pme->nkx;
    ny  = pme->nky;
    nz  = pme->nkz;    

    rxx = pme->recipbox[XX][XX];
    ryx = pme->recipbox[YY][XX];
    ryy = pme->recipbox[YY][YY];
    rzx = pme->recipbox[ZZ][XX];
    rzy = pme->recipbox[ZZ][YY];
    rzz = pme->recipbox[ZZ][ZZ];

    spline = &atc->spline[0];

    pmegrid = &grids->grid;

    pnx = pmegrid->s[XX];
    pny = pmegrid->s[YY];
    pnz = pmegrid->s[ZZ];

    offx = pmegrid->offset[XX];
    offy = pmegrid->offset[YY];
    offz = pmegrid->offset[ZZ];

    ndatatot = pnx*pny*pnz;
    grid = grids->grid.grid;
    order = pmegrid->order;

    sgdx = SUB_GRID_X;
    sgdy = SUB_GRID_Y;
    sgdz = SUB_GRID_Z;

    sub_grid_sub_order[XX] = sgdx - order;
    sub_grid_sub_order[YY] = sgdy - order;
    sub_grid_sub_order[ZZ] = sgdz - order;
    sub_grid = (real *)calloc(sgdx*sgdy*sgdz, sizeof(real));

    // 内存分配时已经使用calloc将grid清0了
    // for (i = 0; i < ndatatot; i++)
    // {
    //     grid[i] = 0;
    // }

    for(i = start; i < end; i++) 
    {
        coefficient = atc->coefficient[i];
        if(0.0 == coefficient) continue;
        xptr   = atc->x[i];
        
        /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
        tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
        ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
        tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );        

        tix = (int)(tx);
        tiy = (int)(ty);
        tiz = (int)(tz);

        sub_grid_base[XX] = pme->nnx[tix];
        sub_grid_base[YY] = pme->nny[tiy];
        sub_grid_base[ZZ] = pme->nnz[tiz];

        sub_boundx = (sub_grid_base[XX] + sgdx <= pnx) ? sgdx : (pnx - sub_grid_base[XX]);
        sub_boundy = (sub_grid_base[YY] + sgdy <= pny) ? sgdy : (pny - sub_grid_base[YY]);
        sub_boundz = (sub_grid_base[ZZ] + sgdz <= pnz) ? sgdz : (pnz - sub_grid_base[ZZ]);

        while(1) {
            coefficient = atc->coefficient[i];
            if(0.0 == coefficient) {
                i++;        // 计算下一个原子
                continue;
            } 

            xptr   = atc->x[i];
            fptr   = atc->fractx[i];
            norder = i*order;

            /* Fractional coordinates along box vectors, add 2.0 to make 100% sure we are positive for triclinic boxes */
            tx = nx * ( xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx + 2.0 );
            ty = ny * (                  xptr[YY] * ryy + xptr[ZZ] * rzy + 2.0 );
            tz = nz * (                                   xptr[ZZ] * rzz + 2.0 );

            tix = (int)(tx);
            tiy = (int)(ty);
            tiz = (int)(tz);

            idxptr_offset[XX] = pme->nnx[tix] - sub_grid_base[XX];
            idxptr_offset[YY] = pme->nny[tiy] - sub_grid_base[YY];
            idxptr_offset[ZZ] = pme->nnz[tiz] - sub_grid_base[ZZ];

            if(idxptr_offset[XX] > sub_grid_sub_order[XX] ||
               idxptr_offset[YY] > sub_grid_sub_order[YY] ||
               idxptr_offset[ZZ] > sub_grid_sub_order[ZZ] ||
               idxptr_offset[XX] < 0 ||
               idxptr_offset[YY] < 0 ||
               idxptr_offset[ZZ] < 0)
            {
                i--;        // 外层循环i++，为了重新计算该原子所以有此操作
                break;
            }

            /* Because decomposition only occurs in x and y,
            * we never have a fraction correction in z.
            */
            fptr[XX] = tx - tix + pme->fshx[tix];
            fptr[YY] = ty - tiy + pme->fshy[tiy];
            fptr[ZZ] = tz - tiz;

            CALC_SPLINE_NEW(4);

            // idxptr[XX] - sub_grid_base[XX] 可以简化
            i0   = idxptr_offset[XX] - offx;
            j0   = idxptr_offset[YY] - offy;
            k0   = idxptr_offset[ZZ] - offz;
            
            //printf("i0 = %d, j0 = %d, k0 = %d\n", i0, j0, k0);

            thx = spline->theta[XX] + norder;
            thy = spline->theta[YY] + norder;
            thz = spline->theta[ZZ] + norder;
            
#ifdef PME_SIMD4_SPREAD_GATHER
#ifdef PME_SIMD4_UNALIGNED
#define PME_SPREAD_SIMD4_ORDER4_V1
#else
#define PME_SPREAD_SIMD4_ALIGNED
#define PME_ORDER 4
#endif
#include "pme-simd4.h"
#else
                    DO_BSPLINE_v1(4);
#endif
            i++;

            if(i == end - 1) break;
        }   // end of while


        // TODO: 用sub_grid更新grid
        for(ii = 0; ii < sub_boundx; ii++) {
            for(jj = 0; jj < sub_boundy; jj++) {
                for(kk = 0; kk < sub_boundz; kk++) {
                    g_index = (ii+sub_grid_base[XX])*pny*pnz + (jj+sub_grid_base[YY])*pnz + kk+sub_grid_base[ZZ];
                    sg_index = ii*sgdy*sgdz + jj*sgdz + kk;
                    grid[g_index] += sub_grid[sg_index];
                    // TODO： sub_grid清零
                    sub_grid[sg_index] = 0.0;
                }
            }
        }

    }   // end of for

}