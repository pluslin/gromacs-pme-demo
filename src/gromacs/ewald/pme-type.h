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
/*! \internal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

/* TODO This file is a temporary holding area for stuff local to the
 * PME code, before it acquires some more normal ewald/file.c and
 * ewald/file.h structure.  In future clean up, get rid of this file,
 * to build more normal. */

#ifndef GMX_EWALD_PME_INTERNAL_H
#define GMX_EWALD_PME_INTERNAL_H

#include "config.h"

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/math/vectypes.h"

//@{
//! Grid indices for A state for charge and Lennard-Jones C6
#define PME_GRID_QA    0
#define PME_GRID_C6A   2
//@}

//@{
/*! \brief Flags that indicate the number of PME grids in use */
#define DO_Q           2 /* Electrostatic grids have index q<2 */
#define DO_Q_AND_LJ    4 /* non-LB LJ grids have index 2 <= q < 4 */
#define DO_Q_AND_LJ_LB 9 /* With LB rules we need a total of 2+7 grids */
//@}

/*! \brief Pascal triangle coefficients scaled with (1/2)^6 for LJ-PME with LB-rules */
static const real lb_scale_factor[] = {
    1.0/64, 6.0/64, 15.0/64, 20.0/64,
    15.0/64, 6.0/64, 1.0/64
};

/*! \brief Pascal triangle coefficients used in solve_pme_lj_yzx, only need to do 4 calculations due to symmetry */
static const real lb_scale_factor_symm[] = { 2.0/64, 12.0/64, 30.0/64, 20.0/64 };

/*! \brief We only define a maximum to be able to use local arrays without allocation.
 * An order larger than 12 should never be needed, even for test cases.
 * If needed it can be changed here.
 */
#define PME_ORDER_MAX 12

/*! \brief As gmx_pme_init, but takes most settings, except the grid, from pme_src */
// int gmx_pme_reinit(struct gmx_pme_t **pmedata,
//                    t_commrec *        cr,
//                    struct gmx_pme_t * pme_src,
//                    const t_inputrec * ir,
//                    ivec               grid_size);

typedef real *splinevec[DIM];


typedef struct {
    int       n;
    int      *ind;
    splinevec theta;
    splinevec dtheta;
} splinedata_t;

typedef struct {
    int      n;
    rvec    *x;
    ivec    *idx;
    rvec    *fractx;
    real    *coefficient;
    splinedata_t   *spline;
}pme_atomcomm_t;

typedef struct {
    ivec  ci;     /* The spatial location of this grid         */
    ivec  n;      /* The used size of *grid, including order-1 */
    ivec  offset; /* The grid offset from the full node grid   */
    int   order;  /* PME spreading order                       */
    ivec  s;      /* The allocated size of *grid, s >= n       */
    real *grid;   /* The grid local thread, size n             */
} pmegrid_t;

typedef struct {
    pmegrid_t  grid;         /* The full node grid (non thread-local)            */
    int        nthread;      /* The number of threads operating on this grid     */
    ivec       nc;           /* The local spatial decomposition over the threads */
    pmegrid_t *grid_th;      /* Array of grids for each thread                   */
    real      *grid_all;     /* Allocated array for the grids in *grid_th        */
    int      **g2t;          /* The grid to thread index                         */
    ivec       nthread_comm; /* The number of threads to communicate with        */
} pmegrids_t;

typedef struct pme_spline_work{
    int dummy;
} spline_work;

typedef struct gmx_pme_t{
    int                      nkx, nky, nkz;
    matrix                   recipbox;
    real                     *fshx, *fshy;
    int                      *nnx, *nny, *nnz;
    int                      bUseThreads;
    int                      pme_order;
    struct pme_spline_work   *spline_work;
    pmegrids_t               pmegrid[DO_Q_AND_LJ_LB];
    int                      nthread;
    pme_atomcomm_t           atc[2];
}tgmx_pme_t;

#endif
