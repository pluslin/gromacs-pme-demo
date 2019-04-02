#include<stdio.h>
#include<stdlib.h>

#define N 121687
#define GRID 72
#define GRIDP 75
#define DIM 3
#define DO_Q_AND_LJ_LB 9
#define ORDER 4

typedef float real;
typedef real matrix[DIM][DIM];
typedef real rvec[DIM];
typedef int ivec[DIM];
typedef real *splinevec[DIM];

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

void save_free(const char *name, const char *file, int line, void *ptr)
{
    if (ptr != NULL)
    {
        free(ptr);
    }
}

#define sfree(ptr) save_free(#ptr, __FILE__, __LINE__, (ptr))

void *save_calloc(const char *name, const char *file, int line, size_t nelem, size_t elsize) {
    void *p;

    p = NULL;
    if((nelem == 0) || (elsize == 0))
    {
        p = NULL;
    }
    else
    {
        if((p = calloc((size_t)nelem, (size_t)elsize)) == NULL)
        {
            printf("Not enough memory. Failed to calloc for %s\n\
            (called from file %s, line %d)", name, file, line);
        }
    }
    return p;
}

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

typedef struct {
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
}gmx_pme_t;

void init_recipbox(matrix rbox)
{
    rbox[0][0] = 0.093561;
    rbox[0][1] = 0.0;
    rbox[0][2] = 0.0;
    rbox[1][0] = 0.0;
    rbox[1][1] = 0.093561;
    rbox[1][2] = 0.0;
    rbox[2][0] = 0.0;
    rbox[2][1] = 0.0;
    rbox[2][2] = 0.093561;
}

// nnx/nny/nnz begin from zero
// 这里需要使用到指针引用，否则传进来的指针所开辟的空间传不出去
// 如果指针在函数内开辟空间，需使用指针引用
void init_nnxyz(int * &nnx, int * &nny, int * &nnz) {
    snew(nnx, 5*GRID);
    snew(nny, 5*GRID);
    snew(nnz, 5*GRID);
    int i;
    for(i = 0; i < 5*GRID; i++) {
        nnx[i] = i % GRID;
        nny[i] = i % GRID;
        nnz[i] = i % GRID;
    }
}

void init_pmegrid(pmegrid_t &grid) {
    grid.ci[0] = 0;
    grid.ci[1] = 0;
    grid.ci[2] = 0;
    grid.n[0] = GRIDP;
    grid.n[1] = GRIDP;
    grid.n[2] = GRIDP;
    grid.offset[0] = 0;
    grid.offset[1] = 0;
    grid.offset[2] = 0;
    grid.order = 4;
    grid.s[0] = GRIDP;
    grid.s[1] = GRIDP;
    grid.s[2] = GRIDP;
    snew(grid.grid, GRIDP*GRIDP*GRIDP);
}

void init_g2t(int ** &g2t) {
    int i, j;
    // g2t = (int **)calloc(sizeof(int*)*DIM);
    // for(i = 0; i < DIM; i++) {
    //     g2t[i] = (int *)calloc(sizeof(int)*GRID);
    // }
    snew(g2t, DIM);
    for(i = 0; i < DIM; i++) {
        snew(g2t[i], GRID);
    }
}

void init_pmegrids(pmegrids_t &grids) {
    init_pmegrid(grids.grid);
    grids.nthread = 1;
    grids.nc[0] = 1;
    grids.nc[1] = 1;
    grids.nc[2] = 1;
    grids.grid_th = NULL;
    grids.grid_all = NULL;
    init_g2t(grids.g2t);        // 将一个二维指针初始化为0 -- [DIM][GRID]
    grids.nthread_comm[0] = 1;
    grids.nthread_comm[1] = 1;
    grids.nthread_comm[2] = 1;
}

void init_coordinate(rvec * &x) {
    snew(x, N);
    char fname[100];
    sprintf(fname, "input.atom");
    FILE *fp = fopen(fname, "r");
    int i;
    for(i = 0; i < N; i++) {
        fscanf(fp, "%f %f %f", &x[i][0], &x[i][1], &x[i][2]);
    }
    fclose(fp);
}

void init_idx(ivec * &idx) {
    snew(idx, N);
}

void init_fractx(rvec * &fractx) {
    snew(fractx, N);
}

void init_coefficient(real * &coefficient) {
    snew(coefficient, N);
    char fname[100];
    sprintf(fname, "input.coefficient");
    FILE *fp = fopen(fname, "r");
    int i;
    for(i = 0; i < N; i++) {
        fscanf(fp, "%f", &coefficient[i]);
    }
    fclose(fp);
}

void init_ind(int * &ind) {
    snew(ind, N);
    int i;
    for(i = 0; i < N; i++) {
        ind[i] = i;
    }
}

void init_splinedata(splinedata_t * &spline) {
    snew(spline, 1);
    spline->n = N;
    snew(spline->theta[0], N*ORDER);
    snew(spline->theta[1], N*ORDER);
    snew(spline->theta[2], N*ORDER);
    snew(spline->dtheta[0], N*ORDER);
    snew(spline->dtheta[1], N*ORDER);
    snew(spline->dtheta[2], N*ORDER);
    init_ind(spline->ind);
}

void init_atomcomm(pme_atomcomm_t &atc) {
    atc.n = N;
    init_coordinate(atc.x);
    init_idx(atc.idx);
    init_fractx(atc.fractx);
    init_coefficient(atc.coefficient);
    init_splinedata(atc.spline);
}

gmx_pme_t *pme_init()
{
    int xx;
    gmx_pme_t *pme = NULL;
    snew(pme, 1);
    pme->nkx = GRID;
    pme->nky = GRID;
    pme->nkz = GRID;
    init_recipbox(pme->recipbox);
    snew(pme->fshx, 5*GRID);
    snew(pme->fshy, 5*GRID);
    init_nnxyz(pme->nnx, pme->nny, pme->nnz);
    pme->bUseThreads = 0;
    pme->pme_order = 4;
    pme->spline_work = NULL;        // 此为结构体指针，可以这么进行处理
    init_pmegrids(pme->pmegrid[0]);   // pmegrid作为结构体，其内部空间在snew(pme,1)已经开辟了
    pme->nthread = 1;
    init_atomcomm(pme->atc[0]);
    return pme;
}

void pmegrids_destroy(pmegrids_t &grids) {
    int i;
    sfree(grids.grid.grid);
    for(i = 0; i < DIM; i++) {
        sfree(grids.g2t[i]);
    }
    sfree(grids.g2t);
}

void splinedata_destroy(splinedata_t * &spline) {
    sfree(spline->theta[0]);
    sfree(spline->theta[1]);
    sfree(spline->theta[2]);
    sfree(spline->dtheta[0]);
    sfree(spline->dtheta[1]);
    sfree(spline->dtheta[2]);
    sfree(spline->ind);
    sfree(spline);
    spline=NULL;
}

void atomcomm_destroy(pme_atomcomm_t &atc) {
    sfree(atc.x);
    sfree(atc.idx);
    sfree(atc.fractx);
    sfree(atc.coefficient);
    splinedata_destroy(atc.spline);
}

int gmx_pme_destroy(gmx_pme_t * &pmedata) {
    sfree(pmedata->fshx);
    sfree(pmedata->fshy);
    sfree(pmedata->nnx);
    sfree(pmedata->nny);
    sfree(pmedata->nnz);
    pmegrids_destroy(pmedata->pmegrid[0]);
    atomcomm_destroy(pmedata->atc[0]);
    sfree(pmedata);
    pmedata=NULL;

    return 0;
}

int main() {
    gmx_pme_t *pme = NULL;
    pme = pme_init();

    printf("Begin\n");
    int i,j;
    //printf("%f\n", pme->atc[0].spline->theta[0][111]);
    //printf("%f\n",pme->pmegrid[0].grid.grid[0]);
    //printf("%f\n", pme->atc[0].x[1][2]);
    printf("\nEnd\n");

    gmx_pme_destroy(pme);
    if(pme == NULL)
        printf("free success\n");
    return 0;
}

/*  读文件
    char fname[100];
    sprintf(fname, "input.atom");
    FILE *fp = fopen(fname, "r");
    int i;
    float x[N][3];
    for(i = 0; i < N; i++) {
        fscanf(fp, "%f %f %f", &x[i][0], &x[i][1], &x[i][2]);
    }
    fclose(fp);
*/