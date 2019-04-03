#include "gromacs/ewald/pme-spread.h"
#include "gromacs/ewald/pme-malloc-free.h"
#include "gromacs/ewald/pme-type.h"

void write_grid(real * grid);
void read_grid(real * &grid);
bool check_grid(real *ref_grid, real *com_grid);

// fftgrid并没有在spread_on_grid过程发挥作用
int main() {
    bool is_right;
    real *fftgrid = NULL;
    snew(fftgrid, GRIDP*GRIDP*GRIDP);
    tgmx_pme_t *pme = NULL;
    pme = pme_init();

    //printf("coordinate %f\n", pme->atc[0].x[1][2]);
    //printf("efficient %f\n", pme->atc[0].coefficient[222]);

    spread_on_grid_new(pme, &pme->atc[0], &pme->pmegrid[0], 1, 1, fftgrid, 1, 0);
    
// CHECK
    //write_grid(pme->pmegrid[0].grid.grid);
    read_grid(fftgrid);
    is_right = check_grid(fftgrid, pme->pmegrid[0].grid.grid);
    if(is_right) {
        printf("computed gird is right.\n");
    } else {
        printf("computed gird is false.\n");
    }

    gmx_pme_destroy(pme);
    if(pme == NULL)
        printf("free success\n");
    return 0;
}

bool check_grid(real *ref_grid, real *com_grid) {
    int pnx = GRIDP, pny = GRIDP, pnz = GRIDP;
    int tt, rr, ee, index;
    for(tt = 0; tt < pnx; tt++) {
        for(rr = 0; rr < pny; rr++) {
            for(ee = 0; ee < pnz; ee++) {
                index = tt*pny*pnz+rr*pnz+ee;
                if((ref_grid[index] - com_grid[index]) > 0.0001)
                    return false;
            }
        }
    }
    return true;
}

void read_grid(real * &grid) {
    int pnx = GRIDP, pny = GRIDP, pnz = GRIDP;

    char filename[100];
    FILE *fp;
    int tt, rr, ee;

    sprintf(filename, "output.grid");
    fp = fopen(filename, "r");
    for(tt = 0; tt < pnx; tt++) {
        for(rr = 0; rr < pny; rr++) {
            for(ee = 0; ee < pnz; ee++) {
                fscanf(fp, "%f ", &grid[tt*pny*pnz+rr*pnz+ee]);
            }
        }
    }
    fclose(fp);
}

void write_grid(real * grid) {
    int pnx = GRIDP, pny = GRIDP, pnz = GRIDP;

    char filename[100];
    FILE *fp;
    int tt, rr, ee;

    sprintf(filename, "example.grid");
    fp = fopen(filename, "w");
    for(tt = 0; tt < pnx; tt++) {
        for(rr = 0; rr < pny; rr++) {
            for(ee = 0; ee < pnz; ee++) {
                fprintf(fp, "%f ", grid[tt*pny*pnz+rr*pnz+ee]);
            }
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

    // printf("Begin\n");
    // int i,j;
    // //printf("%f\n", pme->atc[0].spline->theta[0][111]);
    // //printf("%f\n",pme->pmegrid[0].grid.grid[0]);
    // //printf("%f\n", pme->atc[0].x[1][2]);
    // printf("\nEnd\n");
