#include <string.h>
#include <math.h>
#include <mpi.h>

#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include <fastpm/transfer.h>

#include "pmpfft.h"
#include "pmghosts.h"
#include "pm2lpt.h"

void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, FastPMFuncK * growth_rate_func_k, FastPMStore * p, double shift[3], FastPMKernelType type, FastPMICKernelType ic_type)
{
    /* read out values at locations with an inverted shift */
    int potorder, gradorder, deconvolveorder;
    fastpm_kernel_type_get_orders(type, &potorder, &gradorder, &deconvolveorder);

    /* calculate dx1, dx2, for initial fluctuation delta_k.
     * shift: martin has shift = 0.5, 0.5, 0.5.
     * Use shift of 0, 0, 0 if in doublt. 
     *   */
    ptrdiff_t i;
    int d;

    /* It is important to (de-)shift the particles before creating ghosts 
     * Because we will read out from the (de-)shifted positions.
     * Otherwise the IC will have artifacts along the edges of domains. */
    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] -= shift[d];
        }
    }
    FastPMPainter painter[1];
    fastpm_painter_init(painter, pm, FASTPM_PAINTER_CIC, 0);

    PMGhostData * pgd;
    if (p->dv1) {   //could alternatively use if (growth_rate_func_k)
        pgd = pm_ghosts_create(pm, p, p->attributes | COLUMN_DX1 | COLUMN_DX2 | COLUMN_DV1, painter->support);
    } else {
        pgd = pm_ghosts_create(pm, p, p->attributes | COLUMN_DX1 | COLUMN_DX2, painter->support);
    }
    pm_ghosts_send(pgd, COLUMN_POS);

    FastPMFloat * source =  pm_alloc(pm);
    FastPMFloat * workspace = pm_alloc(pm);
    memset(source, 0, sizeof(source[0]) * pm->allocsize);
    memset(workspace, 0, sizeof(workspace[0]) * pm->allocsize);

    FastPMFloat * field[3];

    for(d = 0; d < 3; d++ ) {
        field[d] = pm_alloc(pm);
        memset(field[d], 0, sizeof(field[d][0]) * pm->allocsize);
    }
    FastPMFieldDescr DX1[] = { {COLUMN_DX1, 0}, {COLUMN_DX1, 1}, {COLUMN_DX1, 2}};
    FastPMFieldDescr DX2[] = { {COLUMN_DX2, 0}, {COLUMN_DX2, 1}, {COLUMN_DX2, 2}};
    FastPMFieldDescr DV1[] = { {COLUMN_DV1, 0}, {COLUMN_DV1, 1}, {COLUMN_DV1, 2}};  // this will do nothing throughout the function if !p->v1
    int D1[] = {1, 2, 0};
    int D2[] = {2, 0, 1};

    /* 1LPT */
    for(d = 0; d < 3; d++) {
        /* dx1 */
        /*Check the variable "ic_kernel_type" and decide who to call the functions*/
        switch(ic_type)
        {
            case FASTPM_IC_KERNEL_STANDARD:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d);
            break;
            case FASTPM_IC_KERNEL_GREEN_FUNCTION:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace,0);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d);
            break;
            default:
                fastpm_raise(-1, "Wrong kernel type\n");
        }

        pm_c2r(pm, workspace);

        fastpm_readout_local(painter, workspace, p, p->np, DX1[d]);
        fastpm_readout_local(painter, workspace, pgd->p, pgd->p->np, DX1[d]);

        /* dv1 */
        if (p->dv1) {
            switch(ic_type)
            {
            case FASTPM_IC_KERNEL_STANDARD:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d);
            break;
            case FASTPM_IC_KERNEL_GREEN_FUNCTION:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace,0);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d);
            break;
            default:
                fastpm_raise(-1, "Wrong kernel type\n");
            }
            
            fastpm_apply_any_transfer(pm, workspace, workspace, (fastpm_fkfunc) fastpm_funck_eval2, growth_rate_func_k);

            pm_c2r(pm, workspace);

            fastpm_readout_local(painter, workspace, p, p->np, DV1[d]);
            fastpm_readout_local(painter, workspace, pgd->p, pgd->p->np, DV1[d]);
        }
    }

    /* 2LPT */
    for(d = 0; d< 3; d++) {
        switch(ic_type)
        {
            case FASTPM_IC_KERNEL_STANDARD:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d);
            break;
            case FASTPM_IC_KERNEL_GREEN_FUNCTION:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace,0);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d);

            break;
            default:
                fastpm_raise(-1, "Wrong kernel type\n");
        }

        pm_c2r(pm, field[d]);
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] += field[d1][i] * field[d2][i];
        }
    }

    for(d = 0; d < 3; d++) {
        int d1 = D1[d];
        int d2 = D2[d];
        switch(ic_type)
        {
            case FASTPM_IC_KERNEL_STANDARD:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d1);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d2);
            break;
            case FASTPM_IC_KERNEL_GREEN_FUNCTION:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace,0);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d1);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d2);
            break;
            default:
                fastpm_raise(-1, "Wrong kernel type\n");
        }

        pm_c2r(pm, workspace);
#pragma omp parallel for
        for(i = 0; i < pm->IRegion.total; i ++) {
            source[i] -= workspace[i] * workspace[i];
        }
    } 
    pm_r2c(pm, source, workspace);
    pm_assign(pm, workspace, source);

    for(d = 0; d < 3; d++) {
        /* 
         * We absorb some the negative factor in za transfer to below;
         *
         * */

        switch(ic_type)
       {
        case FASTPM_IC_KERNEL_STANDARD:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace, potorder);
                fastpm_apply_diff_transfer(pm, workspace, workspace, d);
        break;
        case FASTPM_IC_KERNEL_GREEN_FUNCTION:
                fastpm_apply_laplace_transfer(pm, delta_k, workspace,0);
                fastpm_apply_diff_transfer_green(pm, workspace, workspace, d);
        break;
        default:
                fastpm_raise(-1, "Wrong kernel type\n");
        }
            

        pm_c2r(pm, workspace);

        /* this ensures x = x0 + dx1(t) + dx2(t) */
        fastpm_apply_multiply_transfer(pm, workspace, workspace, 3.0 / 7);

        fastpm_readout_local(painter, workspace, p, p->np, DX2[d]);
        fastpm_readout_local(painter, workspace, pgd->p, pgd->p->np, DX2[d]);
    }

    pm_ghosts_reduce(pgd, COLUMN_DX1, FastPMReduceAddFloat, NULL);
    pm_ghosts_reduce(pgd, COLUMN_DX2, FastPMReduceAddFloat, NULL);
    if (p->dv1) pm_ghosts_reduce(pgd, COLUMN_DV1, FastPMReduceAddFloat, NULL);

#ifdef PM_2LPT_DUMP
    fwrite(p->dx1, sizeof(p->dx1[0]), p->np, fopen("dx1.f4x3", "w"));
    fwrite(p->dx2, sizeof(p->dx2[0]), p->np, fopen("dx2.f4x3", "w"));
#endif

    for(i = 0; i < p->np; i ++) {
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += shift[d];
        }
    }

    for(d = 0; d < 3; d ++) {
        pm_free(pm, field[2-d]);
    }
    pm_free(pm, workspace);
    pm_free(pm, source);

    pm_ghosts_free(pgd);

}

// Interpolate position and velocity for snapshot at a=aout
void
pm_2lpt_evolve(double aout, FastPMStore * p, FastPMCosmology * c, int zaonly)
{
    int np = p->np;

    FastPMGrowthInfo gi;
    fastpm_growth_info_init(&gi, aout, c);

    double D1 = gi.D1;
    double D2 = gi.D2;
    double f1 = gi.f1;
    double f2 = gi.f2;
    double E = HubbleEa(aout, c);

    double dv1_prefac = D1 * aout * aout * E;
    double Dv1 = dv1_prefac * f1;
    double Dv2 = D2 * aout * aout * E * f2;

    fastpm_info("2LPT ICs set at z=%g: E=%g D1=%g, D2=%g, f1=%g, f2=%g\n", 1./aout-1, E, D1, D2, f1, f2);

    if(zaonly) {
        D2 = 0;
        Dv2 = 0;
    }

    int i;
#pragma omp parallel for
    for(i=0; i<np; i++) {
        int d;
        for(d = 0; d < 3; d ++) {
            p->x[i][d] += D1 * p->dx1[i][d] + D2 * p->dx2[i][d];

            if(p->v) {
                p->v[i][d] += p->dx2[i][d]*Dv2;
                if (p->dv1) {
                    p->v[i][d] += dv1_prefac * p->dv1[i][d];
                } else {
                    p->v[i][d] += Dv1 * p->dx1[i][d];
                }
            }
        }
    }
    p->meta.a_x = p->meta.a_v = aout;
}

