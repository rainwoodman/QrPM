#include <math.h>
#include <stdlib.h>
#include <fastpm/libfastpm.h>
#include <fastpm/logging.h>
#include "pmpfft.h"

void
fastpm_apply_smoothing_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double sml)
{

#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        int d;
        int i;
        double *kernel[3];
        for(d = 0; d < 3; d ++) {
            kernel[d] = malloc(sizeof(double) * pm->Nmesh[d]);
            for(i = 0; i < pm->Nmesh[d]; i ++) {
                double kk = kiter.kk[d][i];
                kernel[d][i] = exp(- 0.5 * kk * sml * sml);
            }
        }
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            for(dir = 0; dir < 3; dir++)
                smth *= kernel[dir][kiter.iabs[dir]];

            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
        for(d = 0; d < 3; d ++) {
            free(kernel[d]);
        }
    }
}

void
fastpm_apply_lowpass_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double kth)
{
    double kth2 = kth * kth;
#pragma omp parallel 
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            double kk = 0;
            for(dir = 0; dir < 3; dir++) {
                kk += kiter.kk[dir][kiter.iabs[dir]];
            }
            if(kk < kth2) smth = 1;
            else smth = 0;
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
    }
}

static double sinc_unnormed(double x) {
    if(x < 1e-5 && x > -1e-5) {
        double x2 = x * x;
        return 1.0 - x2 / 6. + x2  * x2 / 120.;
    } else {
        return sin(x) / x;
    }
}


void 
fastpm_apply_decic_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to) 
{

#pragma omp parallel 
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        int d;
        int i;
        double *kernel[3];
        for(d = 0; d < 3; d ++) {
            kernel[d] = malloc(sizeof(double) * pm->Nmesh[d]);
            for(i = 0; i < pm->Nmesh[d]; i ++) {
                double w = kiter.k[d][i] * pm->BoxSize[d] / pm->Nmesh[d];
                double cic = sinc_unnormed(0.5 * w);
                /* Watchout: this does divide by sinc, not sinc 2, */
                kernel[d][i] = 1.0 / pow(cic, 2);
            }
        }
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            for(dir = 0; dir < 3; dir++) 
                smth *= kernel[dir][kiter.iabs[dir]];

            /* - i k[d] */
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
        for(d = 0; d < 3; d ++) {
            free(kernel[d]);
        }
    }
}

void
fastpm_apply_diff_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir)
{
    ptrdiff_t * Nmesh = pm_nmesh(pm);
#pragma omp parallel
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k_finite = kiter.k_finite[dir][kiter.iabs[dir]];
            /* i k[d] */
            if(
                kiter.iabs[0] == (Nmesh[0] - kiter.iabs[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - kiter.iabs[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - kiter.iabs[2]) % Nmesh[2]
            ) {
                /* We are at the nyquist and the diff operator shall be zero;
                 * otherwise the force is not real! */
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            } {
                FastPMFloat tmp[2];
                tmp[0] = - from[kiter.ind + 1] * (k_finite);
                tmp[1] =   from[kiter.ind + 0] * (k_finite);
                to[kiter.ind + 0] = tmp[0];
                to[kiter.ind + 1] = tmp[1];
            }
        }
    }
}

void
fastpm_apply_diff_transfer_green(PM * pm, FastPMFloat * from, FastPMFloat * to, int dir)
{
    ptrdiff_t * Nmesh = pm_nmesh(pm);
#pragma omp parallel
    {
        PMKIter kiter;
        for(pm_kiter_init(pm, &kiter);
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            double k = kiter.k[dir][kiter.iabs[dir]];
            /* i k[d] */
            if(
                kiter.iabs[0] == (Nmesh[0] - kiter.iabs[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - kiter.iabs[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - kiter.iabs[2]) % Nmesh[2]
            ) {
                /* We are at the nyquist and the diff operator shall be zero;
                 * otherwise the force is not real! */
                to[kiter.ind + 0] = 0;
                to[kiter.ind + 1] = 0;
            } {
                FastPMFloat tmp[2];
                tmp[0] = - from[kiter.ind + 1] * (k);
                tmp[1] =   from[kiter.ind + 0] * (k);
                to[kiter.ind + 0] = tmp[0];
                to[kiter.ind + 1] = tmp[1];
            }
        }
    }
}

void
fastpm_apply_laplace_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, int order)
{
    /* order = 0: 1 / kk gives an IC that agrees with linear theory better
     * at intermidiate scales (desi-cosmosim email archives).
     * than
     * order = 1: 1 / sinc(kk) used in original FastPM;
     * */

#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        float ** kklist [3] = {kiter.kk, kiter.kk_finite, kiter.kk_finite2};
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int d;
            double kk_finite = 0;
            for(d = 0; d < 3; d++) {
                kk_finite += kklist[order][d][kiter.iabs[d]];
            }
            ptrdiff_t ind = kiter.ind;
            /* 1 / k2 */
            if(LIKELY(kk_finite != 0)) {
                to[ind + 0] = from[ind + 0] * (1 / kk_finite);
                to[ind + 1] = from[ind + 1] * (1 / kk_finite);
            } else {
                to[ind + 0] = 0;
                to[ind + 1] = 0;
            }
        }
    }
}

void
fastpm_apply_any_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, fastpm_fkfunc func, void * data)
{
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double smth = 1.0;
            double kk = 0;
            for(dir = 0; dir < 3; dir++) {
                kk += kiter.kk[dir][kiter.iabs[dir]];
            }
            double k = sqrt(kk);
            smth = func(k, data);
            to[kiter.ind + 0] = from[kiter.ind + 0] * smth;
            to[kiter.ind + 1] = from[kiter.ind + 1] * smth;
        }
    }
}

void
fastpm_apply_multiply_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, double value)
{
    ptrdiff_t i;
#pragma omp parallel for
    for(i = 0; i < pm_allocsize(pm); i ++) {
        to[i] = from[i] * value;
    }
}

void
fastpm_apply_normalize_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to)
{
    double Norm = 0;
#pragma omp parallel reduction(+: Norm)
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            int dir;
            double kk = 0;
            for(dir = 0; dir < 3; dir++) {
                kk += kiter.kk[dir][kiter.iabs[dir]];
            }
            if(kk == 0) {
                Norm += from[kiter.ind + 0];
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &Norm, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    if(Norm == 0) {
        fastpm_raise(-1, "It makes no sense to normalize a field with a mean of zero.");
    }
    fastpm_apply_multiply_transfer(pm, from, to, 1 / Norm);
}

void
fastpm_apply_c2r_weight_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to)
{
    ptrdiff_t * Nmesh = pm_nmesh(pm);

#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            /* usually each mode in the complex array moves its dual mode too, */
            double weight = 2.0;
            /* but if my dual is myself, then I am used once */
            if(
                kiter.iabs[0] == (Nmesh[0] - kiter.iabs[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - kiter.iabs[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - kiter.iabs[2]) % Nmesh[2]
            ) {
                weight = 1.0;
/*                fastpm_ilog(0, "weight == 1 mode found at %td %td %td\n", kiter.iabs[0], kiter.iabs[1], kiter.iabs[2]); */
            }
            to[kiter.ind + 0] = weight * from[kiter.ind + 0];
            to[kiter.ind + 1] = weight * from[kiter.ind + 1];
        }
    }
}

void
fastpm_apply_modify_mode_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, ptrdiff_t * mode, double value)
{
    fastpm_apply_set_mode_transfer(pm, from, to, mode, value, 0);
}

/* method == 0 for override, method == 1 for add;
 * always use get_mode_transfer to check the result -- this doesn't guarentee setting the mode
 * to the value (because not all modes are free.
 * */
void
fastpm_apply_set_mode_transfer(PM * pm, FastPMFloat * from, FastPMFloat * to, ptrdiff_t * mode, double value, int method)
{
    ptrdiff_t * Nmesh = pm_nmesh(pm);
    if(
        mode[0] == (Nmesh[0] - mode[0]) % Nmesh[0] &&
        mode[1] == (Nmesh[1] - mode[1]) % Nmesh[1] &&
        mode[2] == (Nmesh[2] - mode[2]) % Nmesh[2]
    ) {
        if(mode[3] == 1) {
            /* These modes are purely real, thus we can not set the imag part to non-zero */
            method = 0;
            value = 0;
        }
    }
#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {
            to[kiter.ind + 0] = from[kiter.ind + 0];
            to[kiter.ind + 1] = from[kiter.ind + 1];

            if((
                kiter.iabs[0] == mode[0] &&
                kiter.iabs[1] == mode[1] &&
                kiter.iabs[2] == mode[2]
            )) {
                if(method == 0)
                    to[kiter.ind + mode[3]] = value;
                else
                    to[kiter.ind + mode[3]] += value;
            }
            if((
                kiter.iabs[0] == (Nmesh[0] - mode[0]) % Nmesh[0] &&
                kiter.iabs[1] == (Nmesh[1] - mode[1]) % Nmesh[1] &&
                kiter.iabs[2] == (Nmesh[2] - mode[2]) % Nmesh[2]
            )) {
                /* conjugate plane */
                if(method == 0)
                    to[kiter.ind + mode[3]] = value * ((mode[3] == 0)?1:-1);
                else
                    to[kiter.ind + mode[3]] += value * ((mode[3] == 0)?1:-1);
            }
        }
    }
}

double
fastpm_apply_get_mode_transfer(PM * pm, FastPMFloat * from, ptrdiff_t * mode)
{
    double result = 0.0;

#pragma omp parallel
    {
        PMKIter kiter;
        pm_kiter_init(pm, &kiter);
        for(;
            !pm_kiter_stop(&kiter);
            pm_kiter_next(&kiter)) {

            if((
                kiter.iabs[0] == mode[0] &&
                kiter.iabs[1] == mode[1] &&
                kiter.iabs[2] == mode[2]
            )) {
                result = from[kiter.ind + mode[3]];
            }
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, pm_comm(pm));
    return result;
}

