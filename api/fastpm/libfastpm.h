#ifndef __FASTPM_H__
#define __FASTPM_H__

#include <stddef.h>
#include <stdint.h>
#include <mpi.h>

#ifndef FASTPM_BEGIN_DECLS
#ifdef __cplusplus
#define FASTPM_BEGIN_DECLS extern "C" {
#define FASTPM_END_DECLS }
#else
#define FASTPM_BEGIN_DECLS
#define FASTPM_END_DECLS
#endif
#endif

FASTPM_BEGIN_DECLS

typedef struct PM PM;
typedef struct FastPMStore FastPMStore;
typedef struct FastPMPainter FastPMPainter;
typedef struct FastPMTransition FastPMTransition;
typedef struct FastPMCosmology FastPMCosmology;
typedef struct FastPMHorizon FastPMHorizon;

#ifndef FASTPM_FFT_PRECISION
#define FASTPM_FFT_PRECISION 32
#endif

#if FASTPM_FFT_PRECISION == 64
    typedef double FastPMFloat;
#elif FASTPM_FFT_PRECISION == 32
    typedef float FastPMFloat;
#else
    #error FASTPM_FFT_PRECISION must be 32 or 64
#endif

typedef enum { FASTPM_FORCE_FASTPM = 0,
               FASTPM_FORCE_PM,
               FASTPM_FORCE_COLA,
               FASTPM_FORCE_2LPT,
               FASTPM_FORCE_ZA,
} FastPMForceType;

typedef enum { FASTPM_KERNEL_3_4, FASTPM_KERNEL_3_2, FASTPM_KERNEL_5_4,
               FASTPM_KERNEL_1_4,
               FASTPM_KERNEL_GADGET,
               FASTPM_KERNEL_EASTWOOD,
               FASTPM_KERNEL_NAIVE,
            } FastPMKernelType;
typedef enum { FASTPM_SOFTENING_NONE,
               FASTPM_SOFTENING_GAUSSIAN, FASTPM_SOFTENING_GADGET_LONG_RANGE,
               FASTPM_SOFTENING_TWO_THIRD, FASTPM_SOFTENING_GAUSSIAN36 } FastPMSofteningType;
typedef enum { FASTPM_IC_KERNEL_GREEN_FUNCTION, FASTPM_IC_KERNEL_STANDARD } FastPMICKernelType;

void libfastpm_init();
void libfastpm_cleanup();
void libfastpm_set_memory_bound(size_t size);

extern const char * LIBFASTPM_VERSION;

FASTPM_END_DECLS

typedef double (*fastpm_fkfunc)(double k, void * data);
typedef void   (*fastpm_posfunc)(FastPMStore * p, ptrdiff_t index, double pos[3]);
typedef double (*fastpm_kernelfunc)(double x, double hsupport);
#define fastpm_pkfunc fastpm_fkfunc

#include "events.h"
#include "memory.h"
#include "pmapi.h"
#include "store.h"
#include "painter.h"

#include "FDinterp.h"
#include "cosmology.h"
#include "horizon.h"
#include "transfer.h"
#include "initialcondition.h"
#include "pngaussian.h"
#include "powerspectrum.h"
#include "pgdcorrection.h"
#include "solver.h"
#include "gravity.h"
#include "timemachine.h"
#include "lightcone.h"
#include "utils.h"

#include "thermalvelocity.h"

/* Following functions are internal API */

FASTPM_BEGIN_DECLS

FastPMMemory * _libfastpm_get_gmem();

FASTPM_END_DECLS

#endif
