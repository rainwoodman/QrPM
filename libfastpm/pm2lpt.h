
void 
pm_2lpt_init(PM * pm, FastPMStore * p, int Ngrid, double BoxSize, MPI_Comm comm);

void 
pm_2lpt_solve(PM * pm, FastPMFloat * delta_k, FastPMFuncK * growth_rate_func_k, FastPMStore * p, double shift[3], FastPMKernelType type, FastPMICKernelType ic_type);

void 
pm_2lpt_evolve(double aout, FastPMStore * p, FastPMCosmology * c, int zaonly);
