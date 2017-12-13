-- parameter file
------ Size of the simulation -------- 

-- For Testing
nc = 64
boxsize = 1024

-------- Time Sequence ----
-- linspace: Uniform time steps in a
-- time_step = linspace(0.025, 1.0, 39)
-- logspace: Uniform time steps in loga
-- time_step = linspace(0.01, 1.0, 10)
time_step = linspace(0.1, 1, 10)

output_redshifts= {9.0, 0.0}  -- redshifts of output
compute_potential = true

-- Cosmology --
omega_m = 0.307494
h       = 0.6774

-- Start with a power spectrum file
-- Initial power spectrum: k P(k) in Mpc/h units
-- Must be compatible with the Cosmology parameter
read_powerspectrum= "powerspec.txt"
random_seed= 100

-------- Approximation Method ---------------
force_mode = "fastpm"
-- force_mode = "cola"

pm_nc_factor = 2            -- Particle Mesh grid pm_nc_factor*nc per dimension in the beginning

np_alloc_factor= 8.0      -- Amount of memory allocated for particle

-------- Output ---------------

-- Dark matter particle outputs (all particles)
-- write_runpb_snapshot= "nbodykit/tpm"
write_snapshot= "lightcone/fastpm" 
-- 1d power spectrum (raw), without shotnoise correction
-- write_powerspectrum = "pm/powerspec"
write_lightcone = "lightcone/lightcone"
dh_factor = 1.0
fov = 0
glmatrix = { 
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 2100},
        {0, 0, 0, 0},
        }

tiles = {
    {0, 0, 3},
    {0, 0, 2},
    {0, 0, 1},
    {0, 0, 0},
    }
