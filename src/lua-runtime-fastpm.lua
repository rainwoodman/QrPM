----------------------------------------------------
-- This is the main LUA runtime library of FastPM.
--
-- Author: Yu Feng <rainwoodman@gmail.com> 2016
----------------------------------------------------

local config = require("lua-runtime-config")

local _NAME = ... or 'main'

local fastpm = {}

local schema = config.Schema()
schema.declare{name='nc',                type='int', required=true, help="Number of Particles Per side"}
schema.declare{name='boxsize',           type='number', required=true, help="Size of box in Mpc/h"}
schema.declare{name='time_step',         type='array:number',  required=true, help="Scaling factor of steps, can be linspace(start, end, Nsteps)." }
schema.declare{name='output_redshifts',  type='array:number',  required=false, help="Redshifts for outputs" }
schema.declare{name='aout',              type='array:number',  required=false, help='a of redshifts'}

-- set aout from output_redshifts
function schema.aout.action(aout)
    schema.output_redshifts.required = false
end

-- set aout from output_redshifts
function schema.output_redshifts.action(output_redshifts)
    local aout = {}
    if output_redshifts ~= nil then
        for i, z in pairs(output_redshifts) do
            aout[i] = 1.0 / (z + 1.)
        end
        schema.aout.default = aout
    end
end

-- Note: the order of variable declaration is important when applying scheme.variable.action later in this file.
schema.declare{name='omega_m',           type='number', required=false, help='This is depreciated. Please use Omega_m (uppercase O).'}
schema.declare{name='Omega_m',           type='number', required=false, help='Total matter (cdm + baryon + ncdm) density parameter at z=0'}
schema.declare{name='T_cmb',             type='number', required=false, default=0, help="CMB temperature in K, 0 to turn off radiation."}
schema.declare{name='h',                 type='number', required=true, default=0.7, help="Dimensionless Hubble parameter"}
schema.declare{name='Omega_k',           type='number', required=false, default=0, help="Curvature density parameter. (Omega_k > 0 is an open universe.)"}
schema.declare{name='w0',                type='number', required=false, default=-1, help="Dark energy equation of state 0th order parameter: w(a) = w0 + (1-a) wa."}
schema.declare{name='wa',                type='number', required=false, default=0, help="Dark energy equation of state 1st order parameter: w(a) = w0 + (1-a) wa."}
schema.declare{name='N_eff',             type='number', required=false, default=3.046}
schema.declare{name='N_nu',              type='number', required=false, default=0, help="Total number of neutrinos, massive and massless."}
schema.declare{name='m_ncdm',            type='array:number', required=false, default={}, help="Mass of ncdm particles in eV. Enter in descending order."}
schema.declare{name='pm_nc_factor',      type='array:number',  required=true, help="A list of {a, PM resolution}, "}
schema.declare{name='lpt_nc_factor',     type='number', required=false, default=1, help="PM resolution use in lpt and linear density field."}
schema.declare{name='np_alloc_factor',   type='number', required=true, help="Over allocation factor for load imbalance" }
schema.declare{name='compute_potential', type='boolean', required=false, default=false, help="Calculate the gravitional potential."}
schema.declare{name='n_shell',           type='number', required=false, default=10, help="Number of shells of FD distribution for ncdm splitting. Set n_shell=0 for no ncdm particles."}
schema.declare{name='lvk',               type='boolean', required=false, default=true, help="Use the low velocity kernel when splitting FD for ncdm."}
schema.declare{name='n_side',            type='number', required=false, default=3, help="This is N_fib for fibonacci sphere splitting, or number of sides in HEALPix splitting."}
schema.declare{name='every_ncdm',        type='number', required=false, default=4, help="Subsample ncdm from cdm every..."}
schema.declare{name='ncdm_sphere_scheme',type='enum', required=false, default="fibonacci", help="Split sphere with 'fibonacci' or 'healpix'?"}
schema.ncdm_sphere_scheme.choices = {
    healpix = 'FASTPM_NCDM_SPHERE_HEALPIX',
    fibonacci = 'FASTPM_NCDM_SPHERE_FIBONACCI',
}
schema.declare{name='ncdm_matterlike', type='boolean', required=false, default=true, help="Approximate ncdm as matter-like in the background? If true, Omega_ncdm~1/a^3."}
schema.declare{name='ncdm_freestreaming', type='boolean', required=false, default=true, help="Treat ncdm as free-streaming? If true, source terms ~Omega_c; if false, ~Omega_m."}

-- lra parameters
schema.declare{name='ncdm_linearresponse', type='boolean', required=false, default=false, help="Enable the linear response module for the neutrino contribution. Adds neutrinos to the code without needing extra particles."}
schema.declare{name='ncdm_transfer_redshift', type='number', required=false, help="Redshift of the neutrino transfer function for lra. This should be at the starting redshift of the simulation."}
schema.declare{name='ncdm_transfer_nu_file', type='file', required=false, help="File to store initial transfer function for linear response neutrinos. Needs to contain T_nu / T_cdm + b."}

schema.declare{name='growth_mode', type='enum', default='ODE', help="Evaluate growth factors using a Lambda+CDM-only approximation or with the full ODE. " ..
                                                                     "The full ODE is required for accurate results for runs with radiation or varying DE in the background, " ..
                                                                     "and can also be used for Lambda+CDM-only backgrounds. " ..
                                                                     "The LCDM approximation is included for backward compatibility."}
schema.growth_mode.choices = {
    LCDM = 'FASTPM_GROWTH_MODE_LCDM',
    ODE = 'FASTPM_GROWTH_MODE_ODE',
}

-- enforce Omega_m
function schema.omega_m.action (value)
    if value ~= nil then
        error("omega_m is depreciated, please use Omega_m (uppercase O) instead.")
    end
end

-- check for bad input
function schema.time_step.action (time_step)
    function schema.T_cmb.action (T_cmb)
        function schema.growth_mode.action (growth_mode)
            if T_cmb ~= 0 and growth_mode ~= 'ODE' then
                error("For a run with radiation (T_cmb > 0) use growth_mode='ODE' for accurate results.")
            end
            -- enforce a<0.00625 (hard coded) when using ODE growth mode
            if time_step[1] < 0.00625 and growth_mode == 'ODE' then
                error("Cannot start the simulation at a<0.00625 when growth_mode=='ODE'.")
            end
        end

        function schema.m_ncdm.action (m_ncdm)
            if #m_ncdm ~= 0 then
                for i=2, #m_ncdm do
                    if m_ncdm[i] > m_ncdm[1] then
                        error("Please input the heaviest ncdm particle first.")
                    end
                end
                function schema.n_shell.action (n_shell)
                    function schema.ncdm_freestreaming.action (ncdm_freestreaming)
                        function schema.ncdm_linearresponse.action (ncdm_linearresponse)
                            function schema.ncdm_transfer_redshift.action (ncdm_transfer_redshift)
                                if ncdm_freestreaming and n_shell ~= 0 then
                                     error("For free-streaming ncdm use n_shell = 0 to turn off ncdm particles.")
                                end
                                if ncdm_linearresponse and not ncdm_freestreaming then
                                    error("For linear-response neutrino approach must use free-streaming.")
                                end
                                if ncdm_linearresponse and time_step[1] ~= 1./(1+ncdm_transfer_redshift) then
                                     error("Redshift of the neutrino transfer function for lra should be at the starting redshift of the simulation.")
                                end
                            end
                        end
                    end
                end
                function schema.ncdm_matterlike.action (ncdm_matterlike)
                    if not ncdm_matterlike and T_cmb == 0 then
                         error("For a run with exact Omega_ncdm, T_cmb > 0 is required.")
                    end
                end
            end
        end
    end
end

-- Force calculation --
schema.declare{name='painter_type',        type='enum', default='cic', help="Type of painter."}
schema.declare{name='painter_support',     type='int', default=2, help="Support (size) of the painting kernel"}
schema.painter_type.choices = {
    cic = 'FASTPM_PAINTER_CIC',
    linear = 'FASTPM_PAINTER_LINEAR',
    lanczos = 'FASTPM_PAINTER_LANCZOS',
}
function schema.painter_type.action(painter_type)
    if painter_type ~= 'cic' then
        schema.painter_support.required = true
    end
end
schema.declare{name='force_mode',        type='enum', default='fastpm'}
schema.force_mode.choices = {
    cola = 'FASTPM_FORCE_COLA',
    zola = 'FASTPM_FORCE_FASTPM',
    fastpm = 'FASTPM_FORCE_FASTPM',
    pm   = 'FASTPM_FORCE_PM',
}
schema.declare{name='enforce_broadband_kmax',  type='int', default=4}

-- Primordial Non-Gaussianity --
schema.declare{name='f_nl_type', type='enum', default='none'}
schema.f_nl_type.choices = {
    ['local'] = 'FASTPM_FNL_LOCAL',
    ['none']  = 'FASTPM_FNL_NONE',
}
schema.declare{name='f_nl', type='number'}
schema.declare{name='kmax_primordial_over_knyquist', type='number', default=0.666} --0.25->0.667 This is enough to avoid the aliasing associated to phi^2. Using lower values affect the bispectrum at small scales and the clustering of low mass haloes.
schema.declare{name='scalar_amp', type='number'}
schema.declare{name='scalar_pivot', type='number'}
schema.declare{name='scalar_spectral_index', type='number'}
function schema.f_nl_type.action (f_nl_type)
    if f_nl_type ~= 'none' then
        schema.f_nl.required = true
        schema.scalar_amp.required = true
        schema.scalar_pivot.required = true
        schema.scalar_spectral_index.required = true
    end
end

-- Initial condition --
schema.declare{name='read_lineark',        type='string', help='lineark for cdm'}
schema.declare{name='read_powerspectrum', type='file', help='file to read the linear power spectrum for cdm.'}
schema.declare{name='read_linear_growth_rate', type ='file', help='file to read the linear growth rate (f_1) of cdm. If left empty, will use internal f_1.'}
schema.declare{name='linear_density_redshift', type='number', default=0, help='redshift of the input linear cdm density field. '}

schema.declare{name='read_lineark_ncdm', type='string', help='file to read the lineark of ncdm.'}
schema.declare{name='read_powerspectrum_ncdm', type='file', help='file to read the linear power spectrum of ncdm.'}
schema.declare{name='read_linear_growth_rate_ncdm', type ='file', help='file to read the linear growth rate (f_1) of ncdm. If left empty, will use internal f_1.'}
schema.declare{name='linear_density_redshift_ncdm', type='number', default=0, help='redshift of the input linear ncdm density field.'}

schema.declare{name='read_grafic',        type='string'}
schema.declare{name='read_runpbic',       type='string'}
schema.declare{name='read_whitenoisek',         type='string'}

schema.declare{name='sigma8',             type='number', default=0, help='normalize linear power spectrumt to sigma8(z); this shall be sigma8 at linear_density_redshift, not z=0.'}
schema.declare{name='random_seed',         type='int'}
schema.declare{name='shift',             type='boolean', default=false}
schema.declare{name='inverted_ic',             type='boolean', default=false}
schema.declare{name='remove_cosmic_variance',  type='boolean', default=false}

schema.declare{name='ic_kernel_type',       type='enum', default='green', help='Type of kernel used during initial condition generation. The standard ones gives a suppression of the power spectrum at small scales. The green function agrees with linear theory up to knyquist.'}
schema.ic_kernel_type.choices = {
        green = 'FASTPM_IC_KERNEL_GREEN_FUNCTION',
        standard = 'FASTPM_IC_KERNEL_STANDARD',
    }

function schema.read_grafic.action (read_grafic)
    if read_grafic ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = true
    end
end

function schema.read_lineark.action (read_lineark)
    if read_lineark ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = false
    end
end
function schema.read_runpbic.action (read_runpbic)
    if read_runpbic ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = false
    end
end
function schema.read_whitenoisek.action (read_whitenoisek)
    if read_whitenoisek ~= nil then
        schema.random_seed.required = false
        schema.read_powerspectrum.required = true
    end
end

schema.declare{name='write_linearr',         type='string', help='Write linear field ICs in real space'}
schema.declare{name='write_lineark',         type='string', help='Write linear field ICs in Fourier space'}
schema.declare{name='write_whitenoisek',         type='string'}
schema.declare{name='write_runpbic',       type='string'}
schema.declare{name='write_powerspectrum', type='string'}
schema.declare{name='write_snapshot',      type='string'}
schema.declare{name='write_nonlineark',      type='string'}
schema.declare{name='write_runpb_snapshot', type='string'}
schema.declare{name='particle_fraction',    type='number', default=1.0, help='Fraction of particles to save in the snapshot (sub-sampling)'}
schema.declare{name='sort_snapshot',    type='boolean', default=true, help='sort snapshots by ID; very large communication is incurred during snapshots.'}

schema.declare{name='write_fof',      type='string', help='Path to save the fof catalog, will be in the FOF-0.200 dataset. (or other linking length).'}
schema.declare{name='fof_linkinglength',      type='number', default=0.2, help='linking length of FOF; in units of particle mean separation.'}
schema.declare{name='fof_nmin',      type='number', default=20, help='threshold for making into the FOF catalog.'}
schema.declare{name='fof_kdtree_thresh',      type='number', default=8, help='threshold for spliting a kdtree node. KDTree is used in fof. smaller uses more memory but fof runs faster.'}

schema.declare{name='write_rfof',      type='string', help='path to save the RFOF dataset; parameter defualts are fitted for Illustris sep=0.3 Mpc/h, 40 steps.'}
schema.declare{name='rfof_kdtree_thresh',      type='number', default=8, help='threshold for splitting a kdtree node.'}
schema.declare{name='rfof_linkinglength',      type='number', default=0.2, help='base linking length (0.2 in eq 2.4)'}
schema.declare{name='rfof_nmin',      type='number', default=8, help='threshold for making into the R-FOF catalog.'}
schema.declare{name='rfof_l1',      type='number', default=0.25, help=''}
schema.declare{name='rfof_l6',      type='number', default=0.24, help=''}
schema.declare{name='rfof_a1',      type='number', default=0.012, help=''}
schema.declare{name='rfof_a2',      type='number', default=0.06, help=''}
schema.declare{name='rfof_b1',      type='number', default=7.02, help=''}
schema.declare{name='rfof_b2',      type='number', default=6.025, help=''}

schema.declare{name='lc_amin',
            type='number', help='min scale factor for truncation of lightcone.'}
schema.declare{name='lc_amax',
            type='number', help='max scale factor for truncation of lightcone.'}

schema.declare{name='lc_write_usmesh',         type='string', help='file name base for writing the particle lightcone'}
schema.declare{name='lc_usmesh_nslices',   type='number', default=256, help='Number of slices in AemitIndex attribute of lightcone data sets. also controls number of slices in healpix maps.'}

schema.declare{name='lc_usmesh_alloc_factor',     type='number', default=1.0,
                    help='allocation factor for the unstructured mesh, relative to alloc_factor. Use smaller number to ease OOM errors in lightcone.'}

schema.declare{name='lc_usmesh_fof_padding',     type='number', default=10.0,
                    help='padding in the line of sight direction for light cone fof. roughly the size of a halo.'}

schema.declare{name='lc_usmesh_ell_limit',    type='number', default=0,
               help='Subsample particle fraction depending on redshift, to match the ell. 0 to use particle_fraction for all redshifts'}

schema.declare{name='lc_usmesh_healpix_nside',     type='number', default=0, help='nside for healpix map. particle ID is slice_id * npix + ipix.'}

schema.declare{name='lc_usmesh_tiles',     type='array:number',
        default={
            {0, 0, 0},
        },
        help=[[tiling of the simulation box, in units of box edges.
              all tiles will be considered during lightcone construction.
              tiling occurs before the glmatrix.]]
        }

schema.declare{name='dh_factor',    type='number', default=1.0, help='Scale Hubble distance to amplify the lightcone effect'}
schema.declare{name='lc_fov',     type='number', default=0.0, help=' field of view of the sky in degrees. 0 for flat sky and 360 for full sky. The beam is along the z-direction after glmatrix.'}
schema.declare{name='lc_octants',     type='array:number', default={0, 1, 2, 3, 4, 5, 6, 7},
            help='list of octants to include when fov>=360 degrees.'}

schema.declare{name='lc_glmatrix',     type='array:number',
        default={
            {1, 0, 0, 0,},
            {0, 1, 0, 0,},
            {0, 0, 1, 0,},
            {0, 0, 0, 1,},
        },
        help=[[transformation matrix to move simulation coordinate (x, y, z, 1) to the observer coordinate with a left dot product.
               The observer is sitting at z=0 in the observer coordinate. The last column of the matrix is the translation in Mpc/h.
               use the translation and rotation methods provide in the intepreter to build the matrix. ]]}

schema.declare{name='za',                      type='boolean', default=false, help='use ZA initial condition not 2LPT'}

schema.declare{name='kernel_type',             type='enum', default="1_4", help='Force kernel; affects low mass halos 3_4 gives more low mass halos; 1_4 is consistent with fastpm-python.'}
schema.kernel_type.choices = {
    ['1_4'] = 'FASTPM_KERNEL_1_4',  -- consistent with fastpm-python
    ['3_4'] = 'FASTPM_KERNEL_3_4',  -- legacy fastpm
    ['gadget'] = 'FASTPM_KERNEL_GADGET', -- GADGET long range without exp smoothing.
    ['5_4'] = 'FASTPM_KERNEL_5_4',  -- very bad do not use
    ['eastwood'] = 'FASTPM_KERNEL_EASTWOOD',
    ['naive'] = 'FASTPM_KERNEL_NAIVE',
    ['3_2'] = 'FASTPM_KERNEL_3_2',
}
schema.declare{name='force_softening_type',             type='enum', default="none", help='Softening kernel (wipes out small scale force), very little effect)'}
schema.force_softening_type.choices = {
    none = 'FASTPM_SOFTENING_NONE',
    gaussian = 'FASTPM_SOFTENING_GAUSSIAN',
    gadget_long_range = 'FASTPM_SOFTENING_GADGET_LONG_RANGE',
    gaussian36 = 'FASTPM_SOFTENING_GAUSSIAN36',
    twothird = 'FASTPM_SOFTENING_TWO_THIRD',
}

schema.declare{name='constraints',      type='array:number',  help="A list of {x, y, z, peak-sigma}, giving the constraints in MPC/h units. "}
function schema.constraints.action (constraints)
    if constraints == nil then
        return
    end
    for i,v in pairs(constraints) do
        if #v ~= 4 then
            error("contraints must be a list of 4-vectors (x, y, z, peak-sigma)")
        end
    end
end
schema.declare{name='set_mode_method',  type='string', default="override", help="override or relative"}
schema.declare{name='set_mode',         type='array:number', help="A list of {kix, kiy, kiz, ri, value}, set the IC mode at integer k (ri for real and imag) to value"}

function schema.set_mode.action (set_mode)
    if set_mode == nil then
        return
    end
    for i,v in pairs(set_mode) do
        if #v ~= 5 then
            error("set_mode must be a list of 5-vectors (x, y, z, real_or_imag, value)")
        end
        if v[4] ~= 1 and v[4] ~=0 then
            error("the fourth component specifies real or imag part of the mode. must be 0 or 1")
        end
    end
end

schema.declare{name='pgdc',              type='boolean', default=false, help="if enable pgd correction"}
schema.declare{name='pgdc_alpha0',           type='number', default=0.8, help="alpha parameter in pgd correction"}
schema.declare{name='pgdc_A',           type='number', default=4.0, help="alpha parameter in pgd correction"}
schema.declare{name='pgdc_B',           type='number', default=8.0, help="alpha parameter in pgd correction"}
schema.declare{name='pgdc_kl',           type='number', default=2.0, help="filter large scale parameter in pgd correction"}
schema.declare{name='pgdc_ks',           type='number', default=10.0, help="filter small scale parameter in pgd correction"}

function fastpm.translation(dx, dy, dz)
-- generate a translation gl matrix that shifts the coordinates
    return {
        {1, 0, 0, dx},
        {0, 1, 0, dy},
        {0, 0, 1, dz},
        {0, 0, 0, 1 },
        }
end

function fastpm.outerproduct(a, b, c)
-- generate a list that is outer product of elements of a, b, c
    local r = {}
    for i=1, #a do
        for j=1, #b do
            for k=1, #c do
                r[#r + 1] = {a[i], b[j], c[k]}
            end
        end
    end
    return r
end

function fastpm.linspace(a, e, N, endpoint)
-- Similar to numpy.linspace always append the end
--
-- https://mail.scipy.org/pipermail/numpy-discussion/2016-February/075065.html
    if endpoint == nil then
        endpoint = true
    end

    local r = {}

    if endpoint then
        N1 = N - 1
    else
        N1 = N
    end

    for i=1,N do
        r[i] = 1.0 * (e - a) * (i - 1) / N1 + a
    end

    if endpoint then
        r[N] = e
    end
    return r
end

function fastpm.logspace(a, e, N)
-- a and end are in log10.
-- Returns N elements, including e.
    local r
    r = fastpm.linspace(a, e, N)
    for i, j in pairs(r) do
        r[i] = math.pow(10, j)
    end
    return r
end
function fastpm.blendspace(a, e, a1, a2)
    local r = {}
    local a = a
    local i = 1
    while a < e do
        r[i] = a
        dlna = math.pow(math.pow(1/a1, 2) + math.pow(a/a2, 2), -0.5)
        a = math.exp(math.log(a) + dlna)
        i = i + 1
    end
    r[i] = e
    return r
end

function fastpm.loglinspace(a, m, e, Nlog, Nlin)
-- Take Nlog log steps between a and m,
-- then Nlin lin steps between m and e.
-- a, m, ane e are in linear units.
    local r
    local s
    local t = {}
    local n = 0
    r = fastpm.logspace(math.log10(a), math.log10(m), Nlog+1)
    s = fastpm.linspace(m, e, Nlin+1)
    for i=1,#r do n=n+1; t[n]=r[i] end
    for i=2,#s do n=n+1; t[n]=s[i] end  -- ignore duplicate entry on boundary
    return t
end

function fastpm.test()
    ns = {
        __file__ = "standard.lua",
        boxsize = 384.0,
        cola_stdda = true,
        force_softening_type = "none",
        enforce_broadband_kmax = 4,
        f_nl = 0.100000000000000006,
        f_nl_type = "local",
        force_mode = "pm",
        h = 0.677400000000000002,
        inverted_ic = false,
        kernel_type = "3_4",
        nc = 128,
        np_alloc_factor = 4.0,
        omega_m = 0.30749399999999999,
        prefix = "results-za-nongaussian",
        random_seed = 100,
        read_powerspectrum = "powerspec.txt",
        remove_cosmic_variance = false,
        scalar_amp = 0.000000002441,
        scalar_pivot = 0.002,
        scalar_spectral_index = 0.966700000000000004,
        sigma8 = 0,
        write_lineark = "results-za-nongaussian/lineark",
        write_nonlineark = "results-za-nongaussian/nonlineark",
        write_powerspectrum = "results-za-nongaussian/powerspec",
        write_snapshot = "results-za-nongaussian/fastpm",
        write_whitenoisek = "results-za-nongaussian/whitenoisek",
        za = true,
        args = {
            "standard.lua",
            "za",
            "nongaussian",
        },
        change_pm = {
            0,
        },
        output_redshifts = {
            9.0,
            0.0,
        },
        pm_nc_factor = {
            2,
        },
        time_step = {
            1.0,
        },
    }
    schema.print()
    for i,k in pairs(schema.bind(ns)) do
        print(i, k)
    end
end

fastpm.schema = schema

--
-- The main functions, must be a global symbol
--
function _parse_runmain(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    local globals = setmetatable({}, {__index=_G})

    globals.fastpm = fastpm
    globals.logspace = fastpm.logspace
    globals.linspace = fastpm.linspace
    globals.loglinspace = fastpm.loglinspace

    return config.parse(fastpm.schema, filename, true, globals, {...})
end

function _parse(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    local globals = setmetatable({}, {__index=_G})
    globals.fastpm = fastpm
    globals.logspace = fastpm.logspace
    globals.linspace = fastpm.linspace
    globals.loglinspace = fastpm.loglinspace

    return config.parse(fastpm.schema, filename, false, globals, {...})
end

function _help(filename, ...)

    local fastpm = require('lua-runtime-fastpm')
    local config = require('lua-runtime-config')

    return fastpm.schema.format_help()
end

return fastpm
