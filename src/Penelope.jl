module Penelope

include("hdf5.jl")

using Printf

export Particle,
    Interaction

const penelope_so="deps/penelope.so"
const MaxBodies=5000
# MaxBodies must equal both NBV in PENVARED_mod and NB in module PENGEOM_mod, or numerous
# fixed-size arrays in Fortran will be mis-sized in Julia.

"Particle types"
@enum Particle begin
    InvalidPart=0
    Electron=1
    Photon=2
    Positron=3
end
const NParticleTypes=3

"Electron/Positron interaction types"
@enum Interaction begin
    InvalidInt=0
    Hinge=1
    HardElastic=2
    HardInelastic=3
    Brem=4
    InnerIon=5
    Annihilation=6
    Delta=7
    Aux=8
    XSplit=9
end

"Photon interaction types"
@enum PhotonInteraction begin
    Rayleigh=1
    Compton=2
    PhotoElecAbs=3
    PairProduction=4
    XDelta=7
    XAux=8
    Split=9
end
const NInteractionTypes=8

greet() = print("Hello World!")

# Global variables from the Penelope Fortran modules
mutable struct Track
    ptr_kpar::Ptr{Int32}
    ptr_ibody::Ptr{Int32}
    ptr_mat::Ptr{Int32}
    ptr_ilb::Ptr{Int32}
    ptr_e::Ptr{Float64}
    ptr_x::Ptr{Float64}
    ptr_y::Ptr{Float64}
    ptr_z::Ptr{Float64}
    ptr_u::Ptr{Float64}
    ptr_v::Ptr{Float64}
    ptr_w::Ptr{Float64}
    ptr_wt::Ptr{Float64}
    ptr_e0segm::Ptr{Float64}
    ptr_desoft::Ptr{Float64}
    ptr_ssoft::Ptr{Float64}
end
Track() = Track(zeros(Int, 15)...)
gtrack = Track()

energy(t::Track=gtrack) = unsafe_load(t.ptr_e)
weight(t::Track=gtrack) = unsafe_load(t.ptr_wt)
location(t::Track=gtrack) = [unsafe_load(t.ptr_x), unsafe_load(t.ptr_y), unsafe_load(t.ptr_z)]
direction(t::Track=gtrack) = [unsafe_load(t.ptr_u), unsafe_load(t.ptr_v), unsafe_load(t.ptr_w)]
particle(t::Track=gtrack) = Particle(unsafe_load(t.ptr_kpar))
body(t::Track=gtrack) = unsafe_load(t.ptr_ibody)
material(t::Track=gtrack) = unsafe_load(t.ptr_mat)
ilb(t::Track=gtrack) = [unsafe_load(t.ptr_ilb, i) for i=1:5]

function cause(t::Track=gtrack)
    v = [unsafe_load(t.ptr_ilb, i) for i=1:4]
    part = Particle(v[2])
    interact = Interaction(v[3])
    if part==Photon
        interact = PhotonInteraction(v[3])
    end
    v[1], part, interact, v[4]
end


mutable struct VarianceReduction
    ptr_force::Ptr{Float64}
    ptr_wforce::Ptr{Float64}
    ptr_dea::Ptr{Float64}
    ptr_bremsplit::Ptr{Int32}
    wtmin::Array{Float64,2}  # WLOW in Fortran
    wtmax::Array{Float64,2}  # WHIG in Fortran
    forcing::Array{Bool,2} # LFORCE in Fortran
end
VarianceReduction() = VarianceReduction(zeros(Int, 4)..., zeros(Float64, MaxBodies, NParticleTypes),
    ones(Float64, MaxBodies, NParticleTypes), zeros(Bool, MaxBodies, NParticleTypes))
varred = VarianceReduction()

"""
    set_forcing([vr::VarianceReduction], body, particle, interaction,
    forcefactor, [wtmin, [wtmax]])

Set the forcing factor `forcefactor` for the given `body`, `particle`, and `interaction`.
Set the window [`wtmin`, `wtmax`] for the particle weights that are subject to forcing.
If not given, `wtmin`=0 and `wtmax`=1e6.
"""
function set_forcing(vr::VarianceReduction, body::Integer, particle::Particle, interaction::Interaction,
    forcefactor::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)

    # The forcefactor<0 is a "trick" described in the penelope manual.
    if forcefactor < 0
        return set_forcing_per_path(body, particle, interaction, -forcefactor, wtmin, wtmax)
    end

    if body<1 || body>MaxBodies
        @error "body=$body, require 1 ≤ body ≤ $MaxBodies"
    end
    ipart = Int(particle)
    if forcefactor<1
        @error "forcefactor=$forcefactor, require forcefactor ≥ 1"
    end

    vr.wtmin[body, ipart] = wtmin
    vr.wtmax[body, ipart] = wtmax
    vr.forcing[body, ipart] = true
    idx = body + MaxBodies*(ipart-1 + NParticleTypes*(Int(interaction)-1))
    unsafe_store!(vr.ptr_force, Float64(forcefactor), idx)
    nothing
end
set_forcing(body::Integer, particle::Particle, interaction::Interaction,
    forcefactor::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)=set_forcing(
        varred, body, particle, interaction, forcefactor, wtmin, wtmax)
set_forcing(body::Integer, particle::Particle, interaction::PhotonInteraction,
    forcefactor::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)=set_forcing(
        varred, body, particle, Interaction(Int(interaction)), forcefactor, wtmin, wtmax)

"""
    set_forcing_per_path(body, particle, interaction,
    interactperpath, [wtmin, [wtmax]])

Set the interaction forcing factor to create `interactperpath` interactions on average.

Set the forcing factor for the given `body`, `particle`, and `interaction` such that there are
`interactperpath` interactions per photon mean free path (for photons) or per mean distance required
to slow down to rest from Emax. Also set the window [`wtmin`, `wtmax`] for the particle weights that
are subject to forcing. If not given, `wtmin`=0 and `wtmax`=1e6.

This corresponds to a choice of negative forcing factor in Fortran Penelope.
"""
function set_forcing_per_path(body::Integer, particle::Particle, interaction::Interaction,
    interactperpath::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)
    @assert interactperpath ≥ 0
    emax = 15e3 # TODO what about emax??
    imat = 1 # TODO find the material number from the body number
    mean_n_interact = ccall((:avncol_, penelope_so), Cdouble, (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}),
        emax, Int32(particle), imat, Int32(interaction))
    if mean_n_interact > 1e-8
        forcefactor = interactperpath/mean_n_interact
    else
        forcefactor = interactperpath
    end
    forcefactor = max(forcefactor, 1.0)
    set_forcing(body, particle, interaction, forcefactor, wtmin, wtmax)
end
set_forcing_per_path(body::Integer, particle::Particle, interaction::PhotonInteraction,
    interactperpath::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)=set_forcing_per_path(
        varred, body, particle, Interaction(Int(interaction)), interactperpath, wtmin, wtmax)

function __init__()
    gtrack.ptr_kpar = cglobal((:__track_mod_MOD_kpar, penelope_so), Int32)
    gtrack.ptr_ibody = cglobal((:__track_mod_MOD_ibody, penelope_so), Int32)
    gtrack.ptr_mat = cglobal((:__track_mod_MOD_mat, penelope_so), Int32)
    gtrack.ptr_ilb = cglobal((:__track_mod_MOD_ilb, penelope_so), Int32)
    gtrack.ptr_e = cglobal((:__track_mod_MOD_e, penelope_so), Float64)
    gtrack.ptr_u = cglobal((:__track_mod_MOD_u, penelope_so), Float64)
    gtrack.ptr_v = cglobal((:__track_mod_MOD_v, penelope_so), Float64)
    gtrack.ptr_w = cglobal((:__track_mod_MOD_w, penelope_so), Float64)
    gtrack.ptr_x = cglobal((:__track_mod_MOD_x, penelope_so), Float64)
    gtrack.ptr_y = cglobal((:__track_mod_MOD_y, penelope_so), Float64)
    gtrack.ptr_z = cglobal((:__track_mod_MOD_z, penelope_so), Float64)
    gtrack.ptr_wt = cglobal((:__track_mod_MOD_wght, penelope_so), Float64)
    gtrack.ptr_e0segm = cglobal((:__track_mod_MOD_e0segm, penelope_so), Float64)
    gtrack.ptr_desoft = cglobal((:__track_mod_MOD_desoft, penelope_so), Float64)
    gtrack.ptr_ssoft = cglobal((:__track_mod_MOD_ssoft, penelope_so), Float64)

    varred.ptr_force = cglobal((:__penvared_mod_MOD_force, penelope_so), Float64)
    varred.ptr_wforce = cglobal((:__penvared_mod_MOD_wforce, penelope_so), Float64)
    varred.ptr_dea = cglobal((:__penvared_mod_MOD_dea, penelope_so), Float64)
    varred.ptr_bremsplit = cglobal((:__penvared_mod_MOD_ibrspl, penelope_so), Int32)
end

"""
    flattenstrings(s, strlen)

Flatten the array of strings `s` into a single string of length `strlen` per element, padding each with spaces.
If any element of `s` is longer than `strlen` it will be truncated. Return the flattened string."""
function flattenstrings(s::Vector{T}, strlen::Integer) where {T<:AbstractString}
    nstring = length(s)
    for i=1:nstring
        L = length(s[i])
        if L > strlen
            s[i] = s[i][1:strlen]
        elseif L < strlen
            s[i] = s[i] * " "^(strlen-L)
        end
    end
    out = join(s)
end

"""Represent key simulation parameters for a material"""
struct MaterialParams
    ElectronAbs::Float64
    PhotonAbs::Float64
    PositronAbs::Float64
    C1::Float64
    C2::Float64
    CutoffHard::Float64
    CutoffBrem::Float64
    Filename::String
    function MaterialParams(a1, a2, a3, c1, c2, wcc, wcr, f)
        c1>0.2 && error("c1 must be ≤ 0.2")
        c2>0.2 && error("c2 must be ≤ 0.2")
        new(a1, a2, a3, c1, c2, wcc, wcr, f)
    end
end
MaterialParams() = MaterialParams(1000, 1000, 1000, 0.1, 0.1, 1000, 1000, "Cu.mat")

function RavenParams()
    m1 = MaterialParams(1000, 1000, 1000, 0.1, 0.1, 1000, 1000, "Ti.mat")
    m2 = MaterialParams(1000, 1000, 1000, 0.1, 0.1, 1000, 1000, "SiN.mat")
    m3 = MaterialParams(1000, 1000, 1000, 0.1, 0.1, 1000, 1000, "SiO2.mat")
    [m1, m2, m3]
end

"""
    fortranify(p)

Convert `p` (a MaterialParams object or a vector of them) to the vector of floats that
Penelope::MINITW requires."""
function fortranify(p::Vector{MaterialParams})
    n = length(p)
    v = Array{Float64}(undef, 7n)
    for (i, sp)= enumerate(p)
        v[7i-6] = sp.ElectronAbs
        v[7i-5] = sp.PhotonAbs
        v[7i-4] = sp.PositronAbs
        v[7i-3] = sp.C1
        v[7i-2] = sp.C2
        v[7i-1] = sp.CutoffHard
        v[7i-0] = sp.CutoffBrem
    end
    v
end
fortranify(p::MaterialParams) = fortranify([p])

mutable struct Control
    e_beam::Float64
    beam_location::Vector{Float64}
    beam_direction::Vector{Float64}
    seed::Int32
    geofile::String
    materials::Vector{MaterialParams}
    n_materials::Int
    n_bodies::Int
    body2material::Vector{Int}
    xrf_split::Vector{Int}
    brem_split::Vector{Int}
    function Control(e, bl, bd, s, gf, m)
        Nmat = length(m)
        new(e, bl, bd, s, gf, m, Nmat, 0, Int[], Int[], Int[])
    end
end
default_control(seed) = Control(15e3, [0, 0, 1], [0, 0, -1], Int32(seed), "epma1.geo", [MaterialParams()])
raven_control(seed, θ=0.0) = Control(15e3, [-sin(θ), 0, cos(θ)], [sin(θ), 0, -cos(θ)], Int32(seed), "epmaRaven.geo", RavenParams())

"""
    setup_penelope(seed)

Set up the Penelope program. Initialize random number with `seed`"""
function setup_penelope(seed::Integer, Emax::Real, MaterialParams::Vector{MaterialParams}, geofile::AbstractString)
    reset_counters()

    materialsFiles = [m.Filename for m in MaterialParams]
    Nmaterials = length(MaterialParams)
    flatfiles = flattenstrings(materialsFiles, 20)
    @assert length(MaterialParams) == Nmaterials

    seed32 = Int32(seed)
    Emax = Float64(Emax)

    ccall((:rand0_, penelope_so), Cvoid, (Ref{Int32},), seed32)
    ccall((:minitw_, penelope_so), Cvoid, (Ref{Int32}, Ptr{Float64}), Nmaterials, fortranify(MaterialParams))

    infolevel = Int32(3)
    ccall((:pinitw_, penelope_so), Cvoid, (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{UInt8}), Emax, Nmaterials, infolevel, flatfiles)

    geominput = zeros(Float64, 10)
    null = Ref{Int32}(0)
    NmaterialsG = Ref{Int32}()
    Nbodies = Ref{Int32}()
    geofile_mod = flattenstrings([geofile], 20)
    ccall((:ginitw_, penelope_so), Cvoid, (Ptr{Float64}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Ref{Int32}),
        geominput, null, geofile_mod, NmaterialsG, Nbodies)
    @printf("GINITW called. Returns %d materials and %d bodies\n", NmaterialsG[], Nbodies[])
    @assert Nmaterials == NmaterialsG[]

    # Tranlation from body number to material number, at least within PenGeom.
    ptr = cglobal((:__pengeom_mod_MOD_mater, Penelope.penelope_so), Int32)
    body2material = [unsafe_load(ptr, i) for i=1:Nbodies[]]

    Nmaterials, Nbodies[], body2material
end
setup_penelope(seed::Integer, Emax::Real, MaterialParams::MaterialParams, geofile::AbstractString) = setup_penelope(seed, Emax, [MaterialParams], geofile)
function setup_penelope(c::Control)
    nm, nb, b2m = setup_penelope(c.seed, c.e_beam, c.materials, c.geofile)
    @assert nm == c.n_materials
    c.n_bodies = nb
    c.body2material = b2m

    if nb > length(c.brem_split)
        append!(c.brem_split, ones(Int, nb-length(c.brem_split)))
    end
    if nb > length(c.xrf_split)
        append!(c.xrf_split, ones(Int, nb-length(c.xrf_split)))
    end
    @assert nb == length(c.brem_split)
    for i=1:nb
        unsafe_store!(varred.ptr_bremsplit, c.brem_split[i], i)
    end
    nm, nb, b2m
end

"""
    initialize_track(Eprim::Real, location::Vector, direction::Vector)

Initialize an electron track with energy `Eprim`, 3d location `location` and direction cosines `direction`.
"""
function initialize_track(Eprim::Real, location::Vector, direction::Vector)
    unsafe_store!(gtrack.ptr_kpar, Int32(Electron))
    unsafe_store!(gtrack.ptr_e, Eprim)
    unsafe_store!(gtrack.ptr_x, location[1])
    unsafe_store!(gtrack.ptr_y, location[2])
    unsafe_store!(gtrack.ptr_z, location[3])
    unsafe_store!(gtrack.ptr_u, direction[1])
    unsafe_store!(gtrack.ptr_v, direction[2])
    unsafe_store!(gtrack.ptr_w, direction[3])
    unsafe_store!(gtrack.ptr_wt, 1.0)
    unsafe_store!(gtrack.ptr_ibody, 0)
    ilb = [1, 0, 0, 0, 0]  # signifies a primary particle.
    for i=1:5
        unsafe_store!(gtrack.ptr_ilb, ilb[i], i)
    end
end
initialize_track(c::Control) = initialize_track(c.e_beam, c.beam_location, c.beam_direction)

"Impose soft energy loss on electron/positron track `t` given `distance` moved through medium."
function softEloss(distance::Real, t::Track=gtrack)
    particle(t) == Photon && return
    e = unsafe_load(t.ptr_e0segm)
    e -= unsafe_load(t.ptr_ssoft)*distance
    unsafe_store!(t.ptr_e, e)
end


sumwtr = zeros(Float64, 140, 100, 3)
sumwtx = zeros(Float64, 160, 100, 3)
espect = zeros(Float64, 150, 3)
xspect = zeros(Float64, 1500, 5)
totalcounts = 0
function reset_counters()
    sumwtr[:,:,:].=0.0
    sumwtx[:,:,:].=0.0
    espect[:,:].=0.0
    xspect[:,:].=0.0
    global totalcounts = 0
    nothing
end

function run_sim(c::Control, Nelec::Integer)
    eabs = 1000.0 .+ zeros(Float64, 3, 5)  # TODO: read this from penelope_mod
    penergies = Array{Float64}(undef, 0)
    MAXTRACKSTEPS=10000
    alltracks = Array{Float64}(undef, MAXTRACKSTEPS, 7)
    trackstep = 1

    for n=1:Nelec
        global totalcounts += 1
        initialize_track(c)
        locate_track()
        mat = material()

        # If the particle starts outside all material bodies (as expected), propagate it forward into one.
        if mat == 0
            stepdist = 1e20
            take_one_step(stepdist)
            mat = material()
            if mat == 0  # The particle travels 10^20 cm and still doesn't enter the system.
                continue
            end
        end

        # Clean the secondary stack so it can store all children of this primary.
        clean_secstack()
        particles_to_do = 1
        while particles_to_do > 0  # Loop over all particles, first the primary, then any secondaries.
            start_medium()
            part = particle()
            ipart = Int(part)

            # Split fluorescence x rays.
            ibody = body()
            if part == Photon && c.xrf_split[ibody] > 1
                split_xrays(c, c.xrf_split[ibody])
            end

            # Print some info about photons
            mat = material()
            e = energy()
            emission_spot = location()
            emission_dir = direction()
            ρ = sqrt(sum(emission_spot[1:2].^2))
            z = emission_spot[3]
            w = emission_dir[3]
            wt = weight()
            ilbx = cause()
            if part == Photon
                # @printf("X %8.2f eV ρ=%5.0f nm z=%5.0f nm cos(z)=%7.4f  wt=%.6f %s\n", e, ρ*1e7, z*1e7, w, wt, ilbx)
            end

            while true  # Loop over all steps taken by this particle
                startloc = location()
                ibody = body()
                wt = weight()
                forcing = varred.forcing[ibody, ipart] &&
                    wt > varred.wtmin[ibody, ipart] &&
                    wt ≤ varred.wtmax[ibody, ipart]

                # Compute distance to next interaction (of any type)
                maxdist = 1e-4  # TODO: this should depend on material
                interaction_dist = compute_jump_length(maxdist, forcing)
                actualdist, ncross = take_one_step(interaction_dist)
                if part == Electron && trackstep < MAXTRACKSTEPS && energy()>4900
                    endloc = location()
                    dsq = sum((endloc .- startloc).^2)
                    if dsq > 1e-8
                        endloc = startloc + direction() .* 1e-4
                    end
                    alltracks[trackstep,1:3] = startloc .* 1e7 # to nm
                    alltracks[trackstep,4:6] = endloc .* 1e7
                    alltracks[trackstep,7] = energy()
                    trackstep += 1
                end
                if ncross[] > 0  # Crossed from one material to another.
                    softEloss(actualdist[])
                    mat = material()
                    # Tally electron spectra
                    e = energy()
                    if part == Electron && e>0
                        # @show e, direction()[3]
                        ispect = 1
                        if direction()[3]>0; ispect=2; end
                        if mat != 0; ispect=3; end
                        ebin = round(Int, e/100.0)
                        if ebin≥1 && ebin≤150
                            espect[ebin, ispect] += weight()
                        end
                    end

                    if mat == 0  # Particle left enclosure
                        break
                    end
                    if energy() < eabs[ipart, mat]  # Particle below E threshold was absorbed.
                        break
                    end
                    start_medium()
                    continue
                end

                # Simulate next interaction
                de, icol = knock(forcing)
                part = particle()
                ipart = Int(part)
                mat = material()
                e = energy()
                if e < eabs[ipart, mat]
                    break
                end
            end
            # Done tracking 1 particle

            # Tally facts about that particle
            if part == Photon && e > 0 && w<0
                # @show e, emission_spot, emission_dir
                push!(penergies, e)

                ρ = sqrt(sum(emission_spot[1:2].^2))
                z = emission_spot[3]
                w = emission_dir[3]
                gen, ppart, icol, xrfcode = cause()

                map=0
                if div(xrfcode, 1000000) == 22
                    map=1
                elseif e≤4000
                    map=2
                elseif e>4000
                    map=3
                end
                if map > 0
                    rbin = ceil(Int, ρ*2e6)
                    zbin = ceil(Int, -z*2e6)
                    if rbin<1; rbin=1; end
                    if zbin<1; zbin=1; end
                    if rbin ≤ 140 && zbin ≤ 100
                        sumwtr[rbin, zbin, map] += wt
                    end
                    xbin = round(Int, emission_spot[1]*2e6+80.5)
                    if xbin ≥ 1 && xbin ≤ 160 && zbin ≤ 100
                        sumwtx[xbin, zbin, map] += wt
                    end
                end
                ebin = round(Int, e/10.0)
                if w<0 && ebin≥1 && ebin≤1500
                    if xrfcode == 0
                        xspect[ebin, 1] += wt
                        if z > -70e-7
                            xspect[ebin, 3] += wt
                        elseif z > -2070e-7
                            xspect[ebin, 4] += wt
                        else
                            xspect[ebin, 5] += wt
                        end
                    else
                        xspect[ebin, 2] += wt
                    end
                end
            end

            particles_to_do = number_secondaries()
        end
        # Done tracking all particles from this primary's secondary stack.

    end
    # Done tracking Nprimaries primary particles.

    alltracks[1:trackstep,:]
end

"Generate isotropic random direction cosines in 3D just as Penelope does"
function isotropic_random()
    # WS=-1.0D0+2.0D0*RAND(10.0D0)
    # SDTS=SQRT(1.0D0-WS*WS)
    # DF=TWOPI*RAND(11.0D0)
    # US=COS(DF)*SDTS
    # VS=SIN(DF)*SDTS
    dummy = 10.0
    r1 = ccall((:rand_, penelope_so), Cdouble, (Ref{Float64}, ), dummy)
    r2 = ccall((:rand_, penelope_so), Cdouble, (Ref{Float64}, ), dummy)
    cosθ = 2r1-1.0
    ϕ = 2π*r2
    sinθ = sqrt(1-cosθ^2)
    cos(ϕ)*sinθ, sin(ϕ)*sinθ, cosθ
end


function split_xrays(c::Control, splitfactor::Int)
    part = particle()
    generation, parent, interaction, xrftype = cause()
    # Particle must be a 2nd generation XRF photon, not the result of splitting
    if (part == Photon &&
        splitfactor > 1 &&
        generation ==2 &&
        xrftype > 0 &&
        interaction != XSplit)
        wt = weight()/splitfactor
        ilbx = ilb()
        ilbx[3] = Int32(Split)
        unsafe_store!(gtrack.ptr_wt, wt)
        unsafe_store!(gtrack.ptr_ilb, Int32(Split), 3)
        x, y, z = location()
        for i=2:splitfactor
            u, v, w = isotropic_random()
            ccall((:stores_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int32}, Ptr{Int32}, Ref{Int32}),
            energy(), x, y, z, u, v, w, wt, Int32(Photon), ilbx, 0)
        end
    end
end

"""Simulate an interaction (by wrapping KNOCK). Return (`de`, `icol`)=(deposited energy, kind of event).
See Penelope Table 7.5 for the `icol` convention."""
function knock(forcing::Bool)
    de = Ref{Float64}()
    icol = Ref{Int32}()
    if forcing
        ccall((:knockf_, penelope_so), Cvoid, (Ref{Float64}, Ref{Int32}), de, icol)
    else
        ccall((:knock_, penelope_so), Cvoid, (Ref{Float64}, Ref{Int32}), de, icol)
    end
    de[], Interaction(icol[])
end

"""Return the number of secondaries remaining on the secondary stack (wrap SECPAR)."""
function number_secondaries()
    nleft = Ref{Int32}()
    ccall((:secpar_, penelope_so), Cvoid, (Ref{Int32},), nleft)
    # TODO: stop at 6th generation??
    nleft[]
end

"""
    compute_jump_length(maxdist, forcing)

Compute length of one forward step for the tracked particle, up to a maximum distance of `maxdist`,
with or without interaction forcing according to `forcing` (wrap JUMP and JUMPF).
Return `actual`, the actual distance traveled.

Contains the interaction physics, but assumes materials of infinite size. Use `take_one_step` to consider
the geometrical limitations of the PenGeom setup."""
function compute_jump_length(maxdist::Real, forcing::Bool)
    actualdist = Ref{Float64}(0)
    if forcing
        ccall((:jumpf_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}),
            Float64(maxdist), actualdist)
    else
        ccall((:jump_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}),
            Float64(maxdist), actualdist)
    end
    actualdist[]
end

"""
    take_one_step(maxdist)

Take one forward step for the tracked particle, up to a maximum distance of `maxdist` or until crossing
an interface from one body to anther (wrap STEP). Return (`actual`, `ncross`) where `actual` is the actual distance
and `ncross` is the number of interface crossings (should be 0 or 1).

Contains no interaction physics, but knows the geometrical limitations of the PenGeom setup.
Use `compute_jump_length` before this to capture the interaction physics."""
function take_one_step(maxdist::Real)
    actualdist = Ref{Float64}(0)
    ncross = Ref{Int32}(0)
    ccall((:step_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Int32}),
        Float64(maxdist), actualdist, ncross)
    actualdist[], ncross[]
end


"Clean the stack of unprocessed secondary particles (wrap CLEANS)."
clean_secstack() = ccall((:cleans_, penelope_so), Cvoid, ())

"Start sim in a new medium (wrap START)."
start_medium() = ccall((:start_, penelope_so), Cvoid, ())

"Determine body and material for a track, in PenGeom (wrap LOCATE)."
function locate_track()
    ccall((:locate_, penelope_so), Cvoid, ())
    unsafe_load(gtrack.ptr_ibody), unsafe_load(gtrack.ptr_mat)
end
function locate_track(x::Real, y::Real, z::Real)
    unsafe_store!(gtrack.ptr_x, x)
    unsafe_store!(gtrack.ptr_y, y)
    unsafe_store!(gtrack.ptr_z, z)
    locate_track()
end
locate_track(xyz) = locate_track(xyz...)

function example1(seed::Int)
    c = default_control(seed)
    c.brem_split[1] = 1
    c.xrf_split[1] = 1
    setup_penelope(c)
    set_forcing_per_path(1, Electron, Brem, 5, 0.9, 1.0)
    set_forcing_per_path(1, Electron, InnerIon, 50, 0.9, 1.0)
    set_forcing(1, Photon, Compton, 10, 1e-3, 1.0)
    set_forcing(1, Photon, PhotoElecAbs, 10, 1e-3, 1.0)
    run_sim(c, 50)
    c
end

function example(seed::Int, θ=0.0)
    c = raven_control(seed, θ)
    c.brem_split = ones(Int, 5)
    c.xrf_split = ones(Int, 5)
    c.brem_split[1] = 1
    c.xrf_split[1] = 1
    setup_penelope(c)
    for body=1:3
        set_forcing_per_path(body, Electron, Brem, 5, 0.9, 1.0)
        set_forcing_per_path(body, Electron, InnerIon, 10, 0.9, 1.0)
        set_forcing(body, Photon, Compton, 10, 1e-3, 1.0)
        set_forcing(body, Photon, PhotoElecAbs, 10, 1e-3, 1.0)
    end
    # set_forcing_per_path(1, Electron, InnerIon, 50, 0.1, 1.0)
    run_sim(c, 1000)
    c
end

end # module
