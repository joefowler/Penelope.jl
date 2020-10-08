module Penelope

using Printf

const penelope_so="deps/penelope.so"
const MaxBodies=5000
# MaxBodies must equal both NBV in PENVARED_mod and NB in module PENGEOM_mod.

# Particle types
const Electron=1
const Photon=2
const Positron=3
const NParticleTypes=3

# e+/e- interactions
const Hinge=1
const HardElastic=2
const HardInelastic=3
const Brem=4
const InnerIon=5
const Annihilation=6
# photon interactions
const Rayleigh=1
const Compton=2
const PhotoElecAbs=3
const PairProduction=4
# All particle interactions
const Delta=7
const AuxInteract=8
const NInteractionTypes=8
# TODO: should the above const values be @enums?

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
track = Track()

energy(t::Track) = unsafe_load(t.ptr_e)
weight(t::Track) = unsafe_load(t.ptr_wt)
location(t::Track) = [unsafe_load(t.ptr_x), unsafe_load(t.ptr_y), unsafe_load(t.ptr_z)]
direction(t::Track) = [unsafe_load(t.ptr_u), unsafe_load(t.ptr_v), unsafe_load(t.ptr_w)]
particle(t::Track) = unsafe_load(t.ptr_kpar)
body(t::Track) = unsafe_load(t.ptr_ibody)
material(t::Track) = unsafe_load(t.ptr_mat)
ilb(t::Track) = [unsafe_load(t.ptr_ilb, i) for i=1:5]


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
function set_forcing(vr::VarianceReduction, body::Integer, particle::Integer, interaction::Integer,
    forcefactor::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)
    if body<1 || body>MaxBodies
        @error "body=$body, require 1 ≤ body ≤ $MaxBodies"
    end
    if particle<1 || particle>NParticleTypes
        @error "particle=$particle, require 1 ≤ particle ≤ $NParticleTypes"
    end
    if interaction<1 || interaction>NInteractionTypes
        @error "interaction=$interaction, require 1 ≤ interaction ≤ $NInteractionTypes"
    end
    if forcefactor<1
        @error "forcefactor=$forcefactor, require forcefactor ≥ 1"
    end
    # TODO: handle negative forcing factors somehow

    vr.wtmin[body, particle] = wtmin
    vr.wtmax[body, particle] = wtmax
    vr.forcing[body, particle] = true
    idx = body + MaxBodies*(particle-1 + NParticleTypes*(interaction-1))
    unsafe_store!(vr.ptr_force, Float64(forcefactor), idx)
    nothing
end
set_forcing(body::Integer, particle::Integer, interaction::Integer,
    forcefactor::Real, wtmin::Real=0.0, wtmax::Real=1.0e6)=set_forcing(varred, body, particle, interaction, forcefactor, wtmin, wtmax)

function __init__()
    track.ptr_kpar = cglobal((:__track_mod_MOD_kpar, penelope_so), Int32)
    track.ptr_ibody = cglobal((:__track_mod_MOD_ibody, penelope_so), Int32)
    track.ptr_mat = cglobal((:__track_mod_MOD_mat, penelope_so), Int32)
    track.ptr_ilb = cglobal((:__track_mod_MOD_ilb, penelope_so), Int32)
    track.ptr_e = cglobal((:__track_mod_MOD_e, penelope_so), Float64)
    track.ptr_u = cglobal((:__track_mod_MOD_u, penelope_so), Float64)
    track.ptr_v = cglobal((:__track_mod_MOD_v, penelope_so), Float64)
    track.ptr_w = cglobal((:__track_mod_MOD_w, penelope_so), Float64)
    track.ptr_x = cglobal((:__track_mod_MOD_x, penelope_so), Float64)
    track.ptr_y = cglobal((:__track_mod_MOD_y, penelope_so), Float64)
    track.ptr_z = cglobal((:__track_mod_MOD_z, penelope_so), Float64)
    track.ptr_wt = cglobal((:__track_mod_MOD_wght, penelope_so), Float64)
    track.ptr_e0segm = cglobal((:__track_mod_MOD_e0segm, penelope_so), Float64)
    track.ptr_desoft = cglobal((:__track_mod_MOD_desoft, penelope_so), Float64)
    track.ptr_ssoft = cglobal((:__track_mod_MOD_ssoft, penelope_so), Float64)

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
struct SimParams
    ElectronAbs::Float64
    PhotonAbs::Float64
    PositronAbs::Float64
    C1::Float64
    C2::Float64
    CutoffHard::Float64
    CutoffBrem::Float64
    function SimParams(a1, a2, a3, c1, c2, wcc, wcr)
        c1>0.2 && error("c1 must be ≤ 0.2")
        c2>0.2 && error("c2 must be ≤ 0.2")
        new(a1, a2, a3, c1, c2, wcc, wcr)
    end
end
SimParams() = SimParams(1000, 1000, 1000, 0.1, 0.1, 1000, 1000)

"""
    fortranify(p)

Convert `p` (a SimParams object or a vector of them) to the vector of floats that
Penelope::MINITW requires."""
function fortranify(p::Vector{SimParams})
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
fortranify(p::SimParams) = fortranify([p])

"""
    setup_penelope(seed)

Set up the Penelope program. Initialize random number with `seed`"""
function setup_penelope(seed::Integer, Emax::Real, simparams::Vector{SimParams})

    materialsFiles = ["Cu.mat"]
    flatfiles = flattenstrings(materialsFiles, 20)
    Nmaterials = length(materialsFiles)
    @assert length(simparams) == Nmaterials

    seed32 = Int32(seed)
    Emax = Float64(Emax)

    ccall((:rand0_, penelope_so), Cvoid, (Ref{Int32},), seed32)
    ccall((:minitw_, penelope_so), Cvoid, (Ref{Int32}, Ptr{Float64}), Nmaterials, fortranify(simparams))

    infolevel = Int32(3)
    ccall((:pinitw_, penelope_so), Cvoid, (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ptr{UInt8}), Emax, Nmaterials, infolevel, flatfiles)

    geominput = zeros(Float64, 10)
    null = Ref{Int32}(0)
    NmaterialsG = Ref{Int32}()
    Nbodies = Ref{Int32}()
    ccall((:ginitw_, penelope_so), Cvoid, (Ptr{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}), geominput, null, NmaterialsG, Nbodies)
    Nmaterials, Nbodies[]
end
setup_penelope(seed::Integer, Emax::Real, simparams::SimParams) = setup_penelope(seed, Emax, [simparams])

"""
    initialize_track(Eprim::Real, location::Vector, direction::Vector)

Initialize an electron track with energy `Eprim`, 3d location `location` and direction cosines `direction`.
"""
function initialize_track(Eprim::Real, location::Vector, direction::Vector)
    unsafe_store!(track.ptr_kpar, Electron)
    unsafe_store!(track.ptr_e, Eprim)
    unsafe_store!(track.ptr_x, location[1])
    unsafe_store!(track.ptr_y, location[2])
    unsafe_store!(track.ptr_z, location[3])
    unsafe_store!(track.ptr_u, direction[1])
    unsafe_store!(track.ptr_v, direction[2])
    unsafe_store!(track.ptr_w, direction[3])
    unsafe_store!(track.ptr_wt, 1.0)
    unsafe_store!(track.ptr_ibody, 0)
    ilb = [1, 0, 0, 0, 0]  # signifies a primary particle.
    for i=1:5
        unsafe_store!(track.ptr_ilb, ilb[i], i)
    end
end

"Impose soft energy loss on electron/positron track `t` given `distance` moved through medium."
function softEloss(t::Track, distance::Real)
    particle(t) == Photon && return
    e = unsafe_load(t.ptr_e0segm)
    e -= unsafe_load(t.ptr_ssoft)*distance
    unsafe_store!(t.ptr_e, e)
end



function run_sim(Nelec::Integer)
    Eprim = 15e3
    init_location = [0.0, 0.0, 1.0]
    init_direction = [0.0, 0.0, -1.0]
    eabs = 1000.0 .+ zeros(Float64, 3, 5)  # TODO: read this from penelope_mod
    penergies = Array{Float64}(undef, 0)

    for n=1:Nelec
        initialize_track(Eprim, init_location, init_direction)
        locate_track()
        mat = material(track)

        # If the particle starts outside all material bodies (as expected), propagate it forward into one.
        if mat == 0
            stepdist = 1e20
            take_one_step(stepdist)
            mat = material(track)
            if mat == 0  # The particle travels 10^20 cm and still doesn't enter the system.
                continue
            end
        end

        # Clean the secondary stack so it can store all children of this primary.
        clean_secstack()
        particles_to_do = 1
        while particles_to_do > 0  # Loop over all particles, first the primary, then any secondaries.
            start_medium()
            part = particle(track)
            if part == Photon
                mat = material(track)
                e = energy(track)
                emission_spot = location(track)
                emission_dir = direction(track)
                z = emission_spot[3]
                w = emission_dir[3]
                ilbx = ilb(track)
                @show e, z, w, ilbx
            end

            while true  # Loop over all steps taken by this particle
                ibody = body(track)
                wt = weight(track)
                forcing = varred.forcing[ibody, part] && wt > varred.wtmin[ibody, part] &&
                    wt ≤ varred.wtmax[ibody, part]

                # Compute distance to next interaction (of any type)
                maxdist = 1e-4  # TODO: this should depend on material
                interaction_dist = compute_jump_length(maxdist, forcing)
                actualdist, ncross = take_one_step(interaction_dist)
                if ncross[] > 0  # Crossed from one material to another.
                    softEloss(track, actualdist[])
                    mat = material(track)
                    part = particle(track)
                    if mat == 0  # Particle left enclosure
                        break
                    end
                    if energy(track) < eabs[part, mat]  # Particle below E threshold was absorbed.
                        break
                    end
                    start_medium()
                    continue
                end

                # Simulate next interaction
                de, icol = knock(forcing)
                part = particle(track)
                mat = material(track)
                e = energy(track)
                if e < eabs[part, mat]
                    break
                end
            end
            # Done tracking 1 particle

            e = energy(track)
            if part == Photon && e > 0
                # ilbx = [unsafe_load(ptr_ilb, i) for i=1:4]
                @show energy, emission_spot, emission_dir
                push!(penergies, e)
            end

            particles_to_do = number_secondaries()
        end
        # Done tracking all particles from this primary's secondary stack.

    end
    # Done tracking Nprimaries primary particles.

    penergies
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
    de[], icol[]
end

"""Return the number of secondaries remaining on the secondary stack (wrap SECPAR)."""
function number_secondaries()
    nleft = Ref{Int32}()
    ccall((:secpar_, penelope_so), Cvoid, (Ref{Int32},), nleft)
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
locate_track() = ccall((:locate_, penelope_so), Cvoid, ())

end # module
