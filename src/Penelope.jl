module Penelope

using Printf

const penelope_so="deps/penelope.so"
const Electron=1
const Photon=2
const Positron=3


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

"Impose soft energy loss on track given `distance` moved through medium."
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
        ccall((:locate_, penelope_so), Cvoid, ())
        mat = material(track)

        # If the particle starts outside all material bodies (as expected), propagate it forward
        if mat == 0
            stepdist = 1e20
            actualdist = Ref{Float64}(0)
            ncross = Ref{Int32}(0)
            ccall((:step_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Int32}), stepdist, actualdist, ncross)
            mat = material(track)
            if mat == 0  # The particle travels 10^20 cm and still doesn't enter the system.
                continue
            end
        end

        # Clean the secondary stack so it can store all children of this primary.
        ccall((:cleans_, penelope_so), Cvoid, ())
        while true  # Loop over all particles, first the primary, then any secondaries.
            ccall((:start_, penelope_so), Cvoid, ())  # Start sim in a new medium
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
                maxdist = Ref{Float64}(1e-4)  # TODO: this should depend on material
                actualdist = Ref{Float64}(0)
                ccall((:jump_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}), maxdist, actualdist)
                stepdist = Ref{Float64}(actualdist[])
                ncross = Ref{Int32}(0)
                ccall((:step_, penelope_so), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Int32}), stepdist, actualdist, ncross)
                # ncross[]>0 && @show stepdist[], actualdist[], ncross[]
                if ncross[] > 0  # Crossed from one material to another.
                    softEloss(track, actualdist[])
                    mat = material(track)
                    part = particle(track)
                    if mat == 0  # Particle left enclosure
                        break
                    end
                    if energy(track) < eabs[part, mat]
                        break
                    end
                    ccall((:start_, penelope_so), Cvoid, ())  # Start sim in a new medium
                    continue
                end
                de = Ref{Float64}()
                icol = Ref{Int32}()
                ccall((:knock_, penelope_so), Cvoid, (Ref{Float64}, Ref{Int32}), de, icol)
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

            nleft = Ref{Int32}()
            ccall((:secpar_, penelope_so), Cvoid, (Ref{Int32},), nleft)
            if nleft[] == 0
                break
            end

        end
        # Done tracking all particles from this primary's secondary stack.

    end
    # Done tracking Nprimaries primary particles.

    penergies
end


end # module
