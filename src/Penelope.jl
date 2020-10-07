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
end
Track() = Track(0, 0, 0, 0)
track = Track()

# const ptr_kpar = Ref{Ptr{Int32}}(0)
# const ptr_e = Ref{Ptr{Float64}}(0)
ptr_kpar = Ptr{Int32}(0)
ptr_ibody = Ptr{Int32}(0)
ptr_mat = Ptr{Int32}(0)
ptr_ilb = Ptr{Int32}(0)
ptr_e = Ptr{Float64}(0)
ptr_x = Ptr{Float64}(0)
ptr_y = Ptr{Float64}(0)
ptr_z = Ptr{Float64}(0)
ptr_u = Ptr{Float64}(0)
ptr_v = Ptr{Float64}(0)
ptr_w = Ptr{Float64}(0)
ptr_wt = Ptr{Float64}(0)
ptr_e0segm = Ptr{Float64}(0)
ptr_desoft = Ptr{Float64}(0)
ptr_ssoft = Ptr{Float64}(0)
function __init__()
    global ptr_kpar = cglobal((:__track_mod_MOD_kpar, penelope_so), Int32)
    global ptr_ibody = cglobal((:__track_mod_MOD_ibody, penelope_so), Int32)
    global ptr_mat = cglobal((:__track_mod_MOD_mat, penelope_so), Int32)
    global ptr_ilb = cglobal((:__track_mod_MOD_ilb, penelope_so), Int32)
    track.ptr_kpar = cglobal((:__track_mod_MOD_kpar, penelope_so), Int32)
    track.ptr_ibody = cglobal((:__track_mod_MOD_ibody, penelope_so), Int32)
    track.ptr_mat = cglobal((:__track_mod_MOD_mat, penelope_so), Int32)
    track.ptr_ilb = cglobal((:__track_mod_MOD_ilb, penelope_so), Int32)
    global ptr_e = cglobal((:__track_mod_MOD_e, penelope_so), Float64)
    global ptr_u = cglobal((:__track_mod_MOD_u, penelope_so), Float64)
    global ptr_v = cglobal((:__track_mod_MOD_v, penelope_so), Float64)
    global ptr_w = cglobal((:__track_mod_MOD_w, penelope_so), Float64)
    global ptr_x = cglobal((:__track_mod_MOD_x, penelope_so), Float64)
    global ptr_y = cglobal((:__track_mod_MOD_y, penelope_so), Float64)
    global ptr_z = cglobal((:__track_mod_MOD_z, penelope_so), Float64)
    global ptr_wt = cglobal((:__track_mod_MOD_wght, penelope_so), Float64)
    global ptr_e0segm = cglobal((:__track_mod_MOD_e0segm, penelope_so), Float64)
    global ptr_desoft = cglobal((:__track_mod_MOD_desoft, penelope_so), Float64)
    global ptr_ssoft = cglobal((:__track_mod_MOD_ssoft, penelope_so), Float64)
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
function setup_penelope(seed::Integer, Emax::Float64, simparams::Vector{SimParams})

    materialsFiles = ["Cu.mat"]
    flatfiles = flattenstrings(materialsFiles, 20)
    Nmaterials = length(materialsFiles)
    @assert length(simparams) == Nmaterials

    seed32 = Int32(seed)
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
setup_penelope(seed::Integer, Emax::Float64, simparams::SimParams) = setup_penelope(seed, Emax, [simparams])

"""
    initialize_track(Eprim::Real, location::Vector, direction::Vector)

Initialize an electron track with energy `Eprim`, 3d location `location` and direction cosines `direction`.
"""
function initialize_track(Eprim::Real, location::Vector, direction::Vector)
    ptr_kpar = cglobal((:__track_mod_MOD_kpar, penelope_so), Int32)
    unsafe_store!(ptr_kpar, Electron)
    # unsafe_store!(ptr_e, Eprim)
    # unsafe_store!(ptr_x, location[1])
    # unsafe_store!(ptr_y, location[2])
    # unsafe_store!(ptr_z, location[3])
    # unsafe_store!(ptr_u, direction[1])
    # unsafe_store!(ptr_v, direction[2])
    # unsafe_store!(ptr_w, direction[3])
    # unsafe_store!(ptr_wt, 1.0)
    # unsafe_store!(ptr_ibody, 0)
    # ilb = [1, 0, 0, 0, 0]
    # for i=1:5
    #     unsafe_store!(ptr_ilb, ilb[i], i)
    # end
end


end # module
