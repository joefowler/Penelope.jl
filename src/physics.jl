export
    electronMFP,
    photonMFP,
    particleRange,
    bremYield

"""
    electronMFP(imat[, icol], e)

Return the electron mean free path in material `imat` for interaction `icol` at energy `e`.
If `icol` is omitted, the mean free path for all interactions will be returned.
"""
function electronMFP(imat::Integer, icol::Interaction, e::Real)
    mfp = ccall((:phmfp_, penelope_so), Float64, (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}), Float64(e), Int32(Electron), Int32(imat), Int32(icol))
end
electronMFP(imat::Integer, e::Real) = 1.0/sum([1.0/electronMFP(imat, i, e) for i in (HardElastic, HardInelastic, Brem, InnerIon,Annihilation)])

"""
    photonMFP(imat[, icol], e)

Return the photon mean free path in material `imat` for interaction `icol` at energy `e`.
If `icol` is omitted, the mean free path for all interactions will be returned.
"""
function photonMFP(imat::Integer, icol::PhotonInteraction, e::Real)
    mfp = ccall((:phmfp_, penelope_so), Float64, (Ref{Float64}, Ref{Int32}, Ref{Int32}, Ref{Int32}), Float64(e), Int32(Photon), Int32(imat), Int32(icol))
end
photonMFP(imat::Integer, e::Real) = 1.0/sum([1.0/photonMFP(imat, i, e) for i in (Rayleigh, Compton, PhotoElecAbs)])

"""
    particleRange(imat, ipart::Particle, e)

Compute mean range in material `imat` for particle `ipart` at energy `e`.
Energy `e` may be a real number or a vector."""
function particleRange(imat::Integer, ipart::Particle, e::Real)
    r = ccall((:prange_, penelope_so), Float64, (Ref{Float64}, Ref{Int32}, Ref{Int32}), Float64(e), Int32(ipart), Int32(imat))
end

"""
    bremYield(imat, e)

Return the Bremsstrahlung yield in material `imat` for electrons of energy `e`.
Energy `e` may be a real number or a vector."""
function bremYield(imat::Integer, e::Real)
    r = ccall((:ryield_, penelope_so), Float64, (Ref{Float64}, Ref{Int32}, Ref{Int32}), Float64(e), Int32(Electron), Int32(imat))
end
