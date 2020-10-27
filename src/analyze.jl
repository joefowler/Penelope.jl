using PyPlot, HDF5, Penelope



function plotwhere(wx::AbstractArray)
    clf()
    names=["Fluorescence", "Low-E Brem", "High-E Brem"]
    colors=["r", "purple", "b"]
    x = -397.5:5:400
    z = 2.5:5:500
    for i = 1:2
        subplot(4,2,2i-1)
        for j=1:3
            plot(x, sum(wx[:,:,j], dims=2)[:,1], color=colors[j], label=names[j])
        end
        legend()
        xlabel("X (nm)")
    end
    semilogy()

    for i = 1:2
        subplot(4,2,2i)
        for j=1:3
            s = sum(wx[:,:,j], dims=1)[1,:]
            plot(z[1:14], s[1:14], color=colors[j], label=names[j])
            plot(z[15:end], s[15:end], color=colors[j])
        end
        legend()
        xlabel("Z (nm)")
    end
    semilogy()

    subplot(4,1,3)
    imshow(log10.(wx[:,1:14,1]'), extent=(-400, 400, 70, 0))
    subplot(4,1,4)
    imshow(log10.((wx[:,:,2]+wx[:,:,3])'), extent=(-400, 400, 500, 0))
end

plotwhere(filename::AbstractString) = plotwhere(Penelope.loadall(filename))
plotwhere(d::Dict) = plotwhere(d["sumwtx"] ./ d["totalcounts"])

function plotspect(sx::AbstractArray)
    clf()
    e = 5:10:15000
    plot(e, sx[:,1], "b", label="Brem")
    plot(e, sx[:,2], "r", label="XRF")
    plot(e, sx[:,3], "g", label="Brem in Ti")
    plot(e, sx[:,4], "c", label="Brem in SiN")
    plot(e, sx[:,5], "y", label="Brem in SiO\$_2\$")
    semilogy()
    legend()
    xlabel("Energy (eV)")
    ylabel("Photons per 10 eV per primary electron")
end

plotspect(d::Dict) = plotspect(d["xspect"] ./ d["totalcounts"])
