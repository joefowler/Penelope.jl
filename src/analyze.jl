using PyPlot, HDF5, Penelope



function plotwhere(wx::AbstractArray)
    clf()
    names=["Dot fluorescence", "Other fluorescence"]
    colors=["r", "orange", (0,3,2)][1:2]
    for j=3:16
        push!(names, "$(2*(j-2)) keV")
        push!(colors, ColorMap("viridis")((j-3.0)/13.0))
    end
    names[3:end] .= ""

    x = -1600+5:10:1600
    z = 5:10:7000
    for i = 1:2
        subplot(4,2,2i-1)
        for j=1:16
            pl = sum(wx[:,:,j], dims=2)[:,1]
            plot(x, pl, color=colors[j], label=names[j])
            # rtmeansqx = sqrt(sum(pl.*x.*x)/sum(pl))
            # i==1 && @show j, rtmeansqx
        end
        s = sum(wx[:,:,3:end], dims=(2,3))[:,1,1]
        plot(x, s, color="k", label="All Brem")
        i==1 && legend()
        xlabel("X (nm)")
    end
    semilogy()

    for i = 1:2
        subplot(4,2,2i)
        for j=1:16
            s = sum(wx[:,:,j], dims=1)[1,:]
            plot(z, s, color=colors[j], label=names[j])
            # frac = sum(s[1:14])/sum(s)
            # i==1 && @show j,frac
            # plot(z[1:14], s[1:14], color=colors[j], label=names[j])
            # plot(z[15:end], s[15:end], color=colors[j])
        end
        s = sum(wx[:,:,3:end], dims=(1,3))[1,:,1]
        plot(z, s, color="k")
        xlabel("Z (nm)")
    end
    semilogy()

    subplot(4,1,3)
    imshow(log10.(wx[:,1:10,1]'), extent=(-1600, 1600, 100, 0), vmin=-8, vmax=-5.5) #; colorbar()
    title("Nanodot fluorescence")
    subplot(4,2,7)
    imshow(log10.((wx[:,:,2])'), extent=(-1600, 1600, 7000, 0), vmin=-8, vmax=-5.5) #; colorbar()
    title("All other fluorescence")
    subplot(4,2,8)
    imshow(log10.((sum(wx[:,:,3:end], dims=3)[:,:,1])'), vmin=-8, vmax=-6, extent=(-1600, 1600, 7000, 0)) #; colorbar()
    title("Bremsstrahlung emission")
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
