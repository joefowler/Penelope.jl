using HDF5

function saveall(filename::AbstractString)
    h5open(filename, "w") do file
        @write file sumwtr
        @write file sumwtx
        @write file espect
        @write file xspect
        tcount = [totalcounts]
        @write file tcount
    end
end

function loadall(filename::AbstractString)
    d = Dict()
    h5open(filename, "r") do file
        for name in ("sumwtr", "sumwtx", "espect", "xspect", "tcount")
            d[name] = read(file, name)
        end
    end
    d["totalcounts"] = d["tcount"][1]
    d
end
