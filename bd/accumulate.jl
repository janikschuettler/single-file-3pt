#=
Merging results across independent simulation runs.

`SFBD.bd` returns *unnormalised integer counts*. That is the whole point of this
file: integer histograms from separate runs can be added together exactly, with
no floating-point drift, and normalised once at the very end. So a run of
M = 16,000,000 trajectories is not one job — it is a cluster array of many
independent jobs, each writing into the same HDF5 file, each adding its counts
to whatever is already there.

    C  =  ( sum_j connected_j )  -  ( sum_j Q1_j )( sum_j Q2_j )
          ---------------------     ------------------------------
                 sum_j M_j                  ( sum_j M_j )^2

The functions below are the "add into whatever is already there" part. The
environment (space grid, time grid, coarsening indices) is stored once as file
attributes and checked on every subsequent write, so runs with mismatched
settings cannot silently be merged.

This is the only part of `bd/` that needs a package: HDF5.jl. The simulation
engine itself (`sfbd.jl`) has no dependencies beyond the standard library.

See the cluster job scripts in ../cluster for how this was driven in practice.
=#

using HDF5

"""
True if `path` (a "/"-separated HDF5 group path) already exists in `fid`.
"""
function check_if_path_exists(fid, path)
    path_tokens = split(path, "/")
    if !(path_tokens[1] in keys(fid))
        return false
    else
        p = path_tokens[1]
        for i in 2:length(path_tokens)
            if length(path_tokens[i]) > 0
                if path_tokens[i] in keys(fid["$p"])
                    p = "$p/$(path_tokens[i])"
                else
                    return false
                end
            end
        end
    end
    return true
end

"""
Adds `val` element-wise to the dataset at `path`, creating it if this is the
first run to write there. This is what makes an array job's results accumulate.
"""
function add_to_existing_data(fid, path, val)
    if check_if_path_exists(fid, path)
        write(fid[path], read(fid[path]) .+ val)
    else
        fid[path] = val
    end
end

"""
Same idea for a scalar attribute — used to accumulate the trajectory count M and
the total spent CPU time across all jobs that contributed to a dataset.
"""
function add_to_existing_attribute(fid, path, key, val; func = (a, b) -> a + b)
    if check_if_path_exists(fid, path) && key in keys(attributes(fid[path]))
        val_existing = read(attributes(fid[path])[key])
        delete_attribute(fid[path], key)
        attributes(fid[path])[key] = func(val_existing, val)
    else
        attributes(fid[path])[key] = val
    end
end

"""
Reads back one accumulated run, together with the environment it was run under.
Returns ((connected_nn, propagators_nn), environment); pass those to
SFBD.correlation together with environment["M"] to get the normalised C.
"""
function read_results_h5(filename; type, N, tag, confinement, x, x1)
    run_id = "$type/N=$N/tag=$tag/confinement=$confinement/x=$x"
    environment = Dict()
    fid = h5open("$filename.h5", "r")
    environment["type"] = type
    environment["N"] = N
    environment["tag"] = tag
    environment["confinement"] = confinement
    environment["x"] = x
    environment["x1"] = x1
    for key in keys(attributes(fid))
        environment[key] = read(attributes(fid)[key])
    end
    for key in keys(attributes(fid[run_id]))
        environment[key] = read(attributes(fid[run_id])[key])
    end
    connected_nn = read(fid["$run_id/x1=$x1/connected_nn"])
    g2p_nn = read(fid["$run_id/x1=$x1/g2p_nn"])
    close(fid)
    return (connected_nn, g2p_nn), environment
end

"""
Writes one run's counts into `results/<filename>.h5`, adding to any run already
stored under the same run_id. `env` carries the discretisation, written once as
file attributes and verified on later writes.

    run_id = "<model>/N=<N>/tag=<tag>/confinement=<c>/x=<x>"
"""
function write_results_h5(filename, out, M, execution_time; model, N, tag,
                          confinement, x, x1s, env)
    isdir("results") || mkdir("results")
    path = "results/$filename.h5"

    if isfile(path)
        h5open(path, "r") do fid
            for key in keys(env)
                if key in keys(attributes(fid)) && read(attributes(fid)[key]) != env[key]
                    error("results/$filename.h5 exists but its '$key' differs; " *
                          "refusing to merge runs from different environments")
                end
            end
        end
    end

    h5open(path, "cw") do fid
        run_id = "$model/N=$N/tag=$tag/confinement=$confinement/x=$x"
        for key in keys(env)
            key in keys(attributes(fid)) || (attributes(fid)[key] = env[key])
        end
        for (i_x1, x1) in enumerate(x1s)
            add_to_existing_data(fid, "$run_id/x1=$x1/connected_nn", out[1][i_x1, :, :])
            add_to_existing_data(fid, "$run_id/x1=$x1/g2p_nn", out[2][[1, 1 + i_x1], :])
        end
        add_to_existing_attribute(fid, run_id, "M", M)
        add_to_existing_attribute(fid, run_id, "execution_time", execution_time)
    end

    return path
end
