
"""
    tpt_write(outfile, tpt_result; dir_name, overwrite)

Write `tpt_result` to the file `outfile`, which must be in the `.h5` format.

If `outfile` does not exist, it will be created. Results are written to the directory specified by `dir_name`.

Indices are written to the sub-directory `"indices"` and statistics to the sub-directory `"statistics"`.

### Optional Arguments
- `dir_name`: The name of the directory that `tpt_result` is written to, default `"tpt_homog"`. Directories can be nested, e.g. `dir_name = "trial1/tpt"`.
- `overwrite`: If `true`, the directory `dir_name` will overwrite a directory with the same name if it exists in `outfile`. Default `false`.
"""
function tpt_write(
    outfile::String, 
    tpt_result::AbstractTPTHomogResult;
    dir_name::String = "tpt_homog", 
    overwrite::Bool = false)

    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "cw")

    if dir_name in keys(fout)
        if !overwrite
            @assert !(dir_name in keys(fout)) "This file already has a group with the name: $(dir_name). Pass `overwrite = true` to force a replacement."
        end

        delete_object(fout, dir_name)
    end

    tpt = create_group(fout, dir_name)
    inds = create_group(tpt, "indices")
    stats = create_group(tpt, "statistics")

    for fn in fieldnames(typeof(tpt_result))
        if fn == :sets
            for fns in fieldnames(typeof(tpt_result.sets))
                inds["$fns"] = getfield(tpt_result.sets, fns)
            end
        else
            stats["$fn"] = getfield(tpt_result, fn)
        end
    end

    close(fout)

    @info "TPT results written to $(outfile)."

    return 
end


"""
    tpt_write(outfile, partitions_result; dir_name, overwrite)

Write `partitions_result` to the file `outfile`, which must be in the `.h5` format.

If `outfile` does not exist, it will be created. Results are written to the directory specified by `dir_name`.`.

### Optional Arguments
- `dir_name`: The name of the directory that `partitions_result` is written to, default `"partitions"`. Directories can be nested, e.g. `dir_name = "trial1/partitions"`.
- `overwrite`: If `true`, the directory `dir_name` will overwrite a directory with the same name if it exists in `outfile`. Default `false`.
"""
function tpt_write(
    outfile::String, 
    partitions_result::AbstractPartitionsResult;
    dir_name::String = "partitions", 
    overwrite::Bool = false)

    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "cw")

    if dir_name in keys(fout)
        if !overwrite
            @assert !(dir_name in keys(fout)) "This file already has a group with the name: $(dir_name). Pass `overwrite = true` to force a replacement."
        end

        delete_object(fout, dir_name)
    end

    parts = create_group(fout, dir_name)

    for fn in fieldnames(typeof(partitions_result))
        parts["$fn"] = getfield(partitions_result, fn)
    end

    close(fout)

    @info "Partitions results written to $(outfile)."

    return 
end