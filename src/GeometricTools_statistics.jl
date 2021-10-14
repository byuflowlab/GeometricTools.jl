#=##############################################################################
# DESCRIPTION
    Methods for statistical analysis on geometric datasets
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Oct 2021
  * License   : MIT License
=###############################################################################

# ----------- CARTESIAN OPERATIONS ---------------------------------------------
_agglomerate_mean(name, din, dout, args...) = din[name]
function _agglomerate_variance(name, din, dout, args...; mean_suf="-mean")
    if length(size(din[name]))==1
        return ( (din[name][i]-dout[name*mean_suf][i])^2
                        for i in 1:length(din[name]))
    else
        return ( (din[name][i,j]-dout[name*mean_suf][i,j])^2
                       for i in 1:size(din[name],1), j in 1:size(din[name], 2) )
    end
end

"Correlation of k-th direction with all other directions"
function _agglomerate_correlation(k::Int, name, din, dout, args...; mean_suf="-mean")
    return (
                (din[name][k, j]-dout[name*mean_suf][k, j])*(din[name][(i-1)%3+1, j]-dout[name*mean_suf][(i-1)%3+1, j])
                        for i in 1:size(din[name], 1), j in 1:size(din[name], 2)
           )
end
_agglomerate_correlation1(args...; optargs...) = _agglomerate_correlation(1, args...; optargs...)
_agglomerate_correlation2(args...; optargs...) = _agglomerate_correlation(2, args...; optargs...)
_agglomerate_correlation3(args...; optargs...) = _agglomerate_correlation(3, args...; optargs...)


# ----------- CYLINDRICAL OPERATIONS -------------------------------------------
"Calculate radial and azimuthal direction at each centroid"
function _get_rad_azim_dirs(axial_dir, din)

    keyi = findfirst(x->contains(x, "centroid"), collect(keys(din)))
    if keyi == nothing
        error("Centroids not found."*
        " Make sure that centroids are calculated before calling this function")
    end

    # Fetch position of each element
    centroid_key = collect(keys(din))[keyi]
    centroids = din[centroid_key]
    nelems = size(centroids, 2)

    # Determine radial directions
    radial_dirs_aux = ( ((centroids[i, j] - dot(view(centroids, :, j), axial_dir)*axial_dir[i])
                                            for i in 1:3) for j in 1:nelems )
    radial_dirs = ( collect(R)/norm(R) for R in radial_dirs_aux )

    # Determine azimuthal directions
    azim_dirs = ( cross(axial_dir, R) for R in radial_dirs )

    return radial_dirs, azim_dirs
end

function _agglomerate_cyl_mean(axial_dir::Array{<:Real, 1}, name, din, args...)

    radial_dirs, azim_dirs = _get_rad_azim_dirs(axial_dir, din)

    if length(size(din[name]))==1
        return error("Cylindrical transformation requested on a scalar field!")
    end

    data = din[name]
    ndims = size(data, 1)
    nelems = size(data, 2)
    axial_dirs = (axial_dir for i in 1:nelems)

    return (
                dot( view(data, :, j[1]), j[i+1] )
                    for i in 1:ndims,
                        j in zip(1:nelems, radial_dirs, azim_dirs, axial_dirs)
            )
end

function _agglomerate_cyl_variance(axial_dir::Array{<:Real, 1}, name, din, dout, args...; mean_suf="-cyl-mean")

    radial_dirs, azim_dirs = _get_rad_azim_dirs(axial_dir, din)

    if length(size(din[name]))==1
        return error("Cylindrical transformation requested on a scalar field!")
    end

    Din = din[name]
    Dout = dout[name*mean_suf]
    ndims = size(Din, 1)
    nelems = size(Din, 2)
    axial_dirs = (axial_dir for i in 1:nelems)

    return (
                (dot( view(Din, :, j[1]), j[i+1] ) - Dout[i,j[1]])^2
                    for i in 1:ndims,
                        j in zip(1:nelems, radial_dirs, azim_dirs, axial_dirs)
            )
end

"Correlation of k-th direction with all other directions"
function _agglomerate_cyl_correlation(k::Int, axial_dir::Array{<:Real, 1}, name, din, dout, args...; mean_suf="-cyl-mean")

    radial_dirs, azim_dirs = _get_rad_azim_dirs(axial_dir, din)

    if length(size(din[name]))==1
        return error("Cylindrical transformation requested on a scalar field!")
    end

    Din = din[name]
    Dout = dout[name*mean_suf]
    ndims = size(Din, 1)
    nelems = size(Din, 2)
    axial_dirs = (axial_dir for i in 1:nelems)

    return (
                (dot( view(Din, :, j[1]), j[k+1] ) - Dout[k,j[1]]) * (dot( view(Din, :, j[1]), j[(i-1)%3+1+1] ) - Dout[(i-1)%3+1, j[1]])
                    for i in 1:ndims,
                        j in zip(1:nelems, radial_dirs, azim_dirs, axial_dirs)
            )
end

_agglomerate_cyl_correlation1(args...; optargs...) = _agglomerate_cyl_correlation(1, args...; optargs...)
_agglomerate_cyl_correlation2(args...; optargs...) = _agglomerate_cyl_correlation(2, args...; optargs...)
_agglomerate_cyl_correlation3(args...; optargs...) = _agglomerate_cyl_correlation(3, args...; optargs...)




# ----------- STATISTICAL OPERATIONS ON VTKS -----------------------------------
function calculate_statistics_vtk(vtkfiles;                  # List of files to process
                                    read_path="",            # Folder where to read files from
                                    ites=2,                  # Number of processing iterations
                                                             # Operations to perform; (operation name, operation function, iteration to act on)
                                    operations=[("mean", _agglomerate_mean, 1), ("variance", _agglomerate_variance, 2)],
                                    add_centroid=true,       # Whether to calculate and operate on centroid
                                    verbose=true, v_lvl=0)

    data_out = Dict()
    ncells = nothing
    centroid_cells = nothing

    nfiles = length(vtkfiles)

    for ite in 1:ites

        if verbose
           println("\t"^(v_lvl)*"Iteration $(ite) out of $(ites)")
        end

        for (fi, filename) in enumerate(vtkfiles)

            if verbose && fi%(ceil(Int, nfiles/10))==0
               println("\t"^(v_lvl+1)*"Reading file $(fi) out of $(nfiles)")
            end

            # Read VTK
            points, cells, cell_types, data = read_vtk(filename; path=read_path)


            # Pre-allocate and initialize output data with zeroes
            if ite==1 && fi==1
                ncells = length(cells)
                centroid_cells = zeros(3, ncells)
            end

            # Add centroid to original dataset
            if add_centroid

                # Point centroid
                if !("POINT_DATA" in keys(data))
                    data["POINT_DATA"] = Dict()
                end
                data["POINT_DATA"]["centroid-point"] = points

                # Cell centroid
                if length(cells) != 0

                    centroid_cells .= 0
                    for (celli, cell) in enumerate(cells)
                        for pi in cell
                            centroid_cells[:, celli] .+= view(points, :, pi+1)
                        end
                        centroid_cells[:, celli] ./= ncells
                    end

                    if !("CELL_DATA" in keys(data))
                        data["CELL_DATA"] = Dict()
                    end
                    data["CELL_DATA"]["centroid-cell"] = centroid_cells
                end
            end

            # Pre-allocate and initialize output data with zeroes
            if ite==1 && fi==1

                for (dataparent, subdata) in data
                    data_out[dataparent] = Dict()

                    for (name, this_data) in subdata
                        for (op_name, op) in operations
                            data_out[dataparent][name*"-"*op_name] = zero(this_data)
                        end
                    end
                end

            end

            # Iterate over data sets
            for (dataparent, subdata) in data
                for (name, this_data) in subdata

                    # Apply operation to agglomerate each dataset
                    for (op_name, op, op_ite) in operations
                        if op_ite==ite
                            try
                                data_out[dataparent][name*"-"*op_name] .+= op(name, subdata, data_out[dataparent], points, cells)
                            catch e
                                println(size(data_out[dataparent][name*"-"*op_name]))
                                println(size(collect(op(name, subdata, data_out[dataparent], points, cells))))
                                println()
                                throw(e)
                            end
                        end
                    end

                end
            end

            # Convert the agglomerates into averages
            if fi==nfiles
                for (dataparent, subdata) in data
                    for (name, this_data) in subdata
                        for (op_name, op, op_ite) in operations
                            if op_ite==ite
                                data_out[dataparent][name*"-"*op_name] ./= nfiles
                            end
                        end
                    end
                end
            end

        end

    end

    return data_out
end


function calculate_statistics_vtk(vtkfiles, save_path;
                                    read_path="",
                                    out_filename="statistics",
                                    verbose=true, v_lvl=0,
                                    optargs...)

    data_out = calculate_statistics_vtk(vtkfiles;
                                        read_path=read_path,
                                        verbose=verbose, v_lvl=v_lvl,
                                        optargs...)
    str = nothing

    if save_path != nothing

        if verbose
           println("\t"^(v_lvl)*"Saving statistics under '$(out_filename).vtk'")
        end

        # Read last VTK as default geometry
        points, cells, cell_types, data = read_vtk(vtkfiles[end]; path=read_path)

        # Re-format points
        points = [points[:, pi] for pi in 1:size(points, 2)]

        point_data, cell_data = nothing, nothing

        if "POINT_DATA" in keys(data_out)
            point_data = [
                            Dict(
                                "field_name" => name,
                                "field_type" => length(size(this_data))==1 ? "scalar" : "vector",
                                "field_data" => length(size(this_data))==1 ? this_data : [this_data[:, pi] for pi in 1:size(this_data, 2)]
                            ) for (name, this_data) in data_out["POINT_DATA"]
                        ]
        end

        if "CELL_DATA" in keys(data_out)
            cell_data = [
                            Dict(
                                "field_name" => name,
                                "field_type" => length(size(this_data))==1 ? "scalar" : "vector",
                                "field_data" => length(size(this_data))==1 ? this_data : [this_data[:, celli] for celli in 1:size(this_data, 2)]
                            ) for (name, this_data) in data_out["CELL_DATA"]
                        ]
        end

        # Generate VTK file
        str = generateVTK(out_filename, points;
                                cells=cells,
                                point_data=point_data,
                                cell_data=cell_data,
                                path=save_path)
    end

    return data_out, str
end
