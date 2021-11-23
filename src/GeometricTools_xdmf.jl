#=##############################################################################
# DESCRIPTION
    Functions for writing/reading/processing XDMF data.

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2021
  * License   : MIT License
=###############################################################################


"""
    `generateXDMF_3Dstructured(filename, nodes::Matrix, NDIVS;
path="", num=nothing, time=nothing, compression_optargs=(compress=3, ),
fields=Dict(), debug=false)`

Save a three-dimensional structured grid of nodes `nodes` and dimensions `NDIVS`

# REQUIRED ARGUMENTS
* `filename::String`            : Prefix of files that will be created.
* `nodes::Matrix`               : 3xN matrix containing the 3D position of nodes
                                    where N is the number of nodes.
* `NDIVS::Array{Int, 1}`        : Number of cells in each dimension, for any
                                    number of dimensions. `prod(NDIVS .+ 1)`
                                    must equal N.


# REQUIRED ARGUMENTS
* `path::String`                : Path where to save files.
* `num::Int`                    : File suffix. Saves files as `filename.num.*`.
* `time::Real`                  : Time signature.
* `compression_optargs`         : Optional arguments passed to
                                    `HDF5.create_dataset(...)` for data
                                    compression.
* `fields::Dict`                : Attribute fields that will be saved with the
                                    mesh. This is a dictionary of the form
                                    `Dict[field_name] = Dict(
                                        "entry_type" => "node" or "cell",
                                        "field_type" => "scalar" or "vector",
                                        "field_data" =>  Matrix or Vector object,
                                    )`
"""
function generateXDMF_3Dstructured(filename, nodes::Matrix,
                                    NDIVS::Union{Array{Int, 1}, Tuple};
                                    path="", num=nothing, time=nothing,
                                    compression_optargs=(compress=3, ),
                                    fields=Dict(),
                                    debug=false)

    nnodes = size(nodes, 2)
    ncells = prod(NDIVS)
    dims = Iterators.reverse(n+1 for n in NDIVS)

    dims_str = prod("$(n)"*" "^(n!=NDIVS[1]+1) for n in dims)

    # ERROR CASES
    if size(nodes, 1) != 3
        error("Invalid node matrix. Expected size (3, x); got $(size(nodes)).")
    elseif prod(dims) != nnodes
        error("Invalid structured mesh: NDIVS=$(NDIVS) is inconsistent with number of nodes $(nnodes).")
    end

    # File names
    aux = num!=nothing ? ".$num" : ""
    h5fname = filename*aux*".h5"
    xmffname = filename*aux*".xmf"

    # Create path if needed
    if !ispath(path)
        mkpath(path)
    end

    xdmf_fields = Dict()

    # ------------ Save grid in HDF5 format ---------------------
    HDF5.h5open(joinpath(path, h5fname), "w") do h5

        # Write nodes
        h5["nodes"] = nodes

        # Write fields
        for (key, val) in fields

            # Error case: Unknown entry type
            if val["entry_type"] == "node" || val["entry_type"] == "cell"
                nothing
            elseif val["entry_type"] == "system"
                error("Entry type 'system' not implemented yet!")
            else
                error("Unkown entry type $(entry_type)")
            end

            # Format field in XDMF
            name = key
            center = titlecase(val["entry_type"])
            ftype = titlecase(val["field_type"])
            data = val["field_data"]

            # NOTE: For some reason dset[:, :] .= data doesn't write anything
            #    when data is a Base.Generator, so here I'm restricting it
            #    to only receive explicit arrays and matrices
            if !isa(data, Array)
                error("Data $(typeof(data)) not implemented yet.")
            end

            if ftype=="Scalar"

                fdims = length(data)

                if center=="Node" && length(data)!=nnodes
                    error("Field $(name) expected to have length $(nnodes); got length $(fdims).")
                end

                fdims = (1, fdims)

            elseif ftype=="Vector"

                fdims = size(data)

                if length(fdims)==1       # Case: Array of arrays

                    if length(first(data)) != 3
                        error("Expected array of array.")
                    end

                    # Reformat as matrix
                    # data = (V[i] for i in 1:3, V in data)
                    # NOTE: For some reason, dset .= data works only on arrays
                    data = [V[i] for i in 1:3, V in data]
                    fdims = size(data)

                elseif length(fdims)==2   # Case: Matrix

                    if center=="Node" && fdims != (3, nnodes)
                        error("Field $(name) expected to have size $((3, nnodes)); got $(fdims).")
                    elseif center=="Cell" && fdims != (3, ncells)
                        error("Field $(name) expected to have size $((3, ncells)); got $(fdims).")
                    end

                else
                    error("Field $(name) expected to have 1 or 2 dimensions; got $(length(fdims)).")
                end


            else
                error("Unknown field type $(val["field_type"])")
            end

            # Save field as HDF5 binary
            dtype = HDF5.datatype(Float64)
            dspace = HDF5.dataspace(fdims)
            dset = HDF5.create_dataset(h5, name, dtype, dspace; chunk=fdims, compression_optargs...)

            if debug
                println(key)
                println("\tfdims\t\t= $(fdims)")
                println("\tsize(data)\t= $(size(data))")
                println("\tdspace\t\t= $(dspace)")
                println("\tdset\t\t= $(dset)")
                println("\tsize(dset)\t= $(size(dset))")
            end

            write(dset, data)

            # Save XDMF properties
            xdmf_fields[name] = (center, ftype, fdims, data)

        end
    end

    # ------------ Write XDMF reader file -----------------------
    open(joinpath(path, xmffname), "w") do xmf

        # Open xmf block
        print(xmf, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
        print(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n")
            print(xmf, "\t<Domain>\n")
              print(xmf, "\t\t<Grid Name=\"structuredgrid\" GridType=\"Uniform\">\n")

                  if time != nothing
                      print(xmf, "\t\t\t\t<Time Value=\"", time, "\"/>\n")
                  end

                  # Geometry: Node position
                  print(xmf, "\t\t\t\t<Geometry Type=\"XYZ\">\n")
                    print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                                " Dimensions=\"", 3, " ", nnodes,
                                "\" Format=\"HDF\" Precision=\"8\">",
                                h5fname, ":nodes</DataItem>\n")
                  print(xmf, "\t\t\t\t</Geometry>\n")

                  # Topology: Declare connectivity as structured grid
                  print(xmf, "\t\t\t\t<Topology Type=\"3DSMesh\" Dimensions=\"", dims_str, "\"/>\n")

                  # Attribute: Fields
                  for (name, (center, ftype, fdims, data)) in xdmf_fields

                    # fdims_str = prod("$(n)"*" "^(n!=fdims[end]) for (ni, n) in enumerate(fdims))
                    fdims_str = prod("$(n - 1*(center=="Cell")) " for n in dims)
                    fdims_str *= "$(fdims[1])"


                    print(xmf, "\t\t\t\t<Attribute Name=\"", name,
                                          "\" Center=\"", center, "\" Type=\"", ftype, "\">\n")
                      print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                                    " Dimensions=\"", fdims_str, "\" Format=\"HDF\" Precision=\"8\">",
                                    h5fname, ":", name, "</DataItem>\n")
                    print(xmf, "\t\t\t\t</Attribute>\n")

                    # TODO: Int DataType?
                  end

              print(xmf, "\t\t</Grid>\n")
            print(xmf, "\t</Domain>\n")
        print(xmf, "</Xdmf>\n")
    end

    return xmffname*";"
end
