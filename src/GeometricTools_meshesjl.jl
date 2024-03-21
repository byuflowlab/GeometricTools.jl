#=##############################################################################
# DESCRIPTION
    Interface to Meshes.jl

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Mar 2024
  * License   : MIT License
=###############################################################################

function get_node(mesh::Meshes.SimpleMesh, i::Int)

    if i>length(mesh.vertices)
        error("Requested invalid node index $i; max is $(length(mesh.vertices)).")
    elseif i<1
        error("Invalid index $i (it must be greater than 0).")
    end

    return mesh.vertices[i].coords
end

function get_node(mesh::Meshes.SimpleMesh, coor::Vector{Int})

    for (coori, i) in enumerate(coor)
        if coori == 1
            nothing
        else
            @assert i==0 || i==1 ""*
                "Requested coordinates $(coor) on an unstructured grid but"*
                " only first coordinate can be different than 0 or 1."
        end
    end

    return get_node(mesh, coor[1])
end

function get_cell(mesh::Meshes.SimpleMesh, i::Int)

    if i>length(mesh)
        error("Requested invalid cell index $i; max is $(length(mesh)).")
    elseif i<1
        error("Invalid index $i (it must be greater than 0).")
    end

    return get_cell(mesh.topology, i)
end

get_cell(topology::Meshes.SimpleTopology, i::Int) = topology.connec[i]
# get_cell(topology::Mesh.HalfEdgeTopology, i::Int) =

function get_cell(mesh::Meshes.SimpleMesh, coor::Vector{Int})

    for (coori, i) in enumerate(coor)
        if coori == 1
            nothing
        else
            @assert i==0 || i==1 ""*
                "Requested coordinates $(coor) on an unstructured grid but"*
                " only first coordinate can be different than 0 or 1."
        end
    end

    return get_cell(mesh, coor[1])
end

function get_cell_t!(out, mesh::Meshes.SimpleMesh, coor::AbstractVector{Int})

    for (coori, i) in enumerate(coor)
        if coori == 1
            nothing
        else
            @assert i==0 || i==1 ""*
                "Requested coordinates $(coor) on an unstructured grid but"*
                " only first coordinate can be different than 0 or 1."
        end
    end

    return get_cell_t!(out, mesh, coor[1])
end

# Here we call the not-inplace version since it is non-allocating
get_cell_t!(out, mesh, i::Int) = get_cell_t!(out, get_cell(mesh, i))

function get_cell_t!(out, cell::Meshes.Connectivity{PL, N}) where {PL, N}
    @assert length(out) == N
        "Received storage array of length $(length(out)); expected length $N"

    out .= cell.indices

    return out
end

# function get_cell_t!(out, args1..., mesh::Meshes.SimpleMesh, coor::AbstractVector{Int}, args2...; optargs...)
#     return get_cell_t!(out, mesh, coor)
# end
