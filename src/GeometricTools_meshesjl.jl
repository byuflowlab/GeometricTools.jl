#=##############################################################################
# DESCRIPTION
    Interface to Meshes.jl

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Mar 2024
  * License   : MIT License
=###############################################################################


function get_node(mesh::Meshes.SimpleMesh, coor::Vector{Int})

    _check_simplemesh_coor(coor)

    return get_node(mesh, coor[1])
end

function get_node(mesh::Meshes.SimpleMesh, i::Int)

    if i>length(mesh.vertices)
        error("Requested invalid node index $i; max is $(length(mesh.vertices)).")
    elseif i<1
        error("Invalid index $i (it must be greater than 0).")
    end

    return mesh.vertices[i].coords
end



function get_cell(mesh::Meshes.SimpleMesh, coor::Vector{Int})

    _check_simplemesh_coor(coor)

    return get_cell(mesh, coor[1])
end

function get_cell(mesh::Meshes.SimpleMesh, i::Int)

    if i>length(mesh)
        error("Requested invalid cell index $i; max is $(length(mesh)).")
    elseif i<1
        error("Invalid index $i (it must be greater than 0).")
    end

    return get_cell(mesh.topology, i)
end

# Base case
get_cell(topology::Meshes.SimpleTopology, i::Int) = topology.connec[i].indices
# get_cell(topology::Mesh.HalfEdgeTopology, i::Int) =




function get_cell_t!(out, mesh::Meshes.SimpleMesh, coor, ::Any, args...; optargs...)
    return get_cell_t!(out, mesh, coor)
end

function get_cell_t!(out, mesh::Meshes.SimpleMesh, coor::AbstractVector{Int})

    _check_simplemesh_coor(coor)

    return get_cell_t!(out, mesh, coor[1])
end

# Here we call the not-inplace version (get_cell) since it is non-allocating
get_cell_t!(out, mesh, i::Int) = get_cell_t!(out, get_cell(mesh, i))

# Base case
function get_cell_t!(out, cell::Meshes.Connectivity{PL, N}) where {PL, N}
    @assert length(out) == N
        "Received storage array of length $(length(out)); expected length $N"

    out .= cell.indices

    return out
end


function get_cell_t(mesh::Meshes.SimpleMesh, coor, nodei, ::Any, args...; optargs...)
    return get_cell_t(mesh, coor, nodei)
end

function get_cell_t(mesh::Meshes.SimpleMesh, coor::AbstractVector{Int}, nodei::Int)

    _check_simplemesh_coor(coor)

    return get_cell_t(mesh, coor[1], nodei)
end

get_cell_t(mesh, i::Int, nodei::Int) = get_cell(mesh, i)[nodei]




"Check for the correct coordinate format for an unstructured grid"
function _check_simplemesh_coor(coor)

    for (coori, i) in enumerate(coor)
        if coori == 1
            nothing
        else
            @assert i==0 || i==1 ""*
                "Requested coordinates $(coor) on an unstructured grid but"*
                " only first coordinate can be different than 0 or 1."
        end
    end

    return nothing
end

# function Base.getindex(point::Meshes.Point3, i::Int)
function Base.getindex(point::Meshes.Primitive, i::Int)
    return point.coords[i]
end
