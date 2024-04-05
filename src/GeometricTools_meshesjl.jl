#=##############################################################################
# DESCRIPTION
    Interface to Meshes.jl

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Mar 2024
  * License   : MIT License
=###############################################################################


# function get_node(mesh::Meshes.SimpleMesh, coor::Vector{Int})
#
#     _check_simplemesh_coor(coor)
#
#     return get_node(mesh, coor[1])
# end
#
# function get_node(mesh::Meshes.SimpleMesh, i::Int)
#
#     if i>length(mesh.vertices)
#         error("Requested invalid node index $i; max is $(length(mesh.vertices)).")
#     elseif i<1
#         error("Invalid index $i (it must be greater than 0).")
#     end
#
#     return mesh.vertices[i].coords
# end



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
function get_cell_t!(out, cell::NTuple{N, T}) where {N, T}
    @assert length(out) == N
        "Received storage array of length $(length(out)); expected length $N"

    out .= cell

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

"""
Converts the vertices of a Meshes.jl object (an array of m n-th dimensional
Points) to a collection of nodes (an nxm Matrix).
"""
function vertices2nodes(vertices::AbstractVector{P}
                                    ) where {N, R, P<:Meshes.Primitive{N, R}}

    nodes = zeros(R, N, length(vertices))
    vertices2nodes!(vertices, nodes)

    return nodes
end

function vertices2nodes!(vertices::AbstractVector{P}, nodes::AbstractMatrix
                                                ) where {P<:Meshes.Primitive}
    for (point, node) in zip(vertices, eachcol(nodes))
        node .= point.coords
    end
end

function vertices2nodes(vertices::AbstractVector{P},
                        topology::Meshes.SimpleTopology{Meshes.Connectivity{PL, 2}}
                        ) where {N, R, P<:Meshes.Primitive{N, R}, PL}

    nodes = zeros(R, N, length(vertices))

    ni = 1
    previndices = (-1, -1)
    for (si, elem) in enumerate(topology.elems)

        segment = topology.connec[si]
        # segment = topology.connec[elem]
        indices = segment.indices

        # NOTE: Here we will assume that each node is only used for at most two
        #       segments
        # NOTE: We also assume that the mesh is continuous, meaning that cells
        #       are already contiguous

        if !(indices[2] in previndices) && !(indices[1] in previndices)
            println("$(si)\t$(ni)\tprevindices: $(previndices)\tindices: $(indices)\tnextindices: $(topology.connec[si+1].indices)")
        end

        # Add the first point if it wasn't shared by the previous segment
        if !(indices[1] in previndices)

            nodes[:, ni] .= vertices[indices[1]].coords
            ni += 1
        end

        # Add the second point if it wasn't shared by the previous segment
        if !(indices[2] in previndices)
            nodes[:, ni] .= vertices[indices[2]].coords
            ni += 1
        end

        previndices = indices

    end

    # Remove any extra memory that was not used
    if ni != size(nodes, 2)
        nodes = nodes[:, 1:ni]
    end

    return nodes
end

function vertices2nodes(mesh::Meshes.SimpleMesh)
    return vertices2nodes(mesh.vertices, mesh.topology)
end
