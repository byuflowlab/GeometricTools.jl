struct Field{F}
    field_name::String
    field_type::String
    entry_type::String
    field_data::Vector{F}
end

struct Grid{N,TF}# <: AbstractGrid
    ndivs::NTuple{N, Int}                          # Number of divisions in each dimension
    nodes::Matrix{TF}                              # Nodes stored in each column
    loop_dim::NTuple{N, Bool}                      # Indicates whether first and last point in a given direction is identical
    scalar_fields::Dict{String, Field{TF}}         # Scalar fields
    vector_fields::Dict{String, Field{Vector{TF}}} # Vector fields
end

"""
`Grid(f, spacing, loop_dim)`
"""
function Grid{N,TF}(f, spacing::NTuple{N, <:AbstractVector}, loop_dim::NTuple{N, Bool}) where {N,TF}

    # get number of divisions
    ndivs = length.(spacing) .- 1
    nnodes = get_nnodes(ndivs, loop_dim)
    cidx = CartesianIndices(nnodes)
    total_nnodes = prod(nnodes)

    # initialize nodes
    nodes = zeros(TF, N, total_nnodes)

    # apply transformation to nodes
    for i = 1:total_nnodes
       point = SVector{N}([spacing[j][cidx[i][j]] for j = 1:N])
       nodes[:, i] .= f(point)
    end

    # initialize empty scalar and vector fields
    scalar_fields = Dict{String, Field{TF}}()
    vector_fields = Dict{String, Field{Vector{TF}}}()

    return Grid{N,TF}(ndivs, nodes, loop_dim, scalar_fields, vector_fields)
end

Grid{N,TF}(spacing::NTuple{N, <:AbstractVector}, loop_dim::NTuple{N, Bool}) where {N,TF} =
    Grid{N,TF}((x)->x, spacing, loop_dim)

# --- Required Functions --- #

"""
  `get_node(grid::Grid, i::Integer)`

Returns the position of the i-th node (1-indexed) in the grid
"""
get_node(grid::Grid, i::Integer) = grid.nodes[:, i]

"""
  `get_node(grid::Grid, coor::Vararg{<:Integer, N})`

Returns the position of the node of subscript coordinates `coor` (1-indexed)
"""
get_node(grid::Grid, I::Vararg{<:Integer, N}) where N = get_node(grid, LinearIndices(grid.ndivs)[I...])

"""
  `get_cell(grid, i)`

Returns the nodes indices of i-th cell in the grid (1-indexed)
"""
get_cell(grid::Grid, i::Integer) = get_cell(grid, Tuple(CartesianIndices(dims)[i])...)

"""
  `get_cell(grid, I)`

Returns the nodes indices of i-th cell in the grid (1-indexed)
"""
function get_cell(grid::Grid{N,TF}, I::Vararg{<:Integer, N}) where {N,TF}

    d = ones(Int, N)
    for i = 1:N
        if grid.loop_dim[i] && I[i] == grid.ndivs[i]
            d[i] = 1 - grid.ndivs[i]
        end
    end

    dims = count(grid.ndivs .> 0)

    if dims == 1
        return [LinearIndices(self.nnodes)[I],
                LinearIndices(self.nnodes)[I .+ (d1, 0, 0)[1:N]]]
    elseif dims == 2
        return [LinearIndices(self.nnodes)[I .+ ( 0,  0, 0)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1,  0, 0)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1, d2, 0)[1:N]],
                LinearIndices(self.nnodes)[I .+ ( 0, d2, 0)[1:N]]]
    elseif dims == 3
        return [LinearIndices(self.nnodes)[I .+ ( 0,  0,  0)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1,  0,  0)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1, d2,  0)[1:N]],
                LinearIndices(self.nnodes)[I .+ ( 0, d2,  0)[1:N]],
                LinearIndices(self.nnodes)[I .+ ( 0,  0, d3)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1,  0, d3)[1:N]],
                LinearIndices(self.nnodes)[I .+ (d1, d2, d3)[1:N]],
                LinearIndices(self.nnodes)[I .+ ( 0, d2, d3)[1:N]]]
    else
        error("Definition of $(N)-dimensional cells not implemented yet!")
    end
end

get_scalar_field(grid::Grid{N,<:Any}, field_name::String, i::Integer) where N =
    grid.scalar_field[field_name]["field_data"][i]

"""
`add_scalar_field(grid, field_name, I)`
"""
function get_scalar_field(grid::Grid{N,<:Any}, field_name::String, I::Vararg{<:Integer, N}) where N
    if self.field[field_name]["entry_type"]=="node"
        i = LinearIndices(get_nnodes(grid.ndivs))[I...]
    elseif self.field[field_name]["entry_type"] == "cell"
        i = LinearIndices(grid.ndivs)[I...]
    else
        error("Entry type $(self.field[field_name]["entry_type"]) not implemented.")
    end
    return get_scalar_field(grid, field_name, i)
end

get_vector_field(grid::Grid{N,<:Any}, field_name::String, i::Integer) where N =
    grid.vector_field[field_name]["field_data"][i]

"""
`add_vector_field(grid, field_name, I)`
"""
function get_vector_field(grid::Grid{N,<:Any}, field_name::String, I::Vararg{<:Integer, N}) where N
    if self.field[field_name]["entry_type"]=="node"
        i = LinearIndices(get_nnodes(grid.ndivs))[I...]
    elseif self.field[field_name]["entry_type"] == "cell"
        i = LinearIndices(grid.ndivs)[I...]
    else
        error("Entry type $(self.field[field_name]["entry_type"]) not implemented.")
    end
    return get_vector_field(grid, field_name, i)
end

"""
`add_scalar_field!(grid, field_name, field_data, entry_type; warn_on_overwrite=true)`
"""
function add_scalar_field!(grid::Grid, field_name, field_data, entry_type; warn_on_overwrite=true)

    check_entry_type(grid, field_data, entry_type)

    if field_name in keys(grid.scalar_fields) && warn_on_overwrite
        @warn("Overwritting field $field_name.")
    end

    grid.scalar_field[field_name] = Field(field_name, "scalar", entry_type, field_data)

    return grid
end

"""
`add_vector_field!(grid, field_name, field_data, entry_type; warn_on_overwrite=true)`
"""
function add_vector_field!(grid::Grid, field_name, field_data, entry_type; warn_on_overwrite=true)

    check_entry_type(grid, field_data, entry_type)

    if field_name in keys(grid.scalar_fields) && warn_on_overwrite
        @warn("Overwritting field $field_name.")
    end

    grid.vector_field[field_name] = Field(field_name, "vector", entry_type, field_data)

    return grid
end

"""
`check_entry_type(grid, field_data, entry_type)`
"""
function check_entry_type(grid::Grid, field_data, entry_type)
    if !(entry_type in ["node", "cell", "system"])
        error("Unkown entry type $(entry_type)")
    end

    if entry_type=="node" && size(field_data, 1) != prod(get_nnodes(grid.ndivs))
        error("Invalid node field size."*
                " Expected $(grid.nnodes), got $(size(field_data, 1))")
    end

    if entry_type=="cell" && size(field_data, 1) != get_ncells(grid.ndivs)
        error("Invalid cell field size."*
                " Expected $(grid.nnodes), got $(size(field_data, 1))")
    end

    return nothing
end

"""
`add_scalar_field!(f, grid, field_nname, entry_type; warn_on_overwrite=true)`
"""
function add_scalar_field!(f, grid::Grid, field_name, entry_type; warn_on_overwrite=true)

    if entry_type=="node"
        field_data = [f(get_node(grid, i)) for i in 1:grid.nnodes]
    elseif entry_type=="cell"
        field_data = [f(get_cell(grid, i)) for i in 1:grid.ncells]
    else
        error("Unkown entry type $(entry_type)")
    end

    return add_scalar_field!(grid, field_name, field_data, entry_type, warn_on_overwrite=warn_on_overwrite)
end

"""
`add_vector_field!(f, grid, field_name, entry_type; warn_on_overwrite=true)`
"""
function add_vector_field!(f, grid::Grid, field_name, entry_type; warn_on_overwrite=true)

    if entry_type=="node"
        field_data = [f(get_node(grid, i)) for i in 1:grid.nnodes]
    elseif entry_type=="cell"
        field_data = [f(get_cell(grid, i)) for i in 1:grid.ncells]
    else
        error("Unkown entry type $(entry_type)")
    end

    return add_vector_field!(grid, field_name, field_data, entry_type, warn_on_overwrite=warn_on_overwrite)
end

"""
  `save(grid, filename; args...)`

Outputs a vtk file of this grid. See generateVTK for a descrition of optional
arguments `args...`.
"""
function save(grid::Grid, filename, args...)
    return generateVTK()
end


"""
    get_nnodes(ndivs::NTuple{N, <:Integer}, loop_dim::NTuple{N, Bool})

Calculates number of nodes in each dimension
"""
get_nnodes(ndivs::NTuple{N, <:Integer}, loop_dim::NTuple{N, Bool}) where N = ndivs .- loop_dim .+ 1

"""
    get_ncells(ndivs::NTuple{N, <:Integer})

Calculates number of cells
"""
get_ncells(ndivs::NTuple{N, <:Integer}) where N = prod(max.(ndivs, 1))






"""
  `transform!(grid::Grid, M, T)`

Rotates and translates the grid by the rotation matrix `M` and translation
vector `T` (linear transformation).
"""
function transform!(grid::Grid, M, T; reset_fields=true)

    if reset_fields
        for key in keys(grid.scalar_fields)
            delete!(grid.scalar_fields, key)
        end
        for key in keys(grid.vector_fields)
            delete!(grid.vector_fields, key)
        end
    end

    invM = inv(M)
    for i in 1:grid.nnodes
        grid.nodes[:,i] .= countertransform(grid.nodes[:,i], invM, T)
    end

    return grid
end

"""
  `transform!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
position of every node is given to the function `f`.
"""
function transform!(grid::Grid, f; reset_fields::Bool=true)

    if reset_fields
        for key in keys(grid.scalar_fields)
            delete!(grid.scalar_fields, key)
        end
        for key in keys(grid.vector_fields)
            delete!(grid.vector_fields, key)
        end
    end

    for i in 1:grid.nnodes
        grid.nodes[:,i] .= f(grid.nodes[:,i])
    end

    return grid
end


"""
  `transform2!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
indices of every node is given to the function `f`.
"""
function transform2!(grid::Grid, f; reset_fields=true)

    if reset_fields
        for key in keys(grid.scalar_fields)
            delete!(grid.scalar_fields, key)
        end
        for key in keys(grid.vector_fields)
            delete!(grid.vector_fields, key)
        end
    end

    for i in 1:grid.nnodes
        grid.nodes[:,i] .= f(CartesianIndices(grid.ndivs)[i])
    end
end

"""
  `transform3!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
both the position and indices of every node is given to the function `f`.
"""
function transform3!(grid::Grid, f; reset_fields=true)

    if reset_fields
        for key in keys(grid.scalar_fields)
            delete!(grid.scalar_fields, key)
        end
        for key in keys(grid.vector_fields)
            delete!(grid.vector_fields, key)
        end
    end

    for i in 1:grid.nnodes
        grid.nodes[:,i] .= f(grid.nodes[:,i], CartesianIndices(grid.ndivs)[i])
    end
end
