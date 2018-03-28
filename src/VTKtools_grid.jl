#=##############################################################################
# DESCRIPTION
    Methods for generation and evaluation of grids.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Mar 2018
  * License   : MIT License
=###############################################################################

################################################################################
# GRID
################################################################################
"""
  `Grid(P_min, P_max, NDIVS)`

Generates an n-dimensional grid.

  **Arguments**
  * `P_min::Array{Float64,1}`   : Minimum point of the domain.
  * `P_max::Array{Float64,1}`   : Maximum point of the domain.
  * `NDIVS::Array{Int64,1}`     : Number of divisions in each coordinate.

  **Properties**
  * `dims::Int64`               : Number of dimensions.
  * `nnodes::Int64`             : Number of nodes in the grid.
  * `nodes::Array{Float64,2}`   : Matrix size (`nnodes`, 3) of node position.
  * `field` : Contains calculated fields formated as field[field_name] = Dict(
                            "field_name" => field_name::String,
                            "field_type" => "scalar" or "vector",
                            "field_data" => data
                            )
            where `data` is an array data[i] = [val1, val2, ...] containing
            this field values (scalar or vector) at each node in the grid.

NOTE: All indexing is done linearly, meaning that `nodes` is indexed from 1 to
      `nnodes`, and all data fields follow the same indexing.

NOTE2: `NDIVS` can either be an array of integers with NDIVS[i] indicating the
      number of divisions in the i-th coordinate, or it can be an array of
      sections (see `multidiscretize()` doc) with NDIVS[i] = [sec1, sec2, ...]
      indicating the discretization into sections in the i-th coordinate.
"""
type Grid

  # User inputs
  P_min::Array{T,1} where {T<:Real}   # Minimum point of the domain
  P_max::Array{T,1} where {T<:Real}   # Maximum point of the domain
  NDIVS::Array{T,1} where {T<:Any}    # Number of divisions in each coordinate

  # Properties
  dims::Int64                         # Number of dimensions
  nnodes::Int64                       # Number of nodes
  ncells::Int64                       # Number of cells
  nodes::Array{T,2} where{T<:Real}    # Position of each node
  field::Dict{String, Dict{String, Any}}  # Calculated fields

  # Internal data
  _ndivsnodes::Tuple                  # Number of nodes in each coordinate
  _ndivscells::Tuple                  # Number of cells in each coordinate

  Grid(P_min, P_max, NDIVS,
              dims=_calc_dims(P_min),
                nnodes=_calc_nnodes(NDIVS),
                ncells=_calc_ncells(NDIVS),
                nodes=_generate_grid(P_min, P_max, NDIVS),
                field=Dict{String, Dict{String, Any}}(),
              _ndivsnodes=Tuple(_calc_ndivs(NDIVS)+1),
                _ndivscells=Tuple(_calc_ndivs(NDIVS))
      ) = _check(P_min, P_max, NDIVS) ? new(P_min, P_max, NDIVS,
              dims,
                nnodes,
                ncells,
                nodes,
                field,
              _ndivsnodes,
                _ndivscells
      ) : error("Logic error!")
end


"""
  `get_node(grid, i)`

Returns the position of the i-th node (1-indexed) in the grid
"""
function get_node(self::Grid, i::Int64)
  if i>self.nnodes
    error("Requested invalid node index $i; max is $(self.nnodes).")
  end

  return self.nodes[i, :]
end

"""
  `get_node(grid, coor)`

Returns the position of the node of subscript coordinates `coor` (1-indexed)
"""
function get_node(self::Grid, coor::Array{Int64,1})
  return get_node(self, sub2ind(self._ndivsnodes, coor...))
end

"""
  `get_cell(grid, i)`

Returns the nodes indices of i-th cell in the grid (1-indexed)
"""
function get_cell(self::Grid, i::Int64)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  end
  return get_cell(self, collect(ind2sub(self._ndivscells, i)))
end

"""
  `get_cell(grid, coor)`

Returns the node indices of the cell with subscript coordinates `coor`
(1-indexed). The format corresponds to VTK_HEXAHEDRON (=12) in 3D, VTK_QUAD (=9)
in 2D, or VTK_LINE (=3) in 1D---except that points are 1-indexed instead of
0-indexed.
"""
function get_cell(self::Grid, coor::Array{Int64,1})

  if self.dims==1
    return [sub2ind(self._ndivsnodes, (coor+[0])...),
            sub2ind(self._ndivsnodes, (coor+[1])...)]

  elseif self.dims==2
    return [sub2ind(self._ndivsnodes, (coor+[0,0])...),
            sub2ind(self._ndivsnodes, (coor+[1,0])...),
            sub2ind(self._ndivsnodes, (coor+[1,1])...),
            sub2ind(self._ndivsnodes, (coor+[0,1])...)]

  elseif self.dims==3
    return [sub2ind(self._ndivsnodes, (coor+[0,0,0])...),
            sub2ind(self._ndivsnodes, (coor+[1,0,0])...),
            sub2ind(self._ndivsnodes, (coor+[1,1,0])...),
            sub2ind(self._ndivsnodes, (coor+[0,1,0])...),
            sub2ind(self._ndivsnodes, (coor+[0,0,1])...),
            sub2ind(self._ndivsnodes, (coor+[1,0,1])...),
            sub2ind(self._ndivsnodes, (coor+[1,1,1])...),
            sub2ind(self._ndivsnodes, (coor+[0,1,1])...)]

  else
    error("Definition of $(self.ndims)-dimensional cells not implemented yet!")
  end
end

"""
  `add_field(grid::Grid, field_name::String, field_type::String, field_data)`

Adds a field of data associated to each node.

NOTE: each data entry must be a single value if `field_type==scalar`, or a
      3-element array if `field_type==vector`.
"""
function add_field(grid::Grid, field_name::String, field_type::String,
                                                                field_data)
  # Error cases
  if !(field_type in ["vector", "scalar"])
    error("Unkown field type $(field_type)")
  elseif size(field_data, 1)!=grid.nnodes
    error("Invalid field size."*
            " Expected $(grid.nnodes), got $(size(field_data, 1))")
  end

  if field_name in keys(grid.field)
    warn("Overwritting field $field_name.")
  end

  grid.field[field_name] = Dict("field_name" => field_name,
                                "field_type" => field_type,
                                "field_data" => field_data)
end

"""
  `calculate_field(grid::Grid, f, field_name::String, field_type::String)`

Evaluates the function `f` at each nodes and stores the values as a new field.

NOTE: f must return a single value if `field_type==scalar`, or a 3-element array
      if `field_type==vector`.
"""
function calculate_field(grid::Grid, f, field_name::String, field_type::String)
  field_data = [f(get_node(grid, i)) for i in 1:grid.nnodes]
  add_field(grid, field_name, field_type, field_data)
end

"""
  `lintransform!(grid::Grid, M::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates the grid by the rotation matrix `M` and translation
vector `T` (linear transformation).
"""
function lintransform!(grid::Grid, M::Array{Float64,2}, T::Array{Float64,1};
                        reset_fields::Bool=true)
  # Error cases
  if size(M, 1)!=grid.dims
    error("Invalid rotation matrix `M`."*
            "Expected $(grid.dims) dimensions, got $(size(M,1))")
  elseif length(T)!=grid.dims
    error("Invalid translation vector `T`."*
            "Expected $(grid.dims) dimensions, got $(length(T))")
  end

  if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

  invM = inv(M)
  for i in 1:grid.nnodes
    grid.nodes[i,:] = countertransform(grid.nodes[i,:], invM, T)
  end
end

"""
  `transform!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid.
"""
function transform!(grid::Grid, f; reset_fields::Bool=true)

    if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

    for i in 1:grid.nnodes
      grid.nodes[i,:] = f(grid.nodes[i,:])
    end
end

"""
  `save(grid, filename; args...)`

Outputs a vtk file of this grid
"""
function save(grid::Grid, filename::String; args...)
  # Determins whether to add 0 to points for vtk file
  if grid.dims<=3
    aux1 = zeros(3-grid.dims)
  else
    error("$(grid.dims)-dimensional grids can't be exported as VTKs!")
  end

  points = [vcat(get_node(grid, i), aux1) for i in 1:grid.nnodes]
  cells = [get_cell(grid, i)-1 for i in 1:grid.ncells]
  point_data = length(grid.field)==0 ? nothing : values(grid.field)

  generateVTK(filename, points; cells=cells, point_data=point_data,
                                                  _griddims=grid.dims, args...)
end

##### INTERNAL FUNCTIONS  ######################################################
function _check(P_min::Array{T,1} where {T<:Real},
                P_max::Array{T,1} where {T<:Real} ,
                NDIVS::Array{T,1} where {T<:Any})
  # Error cases
  if size(P_min)!=size(P_max)
    error("`P_min` and `P_max` must have the same dimensions "*
                                          "($(size(P_min))!=$(size(P_max)))")
  elseif typeof(NDIVS)==Array{Int64,1}
    if _calc_dims(P_min)!=length(NDIVS)
      error("Division for each dimension must be given "*
                                          "$(length(P_min))!=$(length(NDIVS))")
    end
    for this_div in NDIVS
      if this_div<=0
        error("Invalid division $this_div")
      end
    end
  end
  for i in 1:_calc_dims(P_min)
    if P_min[i]>=P_max[i]
      error("Values in `P_min` must be less than values in `P_max` "*"
                                          ($(P_min[i])>=$(P_max[i]))")
    end
  end

  return true
end

function _calc_ndivs(NDIVS::Array{T,1} where {T<:Any})
  if typeof(NDIVS)==Array{Int64,1}
    return NDIVS

  elseif typeof(NDIVS)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    return [ sum([sec[2] for sec in secs]) for secs in NDIVS ]

  else
      error("`NDIVS` expected to be type `Array{Int64,1}` or"*
              " `Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}`; got "*
              "$(typeof(NDIVS))")
  end

end


"Calculates number of nodes given NDIVS"
function _calc_nnodes(NDIVS::Array{T,1} where {T<:Any})
  return prod(_calc_ndivs(NDIVS)+1)
end

"Calculates the number of cells given NDIVS"
function _calc_ncells(NDIVS::Array{T,1} where {T<:Any})
  return prod(_calc_ndivs(NDIVS))
end

function _calc_dims(P::Array{T,1} where {T<:Real})
  return length(P)
end

function _generate_grid(P_min::Array{T,1} where {T<:Real},
                        P_max::Array{T,1} where {T<:Real},
                        NDIVS::Array{T,1} where {T<:Any})

  dims = _calc_dims(P_min)
  nnodes = _calc_nnodes(NDIVS)
  nodes = zeros(nnodes, dims)
  ndivs = Tuple(_calc_ndivs(NDIVS))
  ndivsnodes = Tuple(_calc_ndivs(NDIVS)+1)

  # Discretizes each coordinate according to NDIVS
  if typeof(NDIVS)==Array{Int64,1}
    spacing = [linspace(P_min[i], P_max[i], ndivs[i]+1) for i in 1:dims]

  elseif typeof(NDIVS)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    spacing = [multidiscretize(x->x, P_min[i], P_max[i], NDIVS[i]) for i in 1:dims]

  else
    error("Not implemented yet")
  end

  # Creats the grid
  for ind in 1:nnodes
    sub = ind2sub(ndivsnodes, ind)
    p = [spacing[dim][sub[dim]] for dim in 1:dims]
    nodes[ind, :] = p
  end

  return nodes
end






##### END OF GRID ##############################################################
