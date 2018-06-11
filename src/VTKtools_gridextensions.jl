#=##############################################################################
# DESCRIPTION
    Functions common to all grids that have the properties `_ndivsnodes`,
    `_ndivsnodes`, and `_override_vtkcelltype`.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : MIT License
=###############################################################################

"""
  `get_fieldval(grid, field_name, coor)`

Returns the value of node of coordinates `coor` (1-indexed) in the field
'field_name'.
"""
function get_fieldval(self::GridExtentions, field_name::String, coor::Array{Int64,1})
  if !(field_name in keys(self.field))
    error("Field $field_name doesn't exist."*
          " Available fields: $(keys(self.field))")
  end

  if self.field[field_name]["entry_type"]=="node"
    return get_fieldval(self, field_name, sub2ind(self._ndivsnodes, coor...))

  elseif self.field[field_name]["entry_type"]=="cell"
    return get_fieldval(self, field_name, sub2ind(self._ndivscells, coor...))

  else
    error("Entry type $(self.field[field_name]["entry_type"]) not implemented.")
  end
end


"""
  `get_fieldval(grid, field_name, i)`

Returns the value of i-th node (1-indexed) in the field 'field_name'.
"""
function get_fieldval(self::GridExtentions, field_name::String, i::Int64)
  if !(field_name in keys(self.field))
    error("Field $field_name doesn't exist."*
          " Available fields: $(keys(self.field))")
  end

  return self.field[field_name]["field_data"][i]
end

"""
  `add_field(grid::Grid, field_name::String, field_type::String, field_data,
entry_type::String)`

Adds a field of data associated to each node.

NOTE: each data entry must be a single value if `field_type==scalar`, or a
      3-element array if `field_type==vector`.
"""
function add_field(grid::GridExtentions, field_name::String,
                            field_type::String, field_data, entry_type::String)
  # Error cases
  if !(field_type in ["vector", "scalar"])
    error("Unkown field type $(field_type)")
  elseif !(entry_type in ["node", "cell"])
    error("Unkown entry type $(entry_type)")
  elseif entry_type=="node" && size(field_data, 1)!=grid.nnodes
    error("Invalid node field size."*
            " Expected $(grid.nnodes), got $(size(field_data, 1))")
  elseif entry_type=="cell" && size(field_data, 1)!=grid.ncells
    error("Invalid cell field size."*
            " Expected $(grid.nnodes), got $(size(field_data, 1))")
  end


  if field_name in keys(grid.field)
    warn("Overwritting field $field_name.")
  end

  grid.field[field_name] = Dict("field_name" => field_name,
                                "field_type" => field_type,
                                "entry_type" => entry_type,
                                "field_data" => field_data)
end

"""
  `calculate_field(grid::Grid, f, field_name::String, field_type::String)`

Evaluates the function `f` at each nodes and stores the values as a new field.

NOTE: f must return a single value if `field_type==scalar`, or a 3-element array
      if `field_type==vector`.
"""
function calculate_field(grid::GridExtentions, f, field_name::String,
                                        field_type::String, entry_type::String)
  if entry_type=="node"
    field_data = [f(get_node(grid, i)) for i in 1:grid.nnodes]
  elseif entry_type=="cell"
    field_data = [f(get_cell(grid, i)) for i in 1:grid.ncells]
  else
    error("Unkown entry type $(entry_type)")
  end
  add_field(grid, field_name, field_type, field_data, entry_type)
end

"""
  `lintransform!(grid::Grid, M::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates the grid by the rotation matrix `M` and translation
vector `T` (linear transformation).
"""
function lintransform!(grid::GridExtentions, M::Array{Float64,2}, T::Array{Float64,1};
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
function transform!(grid::GridExtentions, f; reset_fields::Bool=true)

    if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

    for i in 1:grid.nnodes
      grid.nodes[i,:] = f(grid.nodes[i,:])
    end
end

"""
  `save(grid, filename; args...)`

Outputs a vtk file of this grid
"""
function save(grid::GridExtentions, filename::String; args...)
  # Determins whether to add 0 to points for vtk file
  if grid.dims<=3
    aux1 = zeros(3-grid.dims)
  else
    error("$(grid.dims)-dimensional grids can't be exported as VTKs!")
  end

  points = [vcat(get_node(grid, i), aux1) for i in 1:grid.nnodes]
  cells = [get_cell(grid, i)-1 for i in 1:grid.ncells]

  point_data = []
  cell_data = []
  for (key, val) in grid.field
    if val["entry_type"] == "node"
      push!(point_data, val)
    elseif val["entry_type"] == "cell"
      push!(cell_data, val)
    else
      error("Unkown entry type $(entry_type)")
    end
  end

  # Checks for quasi-dimensions (NDIVS[i]==0)
  qdims = length( [1 for ndiv in grid._ndivscells if ndiv==0] )
  dims = grid.dims - qdims

  ctype = grid._override_vtkcelltype

  generateVTK(filename, points; cells=cells,
                      point_data=size(point_data,1)!=0 ? point_data : nothing,
                      cell_data=size(cell_data,1)!=0 ? cell_data : nothing,
                      _griddims=dims,
                      override_cell_type=ctype,
                      args...)
end