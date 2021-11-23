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

get_ndivsnodes(self::GridExtentions) = self._ndivsnodes
get_ndivscells(self::GridExtentions) = self._ndivscells

function get_cellnodes(self::GridExtentions, i::Int64)
  return [get_node(self, n) for n in get_cell(self, i)]
end
function get_cellnodes(self::GridExtentions, coor::Array{Int64,1})
  return get_cellnodes(self, sub2ind(self._ndivscells, coor...))
end

"""
  `get_fieldval(grid, field_name, coor)`

Returns the value of node of coordinates `coor` (1-indexed) in the field
'field_name'.
"""
function get_fieldval(self::GridExtentions, field_name::String, coor::Array{Int64,1})

  if self.field[field_name]["entry_type"]=="node"
    return get_fieldval(self, field_name, Base._sub2ind(self._ndivsnodes, coor...)) #TODO: use LinearIndex instead of _sub2ind() (next line too)

  elseif self.field[field_name]["entry_type"]=="cell"
    return get_fieldval(self, field_name, Base._sub2ind(self._ndivscells, coor...))

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
                            field_type::String, field_data, entry_type::String;
                            raise_warn::Bool=true)
  # Error cases
  if !(field_type in ["vector", "scalar"])
    error("Unkown field type $(field_type)")
  elseif !(entry_type in ["node", "cell", "system"])
    error("Unkown entry type $(entry_type)")
elseif entry_type=="node" && length(field_data)!=grid.nnodes
    error("Invalid node field size."*
            " Expected $(grid.nnodes), got $(length(field_data))")
  elseif entry_type=="cell" && length(field_data)!=grid.ncells
    error("Invalid cell field size."*
            " Expected $(grid.nnodes), got $(length(field_data))")
  end


  if field_name in keys(grid.field) && raise_warn
    @warn("Overwritting field $field_name.")
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
  `save(grid, filename; args...)`

Outputs a vtk file of this grid. See generateVTK for a descrition of optional
arguments `args...`.
"""
function save_vtk(grid::GridExtentions, filename::String; O=nothing, Oaxis=nothing,
                                                                        args...)
  # Determines whether to add 0 to points for vtk file
  if grid.dims<=3
    aux1 = zeros(Float64, 3-grid.dims)
  else
    error("$(grid.dims)-dimensional grids can't be exported as VTKs!")
  end

  if O != nothing || Oaxis != nothing

    if grid.dims != 3
      error("Requested space transformation on grid of $(grid.dims) "*
            "dimensions. There is not such implementation yet!")
    end
    O = O==nothing ? zeros(3) : O
    Oaxis = Oaxis==nothing ? Array(1.0I, 3, 3) : Oaxis
    check_coord_sys(Oaxis)
    invOaxis = collect(Oaxis')

    points = [countertransform(get_node(grid, i), invOaxis, O) for i in 1:grid.nnodes]

  else

    points = [vcat(get_node(grid, i), aux1) for i in 1:grid.nnodes]

  end

  cells = [get_cell(grid, i) .- 1 for i in 1:grid.ncells]

  point_data = []
  cell_data = []
  for (key, val) in grid.field
    if val["entry_type"] == "node"
      push!(point_data, val)
    elseif val["entry_type"] == "cell"
      push!(cell_data, val)
    elseif val["entry_type"] == "system"
      nothing
    else
      error("Unkown entry type $(entry_type)")
    end
  end

  # Checks for quasi-dimensions (NDIVS[i]==0)
  qdims = length( [1 for ndiv in grid._ndivscells if ndiv==0] )
  dims = grid.dims - qdims

  ctype = grid._override_vtkcelltype

  return generateVTK(filename, points; cells=cells,
                      point_data=size(point_data,1)!=0 ? point_data : nothing,
                      cell_data=size(cell_data,1)!=0 ? cell_data : nothing,
                      _griddims=dims,
                      override_cell_type=ctype,
                      args...)
end

function save_xdmf(grid::GridExtentions, filename; optargs...)

    # NOTE: loop_dims is currently not supported and this will crash!

    return generateXDMF_3Dstructured(filename, grid.nodes, grid._ndivscells;
                                        fields=grid.field,
                                        optargs...)
end

function save(grid::GridExtentions, filename; format="xdmf", optargs...)

    frmt = lowercase(format)

    if frmt=="xdmf"
        return save_xdmf(grid, filename; optargs...)
    elseif frmt=="vtk"
        return save_vtk(grid, filename; optargs...)
    else
        error("Unknown format $(format). Valid formats are \"xdmf\" and \"vtk\".")
    end
end
