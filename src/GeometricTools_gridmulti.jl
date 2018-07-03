#=##############################################################################
# DESCRIPTION
    Methods for generation and evaluation of multigrids.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : May 2018
  * License   : MIT License
=###############################################################################

################################################################################
# MULTIGRID
################################################################################

type MultiGrid <: AbstractGrid

  # User inputs
  dims::Int64                         # Number of dimensions

  # Properties
  grids::Array{T,1} where {T<:AbstractGrid}     # Grids
  grid_names::Array{String, 1}        # Grid names
  nnodes::Int64                       # Number of nodes
  ncells::Int64                       # Number of cells
  ngrids::Int64                       # Number of grids
  # bbox::Array{Int64, 1}               # Bounding box (cells in each dimension)

  # Internal data
  _grid_origins::Array{Array{Int64,1},1}  # Origin of each grid

  MultiGrid(          dims,
                      grids=AbstractGrid[],
                        grid_names=String[],
                        nnodes=0,
                        ncells=0,
                        ngrids=0,
                        # bbox=zeros(Int64, dims),
                        field=Dict{String, Dict{String, Any}}(),
                      _grid_origins=Array{Int64,1}[]
            ) = new(
                      dims,
                      grids,
                        grid_names,
                        nnodes,
                        ncells,
                        ngrids,
                        # bbox,
                        field,
                      _grid_origins
            )
end

function addgrid(self::MultiGrid, grid_name::String, grid::T) where {T<:AbstractGrid}

  if grid_name in self.grid_names
    error("Grid $grid_name already exists!")
  end

  self.nnodes += grid.nnodes
  self.ncells += grid.ncells

  push!(self.grids, grid)
  push!(self.grid_names, grid_name)
end

"Returns the requested grid"
function get_grid(self::MultiGrid, grid_name::String)

  grids_found = find(x->x==grid_name, self.grid_names)

  if length(grids_found)==0
    error("$grid_name not found. Grids available are $(self.grid_names).")
  elseif length(grids_found)!=1
    error("Logic error!")
  end

  return get_grid(self, grids_found[1])
end

"Returns the requested grid"
function get_grid(self::MultiGrid, grid_index::Int64)
  if grid_index>self.ngrids
    error("Requested invalid grid index $grid_index; max is $(self.ngrids).")
  end
  return self.grids[grid_index]
end


function get_node(self::MultiGrid, i::Int64)
  if i>self.nnodes
    error("Requested invalid node index $i; max is $(self.nnodes).")
  end

  ind = i
  for grid in self.grids
    if ind<=grid.nnodes
      return get_node(grid, ind)
    end
    ind -= grid.nnodes
  end

  error("Logic error!")
end
function get_node(self::MultiGrid, coor::Array{Int64,1})
  error("Not Implemented yet!")
end

function get_cell(self::MultiGrid, i::Int64)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  end

  ind = i
  for grid in self.grids
    if ind<=grid.ncells
      return get_ncell(grid, ind)
    end
    ind -= grid.ncells
  end

  error("Logic error!")

end
function get_cell(self::MultiGrid, coor::Array{Int64,1})
  error("Not Implemented yet!")
end

function get_fieldval(self::MultiGrid, field_name::String, i::Int64)
  if i>self.nnodes
    error("Requested invalid node index $i; max is $(self.nnodes).")
  end

  ind = i
  for grid in self.grids
    if ind<=grid.nnodes
      return get_fieldval(grid, field_name, ind)
    end
    ind -= grid.nnodes
  end

  error("Logic error!")
end
function get_fieldval(self::MultiGrid, field_name::String, coor::Array{Int64,1})
  error("Not Implemented yet!")
end

function add_field(self::MultiGrid, field_name::String, field_type::String,
                                                                        field_data)
  ind = 0
  for grid in self.grids
    add_field(grid, field_name, field_type, field_data[ind+1:ind+grid.nnodes])
    ind += grid.nnodes
  end
end
function calculate_field(self::MultiGrid, args...)
  for grid in self.grids
    calculate_fields(grid, args...)
  end
end
function lintransform!(self::MultiGrid, args...; optargs...)
  for grid in self.grids
    lintransform(grid, args...; optargs...)
  end
end
function transform!(self::MultiGrid, args...; optargs...)
end
function save(self::MultiGrid, args...; optargs...)
  for grid in self.grids
    save(grid, args...; optargs...)
  end
end
function plot(self::MultiGrid; optargs...)
  for grid in self.grids
    save(grid, args...; optargs...)
  end
end
##### INTERNAL FUNCTIONS  ######################################################

##### END OF MULTIGRID #########################################################
