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
  # Optional inputs
  loop_dim::Int64                     # Index of dimension to close in loop

  # Properties
  dims::Int64                         # Number of dimensions
  nnodes::Int64                       # Number of nodes
  ncells::Int64                       # Number of cells
  nodes::Array{T,2} where{T<:Real}    # Position of each node
  field::Dict{String, Dict{String, Any}}  # Calculated fields

  # Internal data
  _ndivsnodes::Tuple                  # Number of nodes in each coordinate
  _ndivscells::Tuple                  # Number of cells in each coordinate

  Grid(P_min, P_max, NDIVS, loop_dim=0,
              dims=_calc_dims(P_min),
                nnodes=_calc_nnodes(NDIVS, loop_dim),
                ncells=_calc_ncells(NDIVS),
                nodes=_generate_grid(P_min, P_max, NDIVS, loop_dim),
                field=Dict{String, Dict{String, Any}}(),
              _ndivsnodes=Tuple(_calc_ndivsnodes(NDIVS, loop_dim)),
                _ndivscells=Tuple(_calc_ndivs(NDIVS))
      ) = _check(P_min, P_max, NDIVS, loop_dim) ? new(P_min, P_max, NDIVS,
              loop_dim,
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
function get_cell(self::Grid, coor_in::Array{Int64,1})

  # Correct 0 coordinate in quasi-dimensions
  coor = [ self._ndivscells[i]==0 && coor_in[i]==0 ? 1 : coor_in[i]
                                                    for i in 1:length(coor_in) ]
  # ERROR CASES
  if length(coor)!=self.dims
    error("$(self.dims)-dimensional grid requires $(self.dims) coordinates,"*
            " got $(length(coor)).")
  end
  for (dim, i) in enumerate(coor)
    if i>self._ndivscells[dim]
      if i==1 && self._ndivscells[dim]==0
        nothing
      else
        error("Requested cell $coor but max cell in"*
              " $dim-dimension is $(self._ndivscells[dim])")
      end
    end
  end

  # Displacements according to VTK convention, or closed loop
  disps = ones(Int64, self.dims)
  lpdm = self.loop_dim
  if lpdm!=0 && coor[lpdm]==self._ndivscells[lpdm]
    disps[lpdm] = -self._ndivscells[lpdm] + 1
  end
  d1, d2, d3 = vcat(disps, ones(Int64, 3-self.dims))

  # Checks for quasi-dimensions (NDIVS[i]==0)
  qdims = length( [1 for ndiv in self._ndivscells if ndiv==0] )
  cdims = length(coor)
  dims = self.dims - qdims
  aux1 = zeros(Int64, qdims)

  if dims==1
    return [sub2ind(self._ndivsnodes, (coor+vcat([0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1], aux1))...)]

  elseif dims==2
    return [sub2ind(self._ndivsnodes, (coor+vcat([0,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,d2], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([0,d2], aux1))...)]

  elseif dims==3
    return [sub2ind(self._ndivsnodes, (coor+vcat([0,0,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,0,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,d2,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([0,d2,0], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([0,0,d3], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,0,d3], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([d1,d2,d3], aux1))...),
            sub2ind(self._ndivsnodes, (coor+vcat([0,d2,d3], aux1))...)]

  else
    error("Definition of $(self.ndims)-dimensional cells not implemented yet!")
  end
end

function get_cellcenter(self::Grid, args...)
  nodes = get_cell(self, args...)
  C = sum([get_node(self, node) for node in nodes])/size(nodes, 1)
  return C
end

"""
  `get_fieldval(grid, field_name, coor)`

Returns the value of node of coordinates `coor` (1-indexed) in the field
'field_name'.
"""
function get_fieldval(self::Grid, field_name::String, coor::Array{Int64,1})
  return get_fieldval(self, field_name, sub2ind(self._ndivsnodes, coor...))
end


"""
  `get_fieldval(grid, field_name, i)`

Returns the value of i-th node (1-indexed) in the field 'field_name'.
"""
function get_fieldval(self::Grid, field_name::String, i::Int64)
  if !(field_name in keys(self.field))
    error("Field $field_name doesn't exist."*
          " Available fields: $(keys(self.field))")
  end

  return self.field[field_name]["field_data"][i]
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

  # Checks for quasi-dimensions (NDIVS[i]==0)
  qdims = length( [1 for ndiv in grid._ndivscells if ndiv==0] )
  dims = grid.dims - qdims

  generateVTK(filename, points; cells=cells, point_data=point_data,
                                                  _griddims=dims, args...)
end

"Plots the grid on PyPlot"
function plot(grid::Grid; fig_name="gridplot", fontsize=15,
                          xlims=nothing, ylims=nothing, zlims=nothing,
                          labelcells=true, labelnodes=false, labelndivs=true,
                          title_str=nothing)

  if grid.dims>3
    error("There is no plotting method for $(grid.dims)-dimensional grids")
  end

  fig = PyPlot.figure(fig_name)
  ax = fig[:gca](projection="3d")
  vectors_to_plot = []

  nc = grid.ncells

  # Iterates over every child plotting them
  for i in 1:nc
    nodes = [vcat(get_node(grid, node), zeros(3-grid.dims))
                                                  for node in get_cell(grid, i)]

    if labelcells
      center = vcat(get_cellcenter(grid, i), zeros(3-grid.dims))
      PyPlot.text3D(center[1], center[2], center[3], "$i", fontsize=fontsize, color="g")
    end


    # Connects the nodes for visualization
    dims = grid.dims - length( [1 for ndiv in grid._ndivscells if ndiv==0] )
    if dims==3
      for j in 1:3
        push!(vectors_to_plot, [nodes[j], -nodes[j]+nodes[j+1]])
        push!(vectors_to_plot, [nodes[4+j], -nodes[4+j]+nodes[4+j+1]])
      end
      push!(vectors_to_plot, [nodes[4], -nodes[4]+nodes[1]])
      push!(vectors_to_plot, [nodes[8], -nodes[8]+nodes[5]])
      for j in 1:4
        push!(vectors_to_plot, [nodes[0+j], -nodes[0+j]+nodes[4+j]])
      end
    elseif dims==2
      for j in 1:3
        push!(vectors_to_plot, [nodes[j], -nodes[j]+nodes[j+1]])
      end
      push!(vectors_to_plot, [nodes[4], -nodes[4]+nodes[1]])
    else
      push!(vectors_to_plot, [nodes[1], -nodes[1]+nodes[2]])
    end
  end

  # Iterates over every dimension labeling the subdivision
  if labelndivs
    for i in 1:grid.dims
      for j in 1:grid._ndivscells[i]
        coor = ones(Int64, grid.dims)
        coor[i] = j
        p1 = get_node(grid, coor)
        if grid.loop_dim==i && j==grid._ndivscells[i]
          coor[i] = 1
        else
          coor[i] += 1
        end
        p2 = get_node(grid, coor)
        center = vcat((p1+p2)/2, zeros(3-grid.dims))
        PyPlot.text3D(center[1], center[2], center[3], "$j", fontsize=fontsize,
                                                                      color="r")
      end
    end
  end

  # Plots
  xyzuvw = [[] for _ in 1:6]
  for vector in vectors_to_plot
      _vector = []
      append!(_vector, vector[1])
      append!(_vector, vector[2])
      for (i,val) in enumerate(_vector)
          push!(xyzuvw[i], val)
      end
  end
  x,y,z,u,v,w = [coor for coor in xyzuvw]
  ax[:quiver](x,y,z, u,v,w, arrow_length_ratio=0.0);

  # Labels nodes
  if labelnodes
    for i in 1:grid.nnodes
      pos = vcat(get_node(grid, i), zeros(3-grid.dims))
      PyPlot.text3D(pos[1], pos[2], pos[3], "$i", fontsize=fontsize, color="k")
    end
  end


  # Format axes
  for (lim, lims, label, lbl) in [(PyPlot.xlim, xlims, PyPlot.xlabel, "x"),
                                  (PyPlot.ylim, ylims, PyPlot.ylabel, "y"),
                                  (PyPlot.zlim, zlims, PyPlot.zlabel, "z")]
    if lims!=nothing; lim(lims); end;
    label(lbl)

  end

  # Legend
  handles = []
  for (flag, clr, lbl) in [(labelnodes, "k", "Node"),
                            (labelcells, "g", "Cell"),
                            (labelndivs, "r", "Coordinate")]
    if flag; push!(handles, patch.Patch(color=clr, label=lbl)); end;
  end
  if size(handles,1)>0; PyPlot.legend(handles=handles); end;

  if title_str!=nothing; PyPlot.title(title_str); end;

end




##### INTERNAL FUNCTIONS  ######################################################
function _check(P_min::Array{T,1} where {T<:Real},
                P_max::Array{T,1} where {T<:Real} ,
                NDIVS::Array{T,1} where {T<:Any},
                loop_dim::Int64)

  # Error cases
  if size(P_min)!=size(P_max)
    error("`P_min` and `P_max` must have the same dimensions "*
                                          "($(size(P_min))!=$(size(P_max)))")
  elseif loop_dim>size(P_min,1) || loop_dim<0
    error("Invalid loop dimension $loop_dim"*
                                      " in $(size(P_min,1))-dimensional grid.")

  elseif typeof(NDIVS)==Array{Int64,1}
    if _calc_dims(P_min)!=length(NDIVS)
      error("Division for each dimension must be given "*
                                          "$(length(P_min))!=$(length(NDIVS))")
    end
    for (i, this_div) in enumerate(NDIVS)
      if this_div<0
        error("Invalid division $this_div in $i-dimension.")
      end
    end
  end
  for i in 1:_calc_dims(P_min)
    if NDIVS[i]==0
      if P_min[i]!=P_max[i]
        error("Value in `P_min` must be equal to `P_max` in quasi-dimension $i"*
                                            "($(P_min[i])!=$(P_max[i]))")
      end
    else
      if P_min[i]>=P_max[i]
        error("Values in `P_min` must be less than values in `P_max` "*"
                                            ($(P_min[i])>=$(P_max[i]))")
      end
    end
  end

  return true
end

function _calc_ndivs(NDIVS::Array{T,1} where {T<:Any})
  if typeof(NDIVS)==Array{Int64,1}
    return deepcopy(NDIVS)

  elseif typeof(NDIVS)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    return [ sum([sec[2] for sec in secs]) for secs in NDIVS ]

  else
      error("`NDIVS` expected to be type `Array{Int64,1}` or"*
              " `Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}`; got "*
              "$(typeof(NDIVS))")
  end

end


"Calculates total number of nodes given NDIVS"
function _calc_nnodes(NDIVS::Array{T,1} where {T<:Any}, loop_dim::Int64)
  ndivs = _calc_ndivsnodes(NDIVS, loop_dim)
  return prod([ elem==0 ? 1 : elem for elem in ndivs])
end

"Calculates number of nodes in each dimension"
function _calc_ndivsnodes(NDIVS::Array{T,1} where {T<:Any}, loop_dim::Int64)
  ndivs = _calc_ndivs(NDIVS)
  if loop_dim!=0
    ndivs[loop_dim] -= 1
  end
  return ndivs+1
end

"Calculates the number of cells given NDIVS"
function _calc_ncells(NDIVS::Array{T,1} where {T<:Any})
  return prod([ elem==0 ? 1 : elem for elem in _calc_ndivs(NDIVS)])
end

function _calc_dims(P::Array{T,1} where {T<:Real})
  return length(P)
end

function _generate_grid(P_min::Array{T,1} where {T<:Real},
                        P_max::Array{T,1} where {T<:Real},
                        NDIVS::Array{T,1} where {T<:Any},
                        loop_dim::Int64)

  dims = _calc_dims(P_min)
  nnodes = _calc_nnodes(NDIVS, loop_dim)
  nodes = zeros(nnodes, dims)
  ndivs = Tuple(_calc_ndivs(NDIVS))

  # Discretizes each coordinate according to NDIVS
  if typeof(NDIVS)==Array{Int64,1}
    spacing = [linspace(P_min[i], P_max[i], ndivs[i]+1) for i in 1:dims]

  elseif typeof(NDIVS)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    spacing = [multidiscretize(x->x, P_min[i], P_max[i], NDIVS[i]) for i in 1:dims]

  else
    error("Not implemented yet")
  end

  # Creates the grid
  ind2 = 1
  ndivsnodes = Tuple(_calc_ndivsnodes(NDIVS, 0))
  for ind in 1:_calc_nnodes(NDIVS, 0)
    sub = ind2sub(ndivsnodes, ind)

    if loop_dim==0 || sub[loop_dim]!=ndivsnodes[loop_dim]
      p = [spacing[dim][sub[dim]] for dim in 1:dims]
      nodes[ind2, :] = p
      ind2 += 1
    end
  end

  return nodes
end






##### END OF GRID ##############################################################
