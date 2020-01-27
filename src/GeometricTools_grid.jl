#=##############################################################################
# DESCRIPTION
    Methods for generation and evaluation of grids.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Mar 2018
  * License   : MIT License
=###############################################################################

struct Field{D}
    field_name::String
    field_type::String
    entry_type::String
    field_data::Vector{D}
end

################################################################################
# GRID
################################################################################
"""
  `Grid(pmin, pmax, ndivs)`

Generates an n-dimensional grid.

  **Arguments**
  * `pmin::Array{Float64,1}`   : Minimum point of the domain.
  * `pmax::Array{Float64,1}`   : Maximum point of the domain.
  * `ndivs::Array{Int64,1}`     : Number of divisions in each coordinate.

  **Properties**
  * `dims::Int64`               : Number of dimensions.
  * `nnodes::Int64`             : Number of nodes in the grid.
  * `nodes::Array{Float64,2}`   : Matrix size (`nnodes`, `dims`) of node position.
  * `field` : Contains calculated fields formated as field[field_name] = Dict(
                            "field_name" => field_name::String,
                            "field_type" => "scalar" or "vector",
                            "entry_type" => "node" or "cell",
                            "field_data" => data
                            )
            where `data` is an array data[i] = [val1, val2, ...] containing
            this field values (scalar or vector) at each node in the grid.

NOTE: All indexing is done linearly, meaning that `nodes` is indexed from 1 to
      `nnodes`, and all data fields follow the same indexing.

NOTE2: `ndivs` can either be an array of integers with ndivs[i] indicating the
      number of divisions in the i-th coordinate, or it can be an array of
      sections (see `multidiscretize()` doc) with ndivs[i] = [sec1, sec2, ...]
      indicating the discretization into sections in the i-th coordinate.
"""
struct Grid{TF,NDIV} <: AbstractGrid
  # User inputs
  pmin::Vector{TF}                       # Minimum point of the domain
  pmax::Vector{TF}                       # Maximum point of the domain
  ndivs::Vector{NDIV}                    # Number of divisions in each coordinate
  # Optional inputs
  loop_dim::Int                          # Index of dimension to close in loop
  # Properties
  dims::Int                              # Number of dimensions
  nnodes::Int                            # Number of nodes
  ncells::Int                            # Number of cells
  nodes::Matrix{TF}                      # Position of each node
  scalar_field::Dict{String, Field{Float64}}         # Calculated scalar fields
  vector_field::Dict{String, Field{Vector{Float64}}} # Calculated vector fields
  # Internal
  _ndivsnodes::Tuple                  # Number of nodes in each coordinate
  _ndivscells::Tuple                  # Number of cells in each coordinate
  _override_vtkcelltype::Int64        # Option for overriding vtk outputs

  Grid(pmin, pmax, ndivs, loop_dim=0,
       dims=_calc_dims(pmin),
       nnodes=_calc_nnodes(ndivs, loop_dim),
       ncells=_calc_ncells(ndivs),
       nodes=_generate_grid(pmin, pmax, ndivs, loop_dim),
       # bbox=_calc_ndivs(ndivs),
       field=Dict{String, Dict{String, Any}}(),
       _ndivsnodes=Tuple(_calc_ndivsnodes(ndivs, loop_dim)),
       _ndivscells=Tuple(_calc_ndivs(ndivs)),
       _override_vtkcelltype=-1
      ) = _check(pmin, pmax, ndivs, loop_dim) ? new(pmin, pmax, ndivs,
              loop_dim,
              dims,
                nnodes,
                ncells,
                nodes,
                # bbox,
                field,
              _ndivsnodes,
                _ndivscells,
                _override_vtkcelltype
      ) : error("Logic error!")
end

"""
  `get_node(grid, i)`

Returns the position of the i-th node (1-indexed) in the grid
"""
get_node(grid::Grid, i::Integer) = grid.nodes[:, i]

"""
  `get_node(grid, coor)`

Returns the position of the node of subscript coordinates `coor` (1-indexed)
"""
get_node(grid::Grid, coor::AbstractVector{Integer}) = get_node(grid, LinearIndices(self._ndivsnodes)[coor...])

"""
  `get_cell(grid, i)`

Returns the nodes indices of i-th cell in the grid (1-indexed)
"""
function get_cell(self::Grid, i::Integer)
  dims = Tuple(d != 0 ? d : 1 for d in self._ndivscells)
  return get_cell(self, collect(Tuple(CartesianIndices(dims)[i])))
end

"""
  `get_cell(grid, coor)`

Returns the node indices of the cell with subscript coordinates `coor`
(1-indexed). The format corresponds to VTK_HEXAHEDRON (=12) in 3D, VTK_QUAD (=9)
in 2D, or VTK_LINE (=3) in 1D---except that points are 1-indexed instead of
0-indexed.
"""
function get_cell(self::Grid, coor_in::AbstractVector{<:Integer})

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

  # Checks for quasi-dimensions (ndivs[i]==0)
  qdims = length( [1 for ndiv in self._ndivscells if ndiv==0] )
  cdims = length(coor)
  dims = self.dims - qdims
  aux1 = zeros(Int64, qdims)

  if dims==1
    return [LinearIndices(self._ndivsnodes)[(coor+vcat([0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1], aux1))...]]

  elseif dims==2
    return [LinearIndices(self._ndivsnodes)[(coor+vcat([0,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,d2], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([0,d2], aux1))...]]

  elseif dims==3
    return [LinearIndices(self._ndivsnodes)[(coor+vcat([0,0,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,0,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,d2,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([0,d2,0], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([0,0,d3], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,0,d3], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([d1,d2,d3], aux1))...],
            LinearIndices(self._ndivsnodes)[(coor+vcat([0,d2,d3], aux1))...]]

  else
    error("Definition of $(self.ndims)-dimensional cells not implemented yet!")
  end
end


"Plots the grid on PyPlot"
function plot(grid::Grid; fig_name="gridplot", fontsize=15,
                          xlims=nothing, ylims=nothing, zlims=nothing,
                          labelcells=true, labelnodes=false, labelndivs=true,
                          title_str=nothing,
                          alpha=1.0)

  if grid.dims>3
    error("There is no plotting method for $(grid.dims)-dimensional grids")
  end

  fig = PyPlot.figure(fig_name)
  ax = fig.gca(projection="3d")
  vectors_to_plot = []

  nc = grid.ncells

  # Iterates over every child plotting them
  for i in 1:nc
    nodes = [vcat(get_node(grid, node), zeros(Float64, 3-grid.dims))
                                                  for node in get_cell(grid, i)]

    if labelcells
      center = vcat(get_cellcenter(grid, i), zeros(Float64, 3-grid.dims))
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
        center = vcat((p1+p2)/2, zeros(Float64, 3-grid.dims))
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
  ax.quiver(x,y,z, u,v,w, arrow_length_ratio=0.0, alpha=alpha);

  # Labels nodes
  if labelnodes
    for i in 1:grid.nnodes
      pos = vcat(get_node(grid, i), zeros(Float64, 3-grid.dims))
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


"""
  `lintransform!(grid::Grid, M::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates the grid by the rotation matrix `M` and translation
vector `T` (linear transformation).
"""
function lintransform!(grid::Grid, M::Array{P,2}, T::Array{P,1};
                        reset_fields::Bool=true) where{P<:Real}
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
    grid.nodes[:,i] = countertransform(grid.nodes[:,i], invM, T)
  end
end

"""
  `transform!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
position of every node is given to the function `f`.
"""
function transform!(grid::Grid, f; reset_fields::Bool=true)

    if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

    for i in 1:grid.nnodes
      grid.nodes[:,i] = f(grid.nodes[:,i])
    end
end


"""
  `transform2!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
indices of every node is given to the function `f`.
"""
function transform2!(grid::Grid, f; reset_fields::Bool=true)

    if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

    for i in 1:grid.nnodes
      grid.nodes[:,i] = f(CartesianIndices(grid._ndivsnodes)[i])
    end
end

"""
  `transform3!(grid::Grid, f)`

Applies the space transformation given by function `f` to the grid, where the
both the position and indices of every node is given to the function `f`.
"""
function transform3!(grid::Grid, f; reset_fields::Bool=true)

    if reset_fields; grid.field = Dict{String, Dict{String, Any}}(); end;

    for i in 1:grid.nnodes
      grid.nodes[:,i] = f(grid.nodes[:,i], CartesianIndices(grid._ndivsnodes)[i])
    end
end


##### INTERNAL FUNCTIONS  ######################################################
function _check(pmin::Array{T,1} where {T<:Real},
                pmax::Array{T,1} where {T<:Real} ,
                ndivs::Array{T,1} where {T<:Any},
                loop_dim::Int64)

  # Error cases
  if size(pmin)!=size(pmax)
    error("`pmin` and `pmax` must have the same dimensions "*
                                          "($(size(pmin))!=$(size(pmax)))")
  elseif loop_dim>size(pmin,1) || loop_dim<0
    error("Invalid loop dimension $loop_dim"*
                                      " in $(size(pmin,1))-dimensional grid.")

  elseif typeof(ndivs)==Array{Int64,1}
    if _calc_dims(pmin)!=length(ndivs)
      error("Division for each dimension must be given "*
                                          "$(length(pmin))!=$(length(ndivs))")
    end
    for (i, this_div) in enumerate(ndivs)
      if this_div<0
        error("Invalid division $this_div in $i-dimension.")
      end
    end
  end
  _ndivs = _calc_ndivs(ndivs)
  for i in 1:_calc_dims(pmin)
    if _ndivs[i]==0
      if pmin[i]!=pmax[i]
        error("Value in `pmin` must be equal to `pmax` in quasi-dimension $i"*
                                            "($(pmin[i])!=$(pmax[i]))")
      end
    else
      if pmin[i]>=pmax[i]
        error("Values in `pmin` must be less than values in `pmax` "*"
                                            ($(pmin[i])>=$(pmax[i]))")
      end
    end
  end

  return true
end

function _calc_ndivs(ndivs::Array{T,1} where {T<:Any})
  if typeof(ndivs)==Array{Int64,1}
    return deepcopy(ndivs)

  elseif typeof(ndivs)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    return [ sum([sec[2] for sec in secs]) for secs in ndivs ]

  else
      error("`ndivs` expected to be type `Array{Int64,1}` or"*
              " `Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}`; got "*
              "$(typeof(ndivs))")
  end

end


"Calculates total number of nodes given ndivs"
function _calc_nnodes(ndivs::Array{T,1} where {T<:Any}, loop_dim::Int64)
  ndivs = _calc_ndivsnodes(ndivs, loop_dim)
  return prod([ elem==0 ? 1 : elem for elem in ndivs])
end

"Calculates number of nodes in each dimension"
function _calc_ndivsnodes(ndivs::Array{T,1} where {T<:Any}, loop_dim::Int64)
  ndivs = _calc_ndivs(ndivs)
  if loop_dim!=0
    ndivs[loop_dim] -= 1
  end
  return ndivs .+ 1
end

"Calculates the number of cells given ndivs"
function _calc_ncells(ndivs::Array{T,1} where {T<:Any})
  return prod([ elem==0 ? 1 : elem for elem in _calc_ndivs(ndivs)])
end

function _calc_dims(P::Array{T,1} where {T<:Real})
  return length(P)
end

function _generate_grid(pmin::Array{T,1} where {T<:Real},
                        pmax::Array{T,1} where {T<:Real},
                        ndivs::Array{T,1} where {T<:Any},
                        loop_dim::Int64)

  dims = _calc_dims(pmin)
  nnodes = _calc_nnodes(ndivs, loop_dim)
  nodes = zeros(Float64, dims, nnodes)
  ndivs = Tuple(_calc_ndivs(ndivs))

  # Discretizes each coordinate according to ndivs
  if typeof(ndivs)==Array{Int64,1}
    spacing = [range(pmin[i], stop=pmax[i], length=ndivs[i]+1) for i in 1:dims]

  elseif typeof(ndivs)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    spacing = [multidiscretize(x->x, pmin[i], pmax[i], ndivs[i]) for i in 1:dims]

  else
    error("Not implemented yet")
  end

  # Creates the grid
  ind2 = 1
  ndivsnodes = Tuple(_calc_ndivsnodes(ndivs, 0))
  for ind in 1:_calc_nnodes(ndivs, 0)
    sub = CartesianIndices(ndivsnodes)[ind]

    if loop_dim==0 || sub[loop_dim]!=ndivsnodes[loop_dim]
      p = [spacing[dim][sub[dim]] for dim in 1:dims]
      nodes[:, ind2] = p
      ind2 += 1
    end
  end

  return nodes
end






##### END OF GRID ##############################################################
