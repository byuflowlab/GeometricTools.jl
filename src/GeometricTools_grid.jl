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

NOTE2: `NDIVS` can either be an array of integers with NDIVS[i] indicating the
      number of divisions in the i-th coordinate, or it can be an array of
      sections (see `multidiscretize()` doc) with NDIVS[i] = [sec1, sec2, ...]
      indicating the discretization into sections in the i-th coordinate.
"""
mutable struct Grid <: AbstractGrid

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
  # bbox::Array{Int64, 1}               # Bounding box (cells in each dimension)
  field::Dict{String, Dict{String, Any}}  # Calculated fields

  # Internal data
  _ndivsnodes::Tuple                  # Number of nodes in each coordinate
  _ndivscells::Tuple                  # Number of cells in each coordinate
  _override_vtkcelltype::Int64        # Option for overriding vtk outputs

  Grid(P_min, P_max, NDIVS, loop_dim=0,
              dims=_calc_dims(P_min),
                nnodes=_calc_nnodes(NDIVS, loop_dim),
                ncells=_calc_ncells(NDIVS),
                nodes=_generate_grid(P_min, P_max, NDIVS, loop_dim),
                # bbox=_calc_ndivs(NDIVS),
                field=Dict{String, Dict{String, Any}}(),
                _ndivsnodes=Tuple(_calc_ndivsnodes(NDIVS, loop_dim)),
                _ndivscells=Tuple(_calc_ndivs(NDIVS)),
                _override_vtkcelltype=-1
      ) = _check(P_min, P_max, NDIVS, loop_dim) ? new(P_min, P_max, NDIVS,
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
function get_node(self::Grid, i::Int64)
  if i>self.nnodes
    error("Requested invalid node index $i; max is $(self.nnodes).")
  elseif i<1
    error("Invalid index $i (it must be greater than 0).")
  end

  return self.nodes[:, i]
end

"""
  `get_node(grid, coor)`

Returns the position of the node of subscript coordinates `coor` (1-indexed)
"""
function get_node(self::Grid, coor::Array{Int64,1})
  return get_node(self, LinearIndices(self._ndivsnodes)[coor...])
end

"""
  `get_cell(grid, i)`

Returns the nodes indices of i-th cell in the grid (1-indexed)
"""
function get_cell(self::Grid, i::Integer)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  elseif i<1
    error("Invalid index $i (it must be greater than 0).")
  end
  dims = Tuple(d != 0 ? d : 1 for d in self._ndivscells)
  # return get_cell(self, collect(ind2sub(self._ndivscells, i)))
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

  # Checks for quasi-dimensions (NDIVS[i]==0)
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

"Same than `get_cell(...)` but with low memory allocation using tuples."
function get_cell_t!(out, self::Grid, coor_in, lin, ndivscells)
  # ERROR CASES: Incorrect number of coordinates
  if length(coor_in)!=self.dims
    error("$(self.dims)-dimensional grid requires $(self.dims) coordinates,"*
            " got $(length(coor_in)).")
  end

  # @inbounds begin

      # Correct 0 coordinate in quasi-dimensions
      # coor = ( c==0 && self._ndivscells[i]==0 ? 1 : c for (i, c) in enumerate(coor_in))
      c1::Int = length(coor_in) >= 1 ? coor_in[1]==0 && ndivscells[1]==0 ? 1 : coor_in[1] : -1
      c2::Int = length(coor_in) >= 2 ? coor_in[2]==0 && ndivscells[2]==0 ? 1 : coor_in[2] : -1
      c3::Int = length(coor_in) >= 3 ? coor_in[3]==0 && ndivscells[3]==0 ? 1 : coor_in[3] : -1

      # ERROR CASE: Coordinates out of bounds
      if (c1 > ndivscells[1] && !(c1==1 && ndivscells[1]==0)) ||
         (c2 > ndivscells[2] && !(c2==1 && ndivscells[2]==0)) ||
         (c3 > ndivscells[3] && !(c3==1 && ndivscells[3]==0)) ||
         c1<=0 || c2<=0 || c3<=0
            error("Requested cell $(coor_in) but max cell is $(ndivscells). $((c1, c2, c3))")
      end

      # Displacements according to VTK convention, or closed loop
      # disps = ones(Int64, self.dims)
      # lpdm = self.loop_dim
      # if lpdm!=0 && coor[lpdm]==ndivscells[lpdm]
      #   disps[lpdm] = -ndivscells[lpdm] + 1
      # end
      # d1, d2, d3 = vcat(disps, ones(Int64, 3-self.dims))

      lpdm = self.loop_dim
      d1::Int = self.dims >= 1 && lpdm==1 && c1==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1
      d2::Int = self.dims >= 2 && lpdm==2 && c2==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1
      d3::Int = self.dims >= 3 && lpdm==3 && c3==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1

      # Checks for quasi-dimensions (NDIVS[i]==0)
      qdims = 0
      for ndiv in ndivscells; if ndiv==0; qdims += 1; end; end;
      cdims = length(coor_in)
      dims = self.dims - qdims

      # aux1 = zeros(Int64, qdims)
      # lin = LinearIndices(self._ndivsnodes)

      if dims==1

        # return [lin[(coor+vcat([0], aux1))...],
        #         lin[(coor+vcat([d1], aux1))...]]

        if cdims==1
            out[1] = lin[c1]
            out[2] = lin[c1+d1]
            return out
        elseif cdims==2
            out[1] = lin[c1, c2]
            out[2] = lin[c1+d1, c2]
            return out
        elseif cdims==3
            out[1] = lin[c1, c2, c3]
            out[2] = lin[c1+d1, c2, c3]
            return out
        else
            error("Logic error: cdims=$cdims")
        end

      elseif dims==2

        # return [lin[(coor+vcat([0,0], aux1))...],
        #         lin[(coor+vcat([d1,0], aux1))...],
        #         lin[(coor+vcat([d1,d2], aux1))...],
        #         lin[(coor+vcat([0,d2], aux1))...]]

        if cdims==2
            out[1] = lin[c1, c2]
            out[2] = lin[c1+d1, c2]
            out[3] = lin[c1+d1, c2+d2]
            out[4] = lin[c1, c2+d2]
            return out
        elseif cdims==3
            out[1] = lin[c1, c2, c3]
            out[2] = lin[c1+d1, c2, c3]
            out[3] = lin[c1+d1, c2+d2, c3]
            out[4] = lin[c1, c2+d2, c3]
            return out
        else
            error("Logic error: cdims=$cdims")
        end

      elseif dims==3
        # return [lin[(coor+vcat([0,0,0], aux1))...],
        #         lin[(coor+vcat([d1,0,0], aux1))...],
        #         lin[(coor+vcat([d1,d2,0], aux1))...],
        #         lin[(coor+vcat([0,d2,0], aux1))...],
        #         lin[(coor+vcat([0,0,d3], aux1))...],
        #         lin[(coor+vcat([d1,0,d3], aux1))...],
        #         lin[(coor+vcat([d1,d2,d3], aux1))...],
        #         lin[(coor+vcat([0,d2,d3], aux1))...]]
        out[1] = lin[c1, c2, c3]
        out[2] = lin[c1+d1, c2, c3]
        out[3] = lin[c1+d1, c2+d2, c3]
        out[4] = lin[c1, c2+d2, c3]
        out[5] = lin[c1, c2, c3+d3]
        out[6] = lin[c1+d1, c2, c3+d3]
        out[7] = lin[c1+d1, c2+d2, c3+d3]
        out[8] = lin[c1, c2+d2, c3+d3]
        return out

      else
        error("Definition of $(self.ndims)-dimensional cells not implemented yet!")
      end

  # end
end

"Same than `get_cell(...)` but with low memory allocation using tuples."
function get_cell_t(self::Grid, coor_in, nodei::Int, lin, ndivscells)
  # ERROR CASES: Incorrect number of coordinates
  if length(coor_in)!=self.dims
    error("$(self.dims)-dimensional grid requires $(self.dims) coordinates,"*
            " got $(length(coor_in)).")
  end

  # @inbounds begin

      # Correct 0 coordinate in quasi-dimensions
      # coor = ( c==0 && self._ndivscells[i]==0 ? 1 : c for (i, c) in enumerate(coor_in))
      c1::Int = length(coor_in) >= 1 ? coor_in[1]==0 && ndivscells[1]==0 ? 1 : coor_in[1] : -1
      c2::Int = length(coor_in) >= 2 ? coor_in[2]==0 && ndivscells[2]==0 ? 1 : coor_in[2] : -1
      c3::Int = length(coor_in) >= 3 ? coor_in[3]==0 && ndivscells[3]==0 ? 1 : coor_in[3] : -1

      # ERROR CASE: Coordinates out of bounds
      if (c1 > ndivscells[1] && !(c1==1 && ndivscells[1]==0)) ||
         (c2 > ndivscells[2] && !(c2==1 && ndivscells[2]==0)) ||
         (c3 > ndivscells[3] && !(c3==1 && ndivscells[3]==0)) ||
         c1<=0 || c2<=0 || c3<=0
            error("Requested cell $(coor_in) but max cell is $(ndivscells). $((c1, c2, c3))")
      end

      # Displacements according to VTK convention, or closed loop
      # disps = ones(Int64, self.dims)
      # lpdm = self.loop_dim
      # if lpdm!=0 && coor[lpdm]==ndivscells[lpdm]
      #   disps[lpdm] = -ndivscells[lpdm] + 1
      # end
      # d1, d2, d3 = vcat(disps, ones(Int64, 3-self.dims))

      lpdm = self.loop_dim
      d1::Int = self.dims >= 1 && lpdm==1 && c1==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1
      d2::Int = self.dims >= 2 && lpdm==2 && c2==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1
      d3::Int = self.dims >= 3 && lpdm==3 && c3==ndivscells[lpdm] ? -ndivscells[lpdm] + 1 : 1

      # Checks for quasi-dimensions (NDIVS[i]==0)
      qdims = 0
      for ndiv in ndivscells; if ndiv==0; qdims += 1; end; end;
      cdims = length(coor_in)
      dims = self.dims - qdims

      # aux1 = zeros(Int64, qdims)
      # lin = LinearIndices(self._ndivsnodes)

      if dims==1

        if nodei < 1 || nodei > 2
            error("Requested invalid node index $(nodei); max is 2.")
        end

        # return [lin[(coor+vcat([0], aux1))...],
        #         lin[(coor+vcat([d1], aux1))...]]

        if cdims==1
            return nodei==1 ?   lin[c1] :
                                lin[c1+d1]
        elseif cdims==2
            return nodei==1 ?   lin[c1, c2] :
                                lin[c1+d1, c2]
        elseif cdims==3
            return nodei==1 ?   lin[c1, c2, c3] :
                                lin[c1+d1, c2, c3]
        else
            error("Logic error: cdims=$cdims")
        end

      elseif dims==2

        if nodei < 1 || nodei > 4
            error("Requested invalid node index $(nodei); max is 4.")
        end

        # return [lin[(coor+vcat([0,0], aux1))...],
        #         lin[(coor+vcat([d1,0], aux1))...],
        #         lin[(coor+vcat([d1,d2], aux1))...],
        #         lin[(coor+vcat([0,d2], aux1))...]]

        if cdims==2
            return  nodei==1 ? lin[c1, c2] :
                    nodei==2 ? lin[c1+d1, c2] :
                    nodei==3 ? lin[c1+d1, c2+d2] :
                               lin[c1, c2+d2]
        elseif cdims==3
            return  nodei==1 ? lin[c1, c2, c3] :
                    nodei==2 ? lin[c1+d1, c2, c3] :
                    nodei==3 ? lin[c1+d1, c2+d2, c3] :
                               lin[c1, c2+d2, c3]
        else
            error("Logic error: cdims=$cdims")
        end

      elseif dims==3
        # return [lin[(coor+vcat([0,0,0], aux1))...],
        #         lin[(coor+vcat([d1,0,0], aux1))...],
        #         lin[(coor+vcat([d1,d2,0], aux1))...],
        #         lin[(coor+vcat([0,d2,0], aux1))...],
        #         lin[(coor+vcat([0,0,d3], aux1))...],
        #         lin[(coor+vcat([d1,0,d3], aux1))...],
        #         lin[(coor+vcat([d1,d2,d3], aux1))...],
        #         lin[(coor+vcat([0,d2,d3], aux1))...]]
        return  nodei==1 ? lin[c1, c2, c3] :
                nodei==2 ? lin[c1+d1, c2, c3] :
                nodei==3 ? lin[c1+d1, c2+d2, c3] :
                nodei==4 ? lin[c1, c2+d2, c3] :
                nodei==5 ? lin[c1, c2, c3+d3] :
                nodei==6 ? lin[c1+d1, c2, c3+d3] :
                nodei==7 ? lin[c1+d1, c2+d2, c3+d3] :
                           lin[c1, c2+d2, c3+d3]
      else
        error("Definition of $(self.ndims)-dimensional cells not implemented yet!")
      end

  # end
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
      grid.nodes[:,i] .= f(view(grid.nodes, :, i))
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
  _ndivs = _calc_ndivs(NDIVS)
  for i in 1:_calc_dims(P_min)
    if _ndivs[i]==0
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
  return ndivs .+ 1
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
  nodes = Array{Float64}(undef, dims, nnodes)
  ndivs = Tuple(_calc_ndivs(NDIVS))

  # Discretizes each coordinate according to NDIVS
  if typeof(NDIVS)==Array{Int64,1}
    spacing = [range(P_min[i], stop=P_max[i], length=ndivs[i]+1) for i in 1:dims]

  elseif typeof(NDIVS)==Array{Array{Tuple{Float64,Int64,Float64,Bool},1},1}
    spacing = [multidiscretize(x->x, P_min[i], P_max[i], NDIVS[i]) for i in 1:dims]

  else
    error("Not implemented yet")
  end

  # Creates the grid
  ind2 = 1
  ndivsnodes = Tuple(_calc_ndivsnodes(NDIVS, 0))
  for ind in 1:_calc_nnodes(NDIVS, 0)
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
