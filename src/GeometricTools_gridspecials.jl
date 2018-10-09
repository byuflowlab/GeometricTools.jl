#=##############################################################################
# DESCRIPTION
    Special grid types
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jun 2018
  * License   : MIT License
=###############################################################################

################################################################################
# TRIANGULAR SURFACE GRID TYPE
################################################################################

"""
  `GridTriangleSurface(orggrid, dimsplit)`

Receives a 3D surface grid (like the one in `Paneled Wing Example` in the
documentation), which by construction is made of nonplanar quadrilateral panels,
and creates a surface grid of planar triangular panels by splitting every
original quadrilateral panel into triangles.

  **Arguments**
  * `orggrid`         : Original quadrilateral surface grid.
  * `dimsplit`        : Dimension along which to split the quadrilaterals.
"""
mutable struct GridTriangleSurface <: AbstractGrid

  # User inputs
  orggrid::Grid                       # Original quadrilateral surface grid
  dimsplit::Int64                     # Dimension to split

  # Properties
  dims::Int64                         # Number of dimensions
  nnodes::Int64                       # Number of nodes
  ncells::Int64                       # Number of cells
  field::Dict{String, Dict{String, Any}}  # Calculated fields

  # Internal data
  _ndivsnodes::Tuple                  # Number of nodes in each coordinate
  _ndivscells::Tuple                  # Number of cells in each coordinate
  _override_vtkcelltype::Int64        # Option for overriding vtk outputs

  GridTriangleSurface(orggrid, dimsplit,
                      dims=orggrid.dims,
                        nnodes=orggrid.nnodes,
                        ncells=2*orggrid.ncells,
                        field=Dict{String, Dict{String, Any}}(),
                      _ndivsnodes=orggrid._ndivsnodes,
                        _ndivscells=_ndivscells(orggrid, dimsplit),
                        _override_vtkcelltype=5
      ) = _checkGridTriangleSurface(orggrid, dimsplit) ? new(orggrid, dimsplit,
                      dims,
                        nnodes,
                        ncells,
                        field,
                      _ndivsnodes,
                        _ndivscells,
                        _override_vtkcelltype
      ) : error("Logic error!")
end


get_node(self::GridTriangleSurface, i::Int64) = get_node(self.orggrid, i)
get_node(self::GridTriangleSurface, coor::Array{Int64,1}) = get_node(self.orggrid, coor)

function get_cell(self::GridTriangleSurface, i::Int64)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  end
  return get_cell(self, collect(ind2sub(self._ndivscells, i)))
end

function get_cell(self::GridTriangleSurface, coor::Array{Int64,1})
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

  # Converts this coordinates to the coordinates of the quadrilateral panel
  quadcoor = Int64[ceil(ind/2^(i==self.dimsplit)) for (i,ind) in enumerate(coor)]

  # Gets the nodes of the quadrilateral panel
  quadnodes = get_cell(self.orggrid, quadcoor)

  # Splits the quadrilateral into a triangle
  if coor[self.dimsplit]%2!=0
    return [quadnodes[1], quadnodes[2], quadnodes[3]]
  else
    return [quadnodes[3], quadnodes[4], quadnodes[1]]
  end
end


"""
  `get_normal(self::GridTriangleSurface, i::Int64)`

Returns the normal vector of the i-th panel.
"""
function get_normal(self::GridTriangleSurface, i::Int64)
  return _calc_normal(get_cellnodes(self, i))
end
function get_normal(self::GridTriangleSurface, coor::Array{Int64,1})
  return get_normal(self, sub2ind(self._ndivsnodes, coor...))
end


"""
  `get_tangent(self::GridTriangleSurface, i::Int64)`

Returns the tangential vector of the i-th panel.
"""
function get_tangent(self::GridTriangleSurface, i::Int64)
  return _calc_tangent(get_cellnodes(self, i))
end
function get_tangent(self::GridTriangleSurface, coor::Array{Int64,1})
  return get_tangent(self, sub2ind(self._ndivsnodes, coor...))
end

"""
  `get_unitvectors(self::GridTriangleSurface, i::Int64)`

Returns the orthogonal system `(t, o, n)` of the i-th panel, with `t` the
tangent vector, `o` the oblique vector (which is also tangent), and `n` the
normal vector.
"""
function get_unitvectors(self::GridTriangleSurface, i::Int64)
  return _calc_unitvectors(get_cellnodes(self, i))
end
function get_unitvectors(self::GridTriangleSurface, coor::Array{Int64,1})
  return get_unitvectors(self, sub2ind(self._ndivsnodes, coor...))
end

function lintransform!(self::GridTriangleSurface, args...; optargs...)
  lintransform!(self.orggrid, args...; optargs...)
end

function transform!(self::GridTriangleSurface, args...; optargs...)
  transform!(self.orggrid, args...; optargs...)
end
function transform2!(self::GridTriangleSurface, args...; optargs...)
  transform!(self.orggrid, args...; optargs...)
end
function transform3!(self::GridTriangleSurface, args...; optargs...)
  transform!(self.orggrid, args...; optargs...)
end


##### INTERNAL FUNCTIONS  ######################################################
function _checkGridTriangleSurface(orggrid::Grid, dimsplit::Int64)
  # Number of quasi-dimensions
  qdims = length( [1 for ndiv in orggrid._ndivscells if ndiv==0] )

  if qdims!=1
    error("Expected one quasi-dimension, found $qdims.")

  elseif orggrid.dims!=3
    error("Expected a three-dimensional grid, got $(orggrid.dims)-dimensional.")

  elseif !(dimsplit in [1,2])
    error("Invalid split dimension $dimsplit")

  end

  return true
end

function _ndivscells(orggrid::Grid, dimsplit::Int64)
  return Tuple([
            2^(i==dimsplit)*divs for (i,divs) in enumerate(orggrid._ndivscells)
              ])
end

function _calc_tangent(nodes::Array{Array{T,1},1}) where{T<:Real}
  t = nodes[2] - nodes[1]
  return t/norm(t)
end

function _calc_normal(nodes::Array{Array{T,1},1}) where{T<:Real}
  n = cross(nodes[2]-nodes[1], nodes[3]-nodes[1])
  return n/norm(n)
end

function _calc_unitvectors(nodes::Array{Array{T,1},1}) where{T<:Real}
  t = _calc_tangent(nodes)
  n = _calc_normal(nodes)
  o = cross(n,t)
  return t, o, n
end
##### END OF TRIANGULAR GRID ###################################################
