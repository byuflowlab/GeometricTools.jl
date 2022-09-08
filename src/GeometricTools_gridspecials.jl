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
  dims = Tuple(d != 0 ? d : 1 for d in self._ndivscells)
  return get_cell(self, collect(Tuple(CartesianIndices(dims)[i])))
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


function get_cell_t!(tri_out, tricoor, quadcoor, quad_out, self::GridTriangleSurface,
                                                i::Int64, lin, ndivscells, cin)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  end

  for j in 1:length(tricoor); tricoor[j] = cin[i][j]; end;
  return get_cell_t!(tri_out, quadcoor, quad_out, self, tricoor, lin, ndivscells)
end

function get_cell_t!(tri_out, quadcoor, quad_out, self::GridTriangleSurface, coor,
                                                                lin, ndivscells)
  # ERROR CASES
  if length(coor)!=self.dims
    error("$(self.dims)-dimensional grid requires $(self.dims) coordinates,"*
            " got $(length(coor)).")
  end
  for (dim, i) in enumerate(coor)
    if i>ndivscells[dim]
      if i==1 && ndivscells[dim]==0
        nothing
      else
        error("Requested cell $coor but max cell in"*
              " $dim-dimension is $(ndivscells[dim])")
      end
    end
  end


  # Converts this coordinates to the coordinates of the quadrilateral panel
  # quadcoor = Int64[ceil(ind/2^(i==self.dimsplit)) for (i,ind) in enumerate(coor)]
  for (i, ind) in enumerate(coor)
      quadcoor[i] = ceil(Int, ind/2^(i==self.dimsplit))
  end

  # Gets the nodes of the quadrilateral panel
  # quadnodes = get_cell(self.orggrid, quadcoor)
  get_cell_t!(quad_out, self.orggrid, quadcoor, lin, ndivscells)

  # Splits the quadrilateral into a triangle
  if coor[self.dimsplit]%2!=0
    # return [quadnodes[1], quadnodes[2], quadnodes[3]]
    tri_out[1] = quad_out[1]
    tri_out[2] = quad_out[2]
    tri_out[3] = quad_out[3]
    return tri_out
  else
    # return [quadnodes[3], quadnodes[4], quadnodes[1]]
    tri_out[1] = quad_out[3]
    tri_out[2] = quad_out[4]
    tri_out[3] = quad_out[1]
    return tri_out
  end
end

"""
    `get_area(self::GridTriangleSurface, i_or_coor::Union{Int, Array{Int,1}})`
Returns the area of the i-th cell.
"""
function get_area(self::GridTriangleSurface, i)
    cell = get_cell(self, i)
    A = get_node(self, cell[1])
    B = get_node(self, cell[2])
    C = get_node(self, cell[3])
    return 0.5*norm(cross(B-A, C-A))
end
function _get_area(nodes, panel)

    p1, p2, p3 = panel

    crss1 = (nodes[2, p2]-nodes[2, p1])*(nodes[3, p3]-nodes[3, p1]) - (nodes[3, p2]-nodes[3, p1])*(nodes[2, p3]-nodes[2, p1])
    crss2 = (nodes[3, p2]-nodes[3, p1])*(nodes[1, p3]-nodes[1, p1]) - (nodes[1, p2]-nodes[1, p1])*(nodes[3, p3]-nodes[3, p1])
    crss3 = (nodes[1, p2]-nodes[1, p1])*(nodes[2, p3]-nodes[2, p1]) - (nodes[2, p2]-nodes[2, p1])*(nodes[1, p3]-nodes[1, p1])

    return 0.5*sqrt(crss1^2 + crss2^2 + crss3^2)
end

"""
    `get_area(self::GridTriangleSurface)`
Returns the area of the entire grid.
"""
function get_area(self::GridTriangleSurface)
    A = 0
    for i in 1:self.ncells
        A += get_area(self, i)
    end
    return A
end

"""
    `get_volume(self::GridTriangleSurface)`
Returns the volume encloused by the grid using Green's theorem. See
https://en.wikipedia.org/wiki/Green%27s_theorem#Area_Calculation
https://en.wikipedia.org/wiki/Polyhedron#Volume
"""
function get_volume(self::GridTriangleSurface)
    V = 0
    for i in 1:self.ncells
        A = get_area(self, i)
        n = get_normal(self, i)
        r = get_cellcenter(self, i)
        V += 1/3 * dot(r, A*n)
    end
    return V
end

"""
    `get_centroid(self::GridTriangleSurface)`
Returns the centroid of the volume encloused by the grid using Green's theorem.
See https://en.wikipedia.org/wiki/Green%27s_theorem#Area_Calculation
"""
function get_centroid(self::GridTriangleSurface)
    R = zeros(self.dims)
    for i in 1:self.ncells
        A = get_area(self, i)
        n = get_normal(self, i)
        r = get_cellcenter(self, i)
        R += 3/4*r * (1/3 * dot(r, A*n))
    end
    return R/get_volume(self)
end

"""
  `get_normal(self::GridTriangleSurface, i::Int64)`

Returns the normal vector of the i-th panel.
"""
function get_normal(self::GridTriangleSurface, i::Int64)
  return _calc_normal(get_cellnodes(self, i))
end
function get_normal(self::GridTriangleSurface, coor::Array{Int64,1})
  return get_normal(self, Base._sub2ind(self._ndivsnodes, coor...))
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

_calc_tnorm(nodes, panel) = sqrt((nodes[1, panel[2]] - nodes[1, panel[1]])^2 + (nodes[2, panel[2]] - nodes[2, panel[1]])^2 + (nodes[3, panel[2]] - nodes[3, panel[1]])^2)
_calc_t1(nodes, panel) = (nodes[1, panel[2]] - nodes[1, panel[1]]) / _calc_tnorm(nodes, panel)
_calc_t2(nodes, panel) = (nodes[2, panel[2]] - nodes[2, panel[1]]) / _calc_tnorm(nodes, panel)
_calc_t3(nodes, panel) = (nodes[3, panel[2]] - nodes[3, panel[1]]) / _calc_tnorm(nodes, panel)

_calc_tnorm(nodes) = sqrt((nodes[2][1] - nodes[1][1])^2 + (nodes[2][2] - nodes[1][2])^2 + (nodes[2][3] - nodes[1][3])^2)
_calc_t1(nodes) = (nodes[2][1] - nodes[1][1]) / _calc_tnorm(nodes)
_calc_t2(nodes) = (nodes[2][2] - nodes[1][2]) / _calc_tnorm(nodes)
_calc_t3(nodes) = (nodes[2][3] - nodes[1][3]) / _calc_tnorm(nodes)
function _calc_tangent(nodes::Array{Arr1,1}) where{Arr1<:AbstractArray}
    # t = nodes[2] - nodes[1]
    # return t/norm(t)
    return [_calc_t1(nodes), _calc_t2(nodes), _calc_t3(nodes)]
end

_calc_n1aux(nodes, panel) = (nodes[2, panel[2]]-nodes[2, panel[1]])*(nodes[3, panel[3]]-nodes[3, panel[1]]) - (nodes[3, panel[2]]-nodes[3, panel[1]])*(nodes[2, panel[3]]-nodes[2, panel[1]])
_calc_n2aux(nodes, panel) = (nodes[3, panel[2]]-nodes[3, panel[1]])*(nodes[1, panel[3]]-nodes[1, panel[1]]) - (nodes[1, panel[2]]-nodes[1, panel[1]])*(nodes[3, panel[3]]-nodes[3, panel[1]])
_calc_n3aux(nodes, panel) = (nodes[1, panel[2]]-nodes[1, panel[1]])*(nodes[2, panel[3]]-nodes[2, panel[1]]) - (nodes[2, panel[2]]-nodes[2, panel[1]])*(nodes[1, panel[3]]-nodes[1, panel[1]])
_calc_nnorm(nodes, panel) = sqrt(_calc_n1aux(nodes, panel)^2 + _calc_n2aux(nodes, panel)^2 + _calc_n3aux(nodes, panel)^2)
_calc_n1(nodes, panel) = _calc_n1aux(nodes, panel) / _calc_nnorm(nodes, panel)
_calc_n2(nodes, panel) = _calc_n2aux(nodes, panel) / _calc_nnorm(nodes, panel)
_calc_n3(nodes, panel) = _calc_n3aux(nodes, panel) / _calc_nnorm(nodes, panel)

_calc_n1aux(nodes) = (nodes[2][2]-nodes[1][2])*(nodes[3][3]-nodes[1][3]) - (nodes[2][3]-nodes[1][3])*(nodes[3][2]-nodes[1][2])
_calc_n2aux(nodes) = (nodes[2][3]-nodes[1][3])*(nodes[3][1]-nodes[1][1]) - (nodes[2][1]-nodes[1][1])*(nodes[3][3]-nodes[1][3])
_calc_n3aux(nodes) = (nodes[2][1]-nodes[1][1])*(nodes[3][2]-nodes[1][2]) - (nodes[2][2]-nodes[1][2])*(nodes[3][1]-nodes[1][1])
_calc_nnorm(nodes) = sqrt(_calc_n1aux(nodes)^2 + _calc_n2aux(nodes)^2 + _calc_n3aux(nodes)^2)
_calc_n1(nodes) = _calc_n1aux(nodes) / _calc_nnorm(nodes)
_calc_n2(nodes) = _calc_n2aux(nodes) / _calc_nnorm(nodes)
_calc_n3(nodes) = _calc_n3aux(nodes) / _calc_nnorm(nodes)
function _calc_normal(nodes::Array{Arr1,1}) where{Arr1<:AbstractArray}
  # n = cross(nodes[2]-nodes[1], nodes[3]-nodes[1])
  # return n/norm(n)
  return [_calc_n1(nodes), _calc_n2(nodes), _calc_n3(nodes)]
end


_calc_o1aux(nodes, panel) = _calc_n2(nodes, panel)*_calc_t3(nodes, panel) - _calc_n3(nodes, panel)*_calc_t2(nodes, panel)
_calc_o2aux(nodes, panel) = _calc_n3(nodes, panel)*_calc_t1(nodes, panel) - _calc_n1(nodes, panel)*_calc_t3(nodes, panel)
_calc_o3aux(nodes, panel) = _calc_n1(nodes, panel)*_calc_t2(nodes, panel) - _calc_n2(nodes, panel)*_calc_t1(nodes, panel)
_calc_onorm(nodes, panel) = sqrt(_calc_o1aux(nodes, panel)^2 + _calc_o2aux(nodes, panel)^2 + _calc_o3aux(nodes, panel)^2)
_calc_o1(nodes, panel) = _calc_o1aux(nodes, panel) / _calc_onorm(nodes, panel)
_calc_o2(nodes, panel) = _calc_o2aux(nodes, panel) / _calc_onorm(nodes, panel)
_calc_o3(nodes, panel) = _calc_o3aux(nodes, panel) / _calc_onorm(nodes, panel)

_calc_o1aux(nodes) = _calc_n2(nodes)*_calc_t3(nodes) - _calc_n3(nodes)*_calc_t2(nodes)
_calc_o2aux(nodes) = _calc_n3(nodes)*_calc_t1(nodes) - _calc_n1(nodes)*_calc_t3(nodes)
_calc_o3aux(nodes) = _calc_n1(nodes)*_calc_t2(nodes) - _calc_n2(nodes)*_calc_t1(nodes)
_calc_onorm(nodes) = sqrt(_calc_o1aux(nodes)^2 + _calc_o2aux(nodes)^2 + _calc_o3aux(nodes)^2)
_calc_o1(nodes) = _calc_o1aux(nodes) / _calc_onorm(nodes)
_calc_o2(nodes) = _calc_o2aux(nodes) / _calc_onorm(nodes)
_calc_o3(nodes) = _calc_o3aux(nodes) / _calc_onorm(nodes)
function _calc_oblique(nodes::Array{Arr1,1}) where{Arr1<:AbstractArray}
  # t = _calc_tangent(nodes)
  # n = _calc_normal(nodes)
  # return cross(n,t)
  return [_calc_o1(nodes), _calc_o2(nodes), _calc_o3(nodes)]
end

function _calc_unitvectors(nodes::Array{Arr1,1}) where{Arr1<:AbstractArray}
  t = _calc_tangent(nodes)
  n = _calc_normal(nodes)
  # o = cross(n,t)
  o = _calc_oblique(nodes)
  return t, o, n
end
##### END OF TRIANGULAR GRID ###################################################
