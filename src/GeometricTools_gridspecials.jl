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
  ndivscells[self.dimsplit] /= 2
  get_cell_t!(quad_out, self.orggrid, quadcoor, lin, ndivscells)
  ndivscells[self.dimsplit] *= 2

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


function get_cell_t(tricoor, quadcoor, self::GridTriangleSurface,
                                                i::Int, nodei::Int, lin, ndivscells, cin)
  if i>self.ncells
    error("Requested invalid cell index $i; max is $(self.ncells).")
  end

  for j in 1:length(tricoor); tricoor[j] = cin[i][j]; end;
  return get_cell_t(quadcoor, self, tricoor, nodei, lin, ndivscells)
end

function get_cell_t(quadcoor, self::GridTriangleSurface, coor, nodei::Int,
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


  # Get the node from the quadrilateral panel
  ndivscells[self.dimsplit] /= 2
  if coor[self.dimsplit]%2!=0

    # return [quadnodes[1], quadnodes[2], quadnodes[3]]
    out = get_cell_t(self.orggrid, quadcoor, nodei, lin, ndivscells)

  else

    # return [quadnodes[3], quadnodes[4], quadnodes[1]]
    out =   nodei == 1 ? get_cell_t(self.orggrid, quadcoor, 3, lin, ndivscells) :
            nodei == 2 ? get_cell_t(self.orggrid, quadcoor, 4, lin, ndivscells) :
                         get_cell_t(self.orggrid, quadcoor, 1, lin, ndivscells)
  end
  ndivscells[self.dimsplit] *= 2

  return out
end

function generate_getcellt_args!(grid)
    # Pre-allocate memory for panel calculation
    lin = LinearIndices(grid._ndivsnodes)
    ndivscells = vcat(grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells)))
    tri_out = zeros(Int, 3)
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)
    quad_out = zeros(Int, 4)

    return tri_out, tricoor, quadcoor, quad_out, lin, ndivscells, cin
end

function generate_getcellt_args(grid)
    # Pre-allocate memory for panel calculation
    lin = LinearIndices(grid._ndivsnodes)
    ndivscells = vcat(grid._ndivscells...)
    cin = CartesianIndices(Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells)))
    tricoor = zeros(Int, 3)
    quadcoor = zeros(Int, 3)

    return tricoor, quadcoor, lin, ndivscells, cin
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





"""
    neighbor!(ncoor::Vector, ni::Int, ci::Int, cin, ccoor,
                                    ndivscells, dimsplit::Int)

Returns the Cartesian coordinates of the `ni`-th neighboor of the `ci`-th cell,
storing the coordinates under `ncoor`. `cin` are the CartesianIndices of a
triangular grid of dimensions `ndivscells` split along dimension `dimsplit`.

> **NOTE:** This assumes that the grid is periodic and closed, meaning that the
neighbors of the starting cells include the end cells and vice versa.

```@example

# Pre-calculations
ndivscells = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
cin = CartesianIndices(ndivscells)
dimsplit = grid.dimsplit

# Identify neighbor
ncoor = ones(Int, 3)               # Store neighboor coordinates here
ni = 3                              # Request the third neighbor
ci = 10                             #   of cell number 10
ccoor = cin[ci]                     # Cartesian indexing of this cell
neighbor!(ncoor, ni, ci, ccoor, ndivscells, dimsplit)

# Fetch normal of such neighbor
normal = get_normal(grid, ncoor)

# Convert Cartesian coordinates to linear indexing
lin = LinearIndices(ndivscells)
ni = lin[ncoor...]
```
"""
function neighbor!(ncoor::AbstractVector{<:Int}, ni::Int, ci::Int,
                    ccoor, ndivscells, dimsplit::Int)

    @assert 1 <= ni <= 3 "Invalid neighbor $(ni); cell has only 3 neighbors."

    if dimsplit==1

        if ci%2==1          # Case: odd cell

            # Neighbor between node 1 and 2
            if ni==1
                ncoor[1] = ccoor[1]+1
                ncoor[2] = ccoor[2] != 1 ? ccoor[2]-1 : ndivscells[2][end]

            # Neighbor between node 2 and 3
            elseif ni==2
                ncoor[1] = ccoor[1]+1 != ndivscells[1][end] ? ccoor[1]+3 : 2
                ncoor[2] = ccoor[2]

            # Neighbor between node 3 and 1
            else
                ncoor[1] = ccoor[1]+1
                ncoor[2] = ccoor[2]
            end

        else                # Case: even cell

            # Neighbor between node 1 and 2
            if ni==1
                ncoor[1] = ccoor[1]-1
                ncoor[2] = ccoor[2] != ndivscells[2][end] ? ccoor[2]+1 : 1

            # Neighbor between node 2 and 3
            elseif ni==2
                ncoor[1] = ccoor[1]-1 != 1 ? ccoor[1]-3 : ndivscells[1][end]-1
                ncoor[2] = ccoor[2]

            # Neighbor between node 3 and 1
            else
                ncoor[1] = ccoor[1]-1
                ncoor[2] = ccoor[2]
            end

        end

    # TODO: Implement dimsplit==2
    # elseif dimsplit==2

    else
        error("Neighborhood with splitting dimension $(dimsplit)"*
                " not implemented yet")
    end

    return ncoor
end

function neighbor(ni::Int, ci::Int, ccoor, ndivscells, dimsplit::Int)
    return neighbor!(ones(Int, 3), ni, ci, ccoor, ndivscells, dimsplit)
end

function neighbor(grid::GridTriangleSurface, ni::Int, ci::Int;
                    preserveEdge::Bool=false)
    # Preserve edge will output [0,0,0] for a non-existent neighbor cell
    # This happens for cells at the edges of the grid

    # Pre-calculations
    ndivscells = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
    cin = CartesianIndices(ndivscells)
    ccoor = cin[ci]

    # Calculate neighbor
    neigh = neighbor(ni, ci, ccoor, ndivscells, grid.dimsplit)

    if preserveEdge
        isFakeNeighbor = false
        inXmin = (isedge(grid, ci; whichedge=1) && ni==2)
        inXmax = (isedge(grid, ci; whichedge=2) && ni==2)
        inYmin = (isedge(grid, ci; whichedge=3) && ni==1)
        inYmax = (isedge(grid, ci; whichedge=4) && ni==1)

        # This has only been tested for dim_split = 1
        if grid.orggrid.loop_dim == 2
            isFakeNeighbor = inXmin || inXmax
        elseif grid.orggrid.loop_dim == 1
            isFakeNeighbor = inYmin || inYmax
        else
            isFakeNeighbor = inXmin || inXmax || inYmin || inYmax
        end
        if isFakeNeighbor; neigh = [0,0,0]; end
    end

    return neigh
end

function neighbor(grid::GridTriangleSurface, ni::Int,
                    ccoor::Union{<:AbstractVector, <:Tuple};
                    preserveEdge::Bool=false)

    # Pre-calculations
    ndivscells = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
    lin = LinearIndices(ndivscells)
    ci = lin[ccoor...]

    # Calculate neighbor
    return neighbor(grid, ni, ci; preserveEdge=preserveEdge)
end

function neighbor(grid::GridTriangleSurface, ni::Int, ccoor::CartesianIndex;
                    preserveEdge::Bool=false)

    # Pre-calculations
    ndivscells = Tuple(collect( 1:(d != 0 ? d : 1) for d in grid._ndivscells))
    lin = LinearIndices(ndivscells)
    ci = lin[ccoor]

    # Calculate neighbor
    return neighbor(grid, ni, ci; preserveEdge=preserveEdge)
end

function isedge(grid::GridTriangleSurface, ci::Int; whichedge::Int=0)
    # whichedge takes values 1,2,3,4 or 0 (default)
    # [1, 2, 3, 4] = [Xmin, Xmax, Ymin, Ymax]
    # If whichedge is used, it checks only that specific boundary


    ret = false
    nx = grid.orggrid.NDIVS[1]
    ny = grid.orggrid.NDIVS[2]

    # Determine if the cell is part of any boundary
    # Xmin, Xmax, Ymin, Ymax are the cell boundaries of the grid
    if grid.dimsplit == 2
        inXmin = (ci-nx-1)%(2*nx) == 0
        inXmax = (ci-nx)%(2*nx) == 0
        inYmin = ci <= nx
        inYmax = ci > 2*nx*ny - nx
    else
        inXmin = (ci-2)%(2*nx) == 0
        inXmax = (ci-(2*nx-1))%(2*nx) == 0
        inYmin = ((ci+1)%2 == 0) && (ci < 2*nx)
        inYmax = ((2*nx*ny-ci)%2 == 0) && (ci > 2*nx*(ny-1)+1)
    end

    if whichedge == 1; ret = inXmin
    elseif whichedge == 2; ret = inXmax
    elseif whichedge == 3; ret = inYmin
    elseif whichedge == 4; ret = inYmax
    else
        if grid.orggrid.loop_dim == 2
            ret = inXmin || inXmax

        elseif grid.orggrid.loop_dim == 1
            ret = inYmin || inYmax

        else  # For loop_dim = 0 and loop_dim > 2
            ret = inXmin || inXmax || inYmin || inYmax

        end
    end

    return ret
end

"""
    get_num_cells_around_node(grid::GridTriangleSurface, ci::CartesianIndex)

Returns number of cells that share a specific node in GridTriangleSurface. The node is input using its CartesianIndex.
This function is useful in the averagine process required for converting a grid with cell-center based values to node-based.
"""
function get_num_cells_around_node(grid::GridTriangleSurface, ci::CartesianIndex)

    # Works for all dim_split values

    nx, ny, nz = get_ndivsnodes(grid)
    loop_dim = grid.orggrid.loop_dim

    # Identify whether the given node lies on any of the boundaries
    isXmin = ci[1] == 1
    isXmax = ci[1] == nx
    isYmin = ci[2] == 1
    isYmax = ci[2] == ny

    # If the grid is looped, some boundaries are not present
    if loop_dim == 1
        isXmin = false
        isXmax = false
    end
    if loop_dim == 2
        isYmin = false
        isYmax = false
    end

    # There are four posibilities for number of cells around a node
    # Corner-1 or 2, Edge-3, Interior-6
    isCorner1 = false
    isCorner2 = false
    isEdge = false

    # Check if node is a part of any boundary
    if (isXmin && isYmin) || (isXmax && isYmax)
        isCorner2 = true
    elseif (isXmin && isYmax) || (isXmax && isYmin)
        isCorner1 = true
    elseif (isXmin || isXmax || isYmin || isYmax)
        isEdge = true
    end

    # Assume the node is an interior node unless otherwise
    num_cells = 6

    if isEdge
        num_cells = 3
    elseif isCorner1
        num_cells = 1
    elseif isCorner2
        num_cells = 2
    end

    return num_cells
end

"""
    get_nodal_data(grid::GridTriangleSurface, field_name::String)

Converts specified cell-centered field data to node-based data
by averaging field values of cells surrounding the node.
"""
function get_nodal_data(grid::GridTriangleSurface, field_name::String)
    nx, ny, nz = get_ndivsnodes(grid)

    # Create an array where each index represents a node
    nodal_data = zeros(grid.nnodes)

    # Parse through each cell and add its field value to that nodal array
    # whose index is the node index. At the end, obtain the average by
    # dividing each element of the nodal array using
    # the number of cells around each node
    vtxs = ones(Int, 3)
    field_val = 0.0
    for i = 1:grid.ncells
        vtxs .= get_cell(grid, i)
        scalar_val = grid.field[field_name]["field_data"][i]
        nodal_data[vtxs[1]] += scalar_val
        nodal_data[vtxs[2]] += scalar_val
        nodal_data[vtxs[3]] += scalar_val
    end

    # Divide each element by number of cells to obtain average
    cart = CartesianIndices((1:nx, 1:ny))
    for i = 1:grid.nnodes
        nodal_data[i] /= get_num_cells_around_node(grid, cart[i])
    end

    return nodal_data
end

"""
    project_3d_2d(p1, p2, p3)

Project 3D vertices of triangle element on to a 2D coordinate system.
Returns coordinates and basis vectors. First coordinate is always origin [0, 0] and hence not returned.
"""
function project_3d_2d(p1, p2, p3)
    # Compute basis vectors of coordinate system
    # using Gram-Schmidt orthogonalization
    a = p2 .- p1
    e1 = a ./ norm(a)

    b = p3 .- p1
    b_dot_e1 = dot(b, e1)
    e2 = b .- b_dot_e1 .* e1
    e2 = e2 ./ norm(e2)

    # Project vertices onto basis vectors to find 2D coordinates
    # t1 = [0.0, 0.0]
    t2 = [dot(a, e1), 0.0]
    t3 = [b_dot_e1, dot(b, e2)]

    return t2, t3, e1, e2
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
