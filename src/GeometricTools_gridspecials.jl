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
    get_tri_gradient!(out, t1, t2, t3, t0, e1, e2, A, b, res)

Returns gradient of a scalar field with values b, provided at the three vertices.
t1, t2 and t3 are x,y coordinates of the triangle in the local 2D coordinate frame.
t0 is centroid (t1 + t2 + t3)/3
e1 and e2 are basis vectors for transforming from 3D to 2D and vice versa.
A is a 3x3 working matrix and b is 3x1 vector of scalar values
"""
function get_tri_gradient!(grad, t1, t2, t3, t0, e1, e2, A, b)
        A[1, 1] = 1.0
        A[1, 2] = t1[1]-t0[1]
        A[1, 3] = t1[2]-t0[2]

        A[2, 1] = 1.0
        A[2, 2] = t2[1]-t0[1]
        A[2, 3] = t2[2]-t0[2]

        A[3, 1] = 1.0
        A[3, 2] = t3[1]-t0[1]
        A[3, 3] = t3[2]-t0[2]

        grad .= A\b
        # dx = grad[2]
        # dy = grad[3]
end

"""
    get_nodal_data(grid::GridTriangleSurface, field_vals; algorithm=2, areas=nothing)

Converts specified cell-centered field data to node-based data
by averaging field values of cells surrounding the node.
`field_vals` could be an array or vector containing the scalar data values.
Algorithms: 1.Averaging 2.Area-weighted
"""
function get_nodal_data(grid::GridTriangleSurface, field_vals; algorithm=3, areas=nothing, centroids=nothing)

    if length(field_vals) != grid.ncells
        error("No. of field_vals do not match no. of cells")
    end

    # Create an array where each index represents a node
    nodal_data = zeros(grid.nnodes)
    net_weight = zeros(grid.nnodes)

    # Preallocate memory for centroid (if needed)
    C = zeros(grid.orggrid.dims)

    # Parse through each cell and add its field value to that nodal array
    # whose index is the node index. At the end, obtain the weighted average by
    # dividing each element of the nodal array using the net weight
    vtxs = ones(Int, 3)
    for i = 1:grid.ncells

        # Fetch indices of this cell's vertices
        vtxs .= get_cell(grid, i)

        # Case 1: Simple average
        if algorithm == 1
            weight = 1.0

            nodal_data[vtxs] .+= weight * field_vals[i]
            net_weight[vtxs] .+= weight

        # Case 2: Area-weighted average
        elseif algorithm == 2
            weight = isnothing(areas) ? get_area(grid, i) : areas[i]

            nodal_data[vtxs] .+= weight * field_vals[i]
            net_weight[vtxs] .+= weight

        # Case 3: Area/distance^2 weighted average
        elseif algorithm == 3

            # This panel's area
            area = isnothing(areas) ? get_area(grid, i) : areas[i]

            # This panel's centroid
            if isnothing(centroids)

                C .= 0
                for vtx in vtxs
                    for j in 1:grid.orggrid.dims
                        C[j] += grid.orggrid.nodes[j, vtx]
                    end
                end
                C ./= length(vtxs)

            else

                C .= view(centroids, :, i)

            end

            # Iterate over vertices
            for vtx in vtxs

                d2 = 0
                for j in 1:grid.orggrid.dims
                    d2 += (C[j] - grid.orggrid.nodes[j, vtx])^2
                end

                weight = area / d2

                nodal_data[vtx] += weight * field_vals[i]
                net_weight[vtx] += weight
            end

        else

            error("Invalid algorithm $(algorithm); valid values are 1, 2, or 3")

        end

    end

    # Divide each element by number of cells to obtain average
    return nodal_data ./ net_weight
end

function get_nodal_data_TEcells(grid::GridTriangleSurface, field_vals, TE_idx, cells_U, cells_L; algorithm=2, areas=nothing)

    if length(field_vals) != grid.ncells
        error("No. of field_vals do not match no. of cells")
    end

    # Create an array where each index represents a TE node
    nodal_data_U = zeros(length(TE_idx))
    nodal_data_L = zeros(length(TE_idx))
    net_weight_U = zeros(length(TE_idx))
    net_weight_L = zeros(length(TE_idx))

    # Weight can be 1 or area. 1 implies simple averaging.
    weight = ones(grid.ncells)
    if algorithm == 2
        # Compute areas if they are not provided
        weight = isnothing(areas) ? [get_area(grid, i) for i in 1:grid.ncells] : areas
    end

    vtxs = ones(Int, 3)

    # Parse through cells and extract field data for TE nodes
    # Cells on upper side wing/body
    for icell in cells_U
        # If vertex is a TE, add the nodal data to that index
        for vtx in get_cell(grid, icell)
            idx = findfirst(isequal(vtx), TE_idx)
            if idx != nothing
                nodal_data_U[idx] += weight[icell] * field_vals[icell]
                net_weight_U[idx] += weight[icell]
            end
        end
    end

    # Cells on lower side of wing/body
    for icell in cells_L
        # If vertex is a TE, add the nodal data to that index
        for vtx in get_cell(grid, icell)
            idx = findfirst(isequal(vtx), TE_idx)
            if idx != nothing
                nodal_data_L[idx] += weight[icell] * field_vals[icell]
                net_weight_L[idx] += weight[icell]
            end
        end
    end

    # Divide each element by number of cells to obtain average
    return nodal_data_U ./ net_weight_U, nodal_data_L ./ net_weight_L
end

"""
    get_nodal_data(grid::GridTriangleSurface; field_name::String;
                    algorithm=2, areas=nothing)

Converts specified cell-centered field data to node-based data
by averaging field values of cells surrounding the node.
Algorithms: 1.Averaging 2.Area-weighted
"""
function get_nodal_data(grid::GridTriangleSurface, field_name::String, args...; optargs...)
    return get_nodal_data(grid, grid.field[field_name]["field_data"], args...; optargs...)
end

"""
    project_3d_2d!(t2, t3, ex, ey, p1, p2, p3)

Project 3D vertices p1, p2, p3 of a triangle element on to a 2D coordinate system.
Returns coordinates and basis vectors.
First coordinate t1 is always origin [0, 0] and hence not returned. t2 and t3 should be 2-element arrays. Basis vectors ex and ey should be 3-element arrays.
"""
function project_3d_2d!(t2, t3, ex, ey, p1, p2, p3)
    # Compute basis vectors of coordinate system
    # using Gram-Schmidt orthogonalization
    a1 = p2[1] - p1[1]
    a2 = p2[2] - p1[2]
    a3 = p2[3] - p1[3]
    a_mag = sqrt(a1^2 + a2^2 + a3^2)

    ex[1] = a1 / a_mag
    ex[2] = a2 / a_mag
    ex[3] = a3 / a_mag

    a_dot_ex = a1*ex[1] + a2*ex[2] + a3*ex[3]

    b1 = p3[1] - p1[1]
    b2 = p3[2] - p1[2]
    b3 = p3[3] - p1[3]

    b_dot_ex = b1*ex[1] + b2*ex[2] + b3*ex[3]
    ey[1] = b1 - b_dot_ex * ex[1]
    ey[2] = b2 - b_dot_ex * ex[2]
    ey[3] = b3 - b_dot_ex * ex[3]
    ey_mag = sqrt(ey[1]^2 + ey[2]^2 + ey[3]^2)

    ey[1] = ey[1] / ey_mag
    ey[2] = ey[2] / ey_mag
    ey[3] = ey[3] / ey_mag

    b_dot_ey = b1*ey[1] + b2*ey[2] + b3*ey[3]

    # Project vertices onto basis vectors to find 2D coordinates
    # t1 = [0.0, 0.0]
    # t2 = [a_dot_ex, 0.0]
    # t3 = [b_dot_ex, b_dot_ey]

    t2[1] = a_dot_ex
    t2[2] = 0.0

    t3[1] = b_dot_ex
    t3[2] = b_dot_ey

    return
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
