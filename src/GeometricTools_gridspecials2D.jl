#=##############################################################################
# DESCRIPTION
    Special 2D grid types
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2024
  * License   : MIT License
=###############################################################################

################################################################################
# 2D SURFACE GRID TYPE
################################################################################

"""
  `GridSurface2D(orggrid)`

Receives a 2D surface grid (i.e., a 2D line) and wraps it for automatic
processing accordingly.

  **Arguments**
  * `orggrid`         : Original 2D line grid.
"""
mutable struct GridSurface2D <: AbstractGrid

    # User inputs
    orggrid::Grid                       # Original quadrilateral surface grid

    # Properties
    dims::Int                           # Number of dimensions
    nnodes::Int                         # Number of nodes
    ncells::Int                         # Number of cells
    field::Dict{String, Dict{String, Any}}  # Calculated fields

    # Internal data
    _nodes::Matrix                      # Position of each node
    _ndivsnodes::Tuple                  # Number of nodes in each coordinate
    _ndivscells::Tuple                  # Number of cells in each coordinate
    _override_vtkcelltype::Int          # Option for overriding vtk outputs

    function GridSurface2D(orggrid;
                            dims=orggrid.dims,
                            nnodes=orggrid.nnodes, ncells=orggrid.ncells,
                            field=Dict{String, Dict{String, Any}}(),
                            _nodes=orggrid.nodes,
                            _ndivsnodes=orggrid._ndivsnodes,
                            _ndivscells=orggrid._ndivscells,
                            _override_vtkcelltype=4
                            )
        return new(orggrid,
                    dims,
                    nnodes, ncells, field,
                    _nodes, _ndivsnodes, _ndivscells,
                    _override_vtkcelltype)
    end

end


##### FUNCTIONS  ###############################################################

get_node(self::GridSurface2D, args...; optargs...) = get_node(self.orggrid, args...; optargs...)
get_cell(self::GridSurface2D, args...; optargs...) = get_cell(self.orggrid, args...; optargs...)
get_cell_t!(self::GridSurface2D, args...; optargs...) = get_cell_t!(self.orggrid, args...; optargs...)
generate_getcellt_args!(self::GridSurface2D, args...; optargs...) = get_cell(self.orggrid, args...; optargs...)
lintransform!(self::GridSurface2D, args...; optargs...) = lintransform!(self.orggrid, args...; optargs...)

get_area(self::GridSurface2D, i) = get_area2D(self._nodes, get_cell(self, i))
function _get_area2D(nodes, panel)
    p1, p2 = panel
    return sqrt( (nodes[1, p2] - nodes[1, p1])^2 + (nodes[2, p2] - nodes[2, p1])^2 )
end

##### INTERNAL FUNCTIONS  ######################################################
_calc_t1_2D(args...) = _calc_t1(args...) / _calc_tnorm_2D(args...)
_calc_t2_2D(args...) = _calc_t2(args...) / _calc_tnorm_2D(args...)
_calc_tnorm_2D(nodes::AbstractMatrix, panel) = sqrt(
                                                        (nodes[1, panel[2]] - nodes[1, panel[1]])^2
                                                        + (nodes[2, panel[2]] - nodes[2, panel[1]])^2
                                                    )

_calc_n1aux_2D(nodes::AbstractMatrix, panel) = -(nodes[2, panel[2]] - nodes[2, panel[1]])
_calc_n2aux_2D(nodes::AbstractMatrix, panel) = nodes[1, panel[2]] - nodes[1, panel[1]]
_calc_nnorm_2D(args...) = sqrt(_calc_n1aux_2D(args...)^2 + _calc_n2aux_2D(args...)^2)
_calc_n1_2D(args...) = _calc_n1aux_2D(args...) / _calc_nnorm_2D(args...)
_calc_n2_2D(args...) = _calc_n2aux_2D(args...) / _calc_nnorm_2D(args...)

##### END OF TRIANGULAR GRID ###################################################
