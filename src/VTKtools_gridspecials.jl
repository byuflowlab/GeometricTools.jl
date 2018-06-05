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
# TRIANGULAR GRID TYPE
################################################################################
"""
  Implementations of AbstractGrid are expected to have the following fields:

  * `dims::Int64`               : Number of dimensions.
  * `nnodes::Int64`             : Number of nodes in the grid.
  * `ncells::Int64`             : Number of cells in the grid.
  <!-- * `bbox::Array{Int64, 1}`     : Bounding box of the grid (cells in each dim). -->

  and the following functions
  ```julia
  function get_node(self::MyGrid, i::Int64)
    # Returns the position of the i-th node (1-indexed) in the grid
  end
  function get_node(self::MyGrid, coor::Array{Int64,1})
    # Returns the position of the node of subscript coordinates `coor`
    # (1-indexed)
  end
  function get_cell(self::MyGrid, i::Int64)
    # Returns the nodes indices of i-th cell in the grid (1-indexed)
  end
  function get_cell(self::MyGrid, coor::Array{Int64,1})
    # Returns the node indices of the cell with subscript coordinates `coor`
    # (1-indexed). The format corresponds to VTK_HEXAHEDRON (=12) in 3D,
    # VTK_QUAD (=9) in 2D, or VTK_LINE (=3) in 1D---except that points are
    # 1-indexed instead of 0-indexed.
  end
  function get_fieldval(self::MyGrid, field_name::String, coor::Array{Int64,1})
    # Returns the value of node of coordinates `coor` (1-indexed) in the field
    # 'field_name'.
  end
  function get_fieldval(self::MyGrid, field_name::String, i::Int64)
    # Returns the value of i-th node (1-indexed) in the field 'field_name'.
  end
  function add_field(self::MyGrid, field_name::String, field_type::String,
                                                                  field_data)
    # Adds a field of data associated to each node.
    #
    # NOTE: each data entry must be a single value if `field_type==scalar`, or a
    #       3-element array if `field_type==vector`.
  end
  function calculate_field(self::MyGrid, f, field_name::String, field_type::String)
    # Evaluates the function `f` at each nodes and stores the values as a new
    # field.
    #
    # NOTE: f must return a single value if `field_type==scalar`, or a 3-element
    #       array if `field_type==vector`.
  end
  function lintransform!(self::MyGrid, M::Array{Float64,2}, T::Array{Float64,1};
                          reset_fields::Bool=true)
    # Rotates and translates the grid by the rotation matrix `M` and translation
    # vector `T` (linear transformation).
  end
  function transform!(self::MyGrid, f; reset_fields::Bool=true)
    # Applies the space transformation given by function `f` to the grid.
  end
  function save(self::MyGrid, filename::String; args...)
    # Outputs a vtk file of this grid
  end
  function plot(self::MyGrid; fig_name="gridplot", fontsize=15,
                            xlims=nothing, ylims=nothing, zlims=nothing,
                            labelcells=true, labelnodes=false, labelndivs=true,
                            title_str=nothing)
    # Plots the grid on PyPlot
  end
  ```
"""
abstract type AbstractGrid end


for header_name in ["grid", "multigrid"]
  include("VTKtools_"*header_name*".jl")
end

# Implementations of AbstractGrid
GridTypes = Union{Grid, MultiGrid}



##### GENERIC FUNCTIONS  #######################################################
"Returns the centroid of the cell"
function get_cellcenter(self::GridTypes, args...)
  nodes = get_cell(self, args...)
  C = sum([get_node(self, node) for node in nodes])/size(nodes, 1)
  return C
end


##### END OF ABSTRACT GRID #####################################################
