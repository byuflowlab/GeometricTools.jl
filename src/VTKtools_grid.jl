#=##############################################################################
# DESCRIPTION
    Methods for generation and evaluation of grids.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Feb 2018
  * License   : MIT License
=###############################################################################


"""
  `Grid(P_min, P_max, NDIVS)`

Generates an n-dimensional grid.

  **Arguments**
  * `P_min::Array{Float64,1}`   : Minimum point of the domain.
  * `P_max::Array{Float64,1}`   : Maximum point of the domain.
  * `NDIVS::Array{Int64,1}`     : Number of divisions in each coordinate.

  **Properties**
  * `sol` : Contains calculated fields formated as sol[field_name] = Dict(
                            "field_name" => field_name::String,
                            "field_type" => "scalar" or "vector",
                            "field_data" => data
                            )
            where `data` is an array data[i] = [val1, val2, ...] containing
            this field values (scalar or vector) at each node in the grid.
"""
type Grid

  # User inputs
  P_min::Array{T,1} where {T<:Real}   # Minimum point of the domain
  P_max::Array{T,1} where {T<:Real}   # Maximum point of the domain
  NDIVS::Array{T,1} where {T<:Signed} # Number of divisions in each coordinate

  # Properties
  dims::Int64                         # Number of dimensions
  nnodes::Int64                       # Number of nodes
  sol::Array{Dict{String, Any}, 1}    # Calculated fields

  # Error cases
  if size(P_min)!=size(P_max)
    error("`P_min` and `P_max` must have the same dimensions "*
                                          "($(size(P_min))!=$(size(P_max)))")
  elseif length(P_min)!=length(NDIVS)
    error("Division for each dimension must be given "*
                                          "$(length(P_min))!=$(length(NDIVS))")
  end
  for i in 1:length(P_min)
    if P_min[i]>=P_max[i]
      error("Elements in `P_min` must be lower than elements in `P_max` "*"
                                          ($(P_min[i])>=$(P_max[i]))")
    end
  end
  for this_div in NDIVS
    if this_div<=0
      error("Invalid division $this_div$")
    end
  end

  Grid(P_min, P_max, NDIVS,
              dims=length(P_min), nnodes=prod(NDIVS), sol=Dict{String, Any}[]
      ) = new(P_min, P_max, NDIVS,
              dims, nnodes, sol
      )
end


# "Returns the position of the i-th node in the grid. Nodes are indexed starting on 1."
# function get_X(self::Grid, i::Int64)
#   if i>
#   for
# end
