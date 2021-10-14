"""
  Methods for the manipulation of geometric data.

  # AUTHORSHIP
    * Author    : Eduardo J Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Aug 2017
    * License   : MIT License
"""
module GeometricTools

using Printf
using LinearAlgebra
using Requires
using Statistics
import PyCall
import Dierckx
import Roots
import QuadGK

const patch = PyCall.PyNULL()
const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Type of multidiscretize input
const multidisctype = Vector{Tuple{Float64,Int64,Float64,Bool}}

function __init__()
    copy!(patch, PyCall.pyimport_conda("matplotlib.patches", "patches"))
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("GeometricTools_plotting.jl")
end

for header_name in ["vtk", "geometry", "misc", "gridabstract", "airfoil",
                    "surfacing", "plot3d", "vtk_shapes", "conics",
                    "statistics", "linearalgebra",
                    "DEPRECATED"]
    include("GeometricTools_"*header_name*".jl")
end

end # END OF MODULE
