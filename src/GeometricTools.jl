"""
  Methods for the manipulation of geometric data.

  # AUTHORSHIP
    * Author    : Eduardo J Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Aug 2017
    * License   : MIT License
"""
module GeometricTools

# Trick for avoiding OpenSSL compatibility issues with PyCall
# https://discourse.julialang.org/t/error-loading-openssl-jll-version-openssl-3-3-0-not-found-when-precompiling-mldatasets/128977/3
using OpenSSL_jll

using Printf
using LinearAlgebra
using Requires
using Statistics
import Dierckx
import FLOWMath
import Roots
import QuadGK
import HDF5

import Meshes

import LaTeXStrings: @L_str

import Requires: @require


const module_path = splitdir(@__FILE__)[1]      # Path to this module
                                                # Type of multidiscretize input
const multidisctype = Vector{Tuple{Float64,Int64,Float64,Bool}}

for header_name in ["vtk", "geometry", "misc", "gridabstract", "airfoil",
                    "discretization",
                    "surfacing", "plot3d", "vtk_shapes", "conics",
                    "statistics", "linearalgebra", "xdmf",
                    "DEPRECATED", 
                    "meshesjl"]
    include("GeometricTools_"*header_name*".jl")
end



function __init__()

    # Conditionally load monitors if PyPlot is available
    try
        @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin

            import .PyPlot
            const plt = PyPlot

          # import PyCall
          # @PyCall.pyimport matplotlib.patches as patch

            for header_name in ["plotting"]
              include("GeometricTools_"*header_name*".jl")
            end

        end

    catch e
        @warn "PyPlot is not available; PyPlot (alias `plt`) will not be loaded"
    end

end

end # END OF MODULE
