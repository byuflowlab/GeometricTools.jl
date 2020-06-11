#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################
# ------------ MODULE IMPORTS ------------------------------------------------
import Statistics
using PyPlot

const module_path,_ = splitdir(@__FILE__);      # Path to this module
const data_path = module_path*"/../docs/data/"       # Path to data folder

import GeometricTools
gt = GeometricTools

# ------------ HEADERS -------------------------------------------------------
for header_name in ["simple", "airfoil", "wing", "blade", "taylor_wing"]
include("example_"*header_name*".jl")
end
