#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

module VTKtools_examples

  # ------------ MODULE IMPORTS ------------------------------------------------
  using PyPlot

  const module_path,_ = splitdir(@__FILE__);      # Path to this module
  const data_path = module_path*"/../docs/data/"       # Path to data folder
  push!(LOAD_PATH, joinpath(module_path,"../src/"))   # Point to VTKtools source

  using VTKtools
  # include("../src/VTKtools.jl")
  vtk = VTKtools

  # ------------ HEADERS -------------------------------------------------------
  for header_name in ["simple", "airfoil", "wing", "blade", "taylor_wing"]
    include("example_"*header_name*".jl")
  end

end
