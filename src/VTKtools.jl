"""
  Methods for the manipulation of geometric data and formatting into VTK
    format.

  # AUTHORSHIP
    * Author    : Eduardo J Alvarez
    * Email     : Edo.AlvarezR@gmail.com
    * Created   : Aug 2017
    * License   : MIT License

  # DEVELOPMENT HISTORY
    * 20170801  : Write a few simple examples of VTK Legacy formatting.
    * 20171124  : Write methods for lofting and surfacing. Package formatting.
"""
module VTKtools

# ------------ GENERIC MODULES -------------------------------------------------
using PyCall
using Dierckx
using Roots
using QuadGK

# ------------ GLOBAL VARIABLES ------------------------------------------------
global module_path; module_path,_ = splitdir(@__FILE__);   # Path to this module
global data_path = module_path*"/../data/"       # Path to data folder

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["vtk", "geometry", "misc"]
  include("VTKtools_"*header_name*".jl")
end

end # END OF MODULE
