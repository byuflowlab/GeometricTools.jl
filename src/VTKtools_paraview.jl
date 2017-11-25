################################################################################
# # DESCRIPTION
#     Wrapper for Python Paraview functions.
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
################################################################################

# Wrap VTKtools_paraview.py
@pyimport imp
(file, filename, data) = imp.find_module("VTKtools_paraview", [module_path])
parapy = imp.load_module("VTKtools_paraview", file, filename, data)
