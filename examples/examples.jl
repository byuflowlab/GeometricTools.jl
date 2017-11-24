#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

module_path,_ = splitdir(@__FILE__);   # Path to this module
push!(LOAD_PATH, joinpath(module_path,"../src/"))

import VTKtools
const vtk = VTKtools

"""
  Example of using lines2vtk with generateVTK. Generates a simple 6-sided box.
"""
function main()
  # Defining points
  p0 = [0,0,0]
  p1 = [1,0,0]
  p2 = [1,1,0]
  p3 = [0,1,0]
  p4 = [0,0,1]
  p5 = [1,0,1]
  p6 = [1,1,1]
  p7 = [0,1,1]

  # Defining sides (cells)
  c0 = [p0, p3, p2, p1] # Bottom
  c1 = [p4, p5, p6, p7] # Top
  c2 = [p0, p4, p7, p3] # Left
  c3 = [p1, p2, p6, p5] # Right
  c4 = [p0, p1, p5, p4] # Front
  c5 = [p2, p3, p7, p6] # Back

  # Formats the cells into vtk format
  points, vtk_cells = vtk.lines2vtk([c0,c1,c2,c3,c4,c5])

  # Generates the vtk file
  vtk.generateVTK("vtk_example", points; cells=vtk_cells)

  # Calls paraview
  run(`paraview --data="temp_vtk_example.vtk;"`)


  # ------ Example of outputting a data field ----------------------------------
  # Data for each point in each cell
  d0 = [0, 3, 2, 1] # Data of cell 0
  d1 = [4, 5, 6, 7]
  d2 = [0, 4, 7, 3]
  d3 = [1, 2, 6, 5]
  d4 = [0, 1, 5, 4]
  d5 = [2, 3, 7, 6]

  # Formats the cells into vtk format
  points, vtk_cells, vtk_values = vtk.lines2vtk([c0,c1,c2,c3,c4,c5];
                                             values=[d0,d1,d2,d3,d4,d5])

  # Formats the data for generateVTK
  data = []
  push!(data, Dict(
                  "field_name" => "Point_index",
                  "field_type" => "scalar",
                  "field_data" => vtk_values
                  )
       )
  # Generates the vtk file
  vtk.generateVTK("vtk_example_w_data", points;
               cells=vtk_cells, point_data=data)
  # Calls paraview
  run(`paraview --data="temp_vtk_example_w_data.vtk;"`)

end
