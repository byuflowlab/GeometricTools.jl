#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

include("../src/VTKtools.jl")
vtk = VTKtools

"""
  Example of using lines2vtk with generateVTK. Generates a simple 6-sided box.
"""
function simple_box(; prompt=true)
  file_name1 = "temp_vtk_example"
  file_name2 = "temp_vtk_example_w_data"
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
  vtk.generateVTK(file_name1, points; cells=vtk_cells)

  # Calls paraview
  run(`paraview --data="$(file_name1).vtk;"`)


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
  vtk.generateVTK(file_name2, points;
               cells=vtk_cells, point_data=data)
  # Calls paraview
  run(`paraview --data="$(file_name2).vtk;"`)

  # Deletes files
  if prompt
    print("Delete vtk files? ([y]/n) ")
    y = readline()
  else
    y = "y"
  end
  if y=="y"; run(`rm -f $(file_name1).vtk $(file_name2).vtk`); end;

end

"""
  Example of using a parametrically meshed surface using `discretize()`,
  `lines2vtkcells()`, and `generateVTK()`.
"""
function parametric_mesh(; prompt=true)
  file_name = "temp_mesh_example"

  # Helicoid function
  revs = 2                        # Number of revolutions
  rho1 = 8                        # Outer radius
  rho2 = rho1/4                   # Inner radius
  k = 1                           # Period
  f1(z) = [rho1*cos(k*z), rho1*sin(k*z), z] # Outer edge
  f2(z) = [rho2*cos(k*z), rho2*sin(k*z), z] # Inner edge

  # Discretize edges
  xlow, xhigh = 0, revs*2*pi      # Bounds of parametric edges
  n = 100                         # Number of cells
  r = 1.0                         # Expansion ratio
  line1 = vtk.discretize(f1, xlow, xhigh, n, r)
  line2 = vtk.discretize(f2, xlow, xhigh, n, r)

  # Dummy point data for good looking visuals
  pd1 = [i for i in 1:size(line1)[1]]
  pd2 = size(line1)[1]+[i for i in 1:size(line2)[1]]

  # Generates cells in VTK Legacy format
  out = vtk.lines2vtkcells(line1, line2; point_data1=pd1, point_data2=pd2)
  points, vtk_cells, point_data = out


  # Formats the point data for generateVTK
  data = []
  push!(data, Dict(
                  "field_name" => "Point_index",
                  "field_type" => "scalar",
                  "field_data" => point_data
                  )
       )


   # Generates the vtk file
   vtk.generateVTK(file_name, points; cells=vtk_cells, point_data=data)

   # Calls paraview
  #  run(`paraview --data="$(file_name).vtk;"`)

   # Deletes files
   if prompt
     print("Delete vtk files? ([y]/n) ")
     y = readline()
   else
     y = "y"
   end
   if y=="y"; run(`rm -f $(file_name).vtk`); end;
end


"""
  Example of using a parametrically meshed surface using `discretize()`,
  `lines2vtkmulticells()`, and `generateVTK()`.
"""
function parametric_mesh2(; prompt=true)
  file_name = "temp_mesh_example2"

  # Helicoid function
  revs = 2                        # Number of revolutions
  rho1 = 8                        # Outer radius
  rho2 = rho1/4                   # Inner radius
  k = 1                           # Period
  f1(z) = [rho1*cos(k*z), rho1*sin(k*z), z] # Outer edge
  f2(z) = [rho2*cos(k*z), rho2*sin(k*z), z] # Inner edge

  # Discretize edges
  xlow, xhigh = 0, revs*2*pi      # Bounds of parametric edges
  n = 100                         # Number of cells
  r = 1.0                         # Expansion ratio
  line1 = vtk.discretize(f1, xlow, xhigh, n, r)
  line2 = vtk.discretize(f2, xlow, xhigh, n, r)

  # Discretize width between lines
  nwidth = 25                     # Number of rows between edges
  rwidth = 10.0                   # Expansion ratio
  central = true                  # Expands about center between edges
  sections = [(1.0, nwidth, rwidth, central)]   # Discretization sections

  # Dummy point data for good looking visuals
  pd1 = [i for i in 1:size(line1)[1]]
  pd2 = size(line1)[1]+[i for i in 1:size(line2)[1]]

  # Generates cells in VTK Legacy format
  out = vtk.lines2vtkmulticells(line1, line2, sections;
                                            point_data1=pd1, point_data2=pd2)
  points, vtk_cells, point_data = out

  println(point_data)


  # Formats the point data for generateVTK
  data = []
  push!(data, Dict(
                  "field_name" => "Point_index",
                  "field_type" => "scalar",
                  "field_data" => point_data
                  )
       )


   # Generates the vtk file
   vtk.generateVTK(file_name, points; cells=vtk_cells, point_data=data)

   # Calls paraview
  #  run(`paraview --data="$(file_name).vtk;"`)

   # Deletes files
   if prompt
     print("Delete vtk files? ([y]/n) ")
     y = readline()
   else
     y = "y"
   end
   if y=="y"; run(`rm -f $(file_name).vtk`); end;
end












#
