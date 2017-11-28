#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

function wing_example(; prompt=true, file_name="temp_wing00")

  # PARAMETERS
  n_up = 40             # Number of sections in upper surface of wing
  n_lower = 20          # Number of sections in lower surface of wing
  z1 = 0.0              # Position of first airfoil
  z2 = 10.0             # Position of second airfoil
  c1 = 2.5              # Chord of first airfoil
  c2 = 1.0              # Chord of second airfoil
  # file_name = "temp_wing00"    # Name for vtk outpout file

  # Reads the original airfoil geometry from airfoiltools.com
  org_x1, org_y1 = readcontour(vtk.data_path*"naca6412.dat"; header_len=1)
  org_x2, org_y2 = readcontour(vtk.data_path*"naca0006.dat"; header_len=1)

  # Separate upper and lower sides to make the contour injective in x
  upper1, lower1 = splitcontour(org_x1, org_y1)
  upper2, lower2 = splitcontour(org_x2, org_y2)

  # Parameterize both sides independently
  fun_upper1 = vtk.parameterize(upper1[1], upper1[2], zeros(upper1[1]); inj_var=1)
  fun_lower1 = vtk.parameterize(lower1[1], lower1[2], zeros(lower1[1]); inj_var=1)
  fun_upper2 = vtk.parameterize(upper2[1], upper2[2], zeros(upper2[1]); inj_var=1)
  fun_lower2 = vtk.parameterize(lower2[1], lower2[2], zeros(lower2[1]); inj_var=1)

  # Upper surface sections
  aux1 = Int(floor(20/53*n_up))
  sec1 = (0.35, aux1, 3.0, true) # 35% of the line has 20 sections in ratio 3.0 around center
  sec2 = (0.65, n_up-aux1, 3.0, true) # 65% of the line has 33 sections in ratio 3.0 around center

  # New discretization for both surfaces
  upper_points1 = vtk.multidiscretize(fun_upper1, 0, 1, [sec1,sec2])
  lower_points1 = vtk.discretize(fun_lower1, 0, 1, n_lower, 8.0; central=true)
  upper_points2 = vtk.multidiscretize(fun_upper2, 0, 1, [sec1,sec2])
  lower_points2 = vtk.discretize(fun_lower2, 0, 1, n_lower, 8.0; central=true)

  # Put both surfaces back together from TE over the top and from LE over the bottom.
  reverse!(upper_points1)                           # Trailing edge over the top
  new_x1 = [point[1] for point in upper_points1]
  new_y1 = [point[2] for point in upper_points1]     # Leading edge over the bottom
  new_x1 = vcat(new_x1, [point[1] for point in lower_points1])
  new_y1 = vcat(new_y1, [point[2] for point in lower_points1])
  reverse!(upper_points2)                           # Trailing edge over the top
  new_x2 = [point[1] for point in upper_points2]
  new_y2 = [point[2] for point in upper_points2]     # Leading edge over the bottom
  new_x2 = vcat(new_x2, [point[1] for point in lower_points2])
  new_y2 = vcat(new_y2, [point[2] for point in lower_points2])

  plot_airfoil(new_x1, new_y1; style="--.k", title_str="NACA6412")
  println("Close figure and press ENTER")
  readline()
  plot_airfoil(new_x2, new_y2; style="--.k", title_str="NACA0006")
  println("Close figure and press ENTER")
  readline()

  npoints = size(new_x1)[1]         # Number of points on each airfoil contour
  # Dimensions the airfoil acording to their chord length and positions them in z
  airfoil1 = [[c1*new_x1[i], c1*new_y1[i], z1] for i in 1:npoints]
  airfoil2 = [[c2*new_x2[i], c2*new_y2[i], z2] for i in 1:npoints]

  # Discretization of the lofting between the airfoils
  nloft = 50                      # Number of sections
  rloft = 10.0                   # Expansion ratio
  central = true                  # Expands about center between edges
  sections = [(1.0, nloft, rloft, central)]   # Discretization sections

  # Dummy point data for good looking visuals
  pd1 = [i for i in 1:size(airfoil1)[1]]
  pd2 = size(airfoil1)[1]+[i for i in 1:size(airfoil2)[1]]

  # Generates cells in VTK Legacy format
  out = vtk.lines2vtkmulticells(airfoil1, airfoil2, sections;
                                          point_data1=pd1, point_data2=pd2)
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
  run(`paraview --data="$(file_name).vtk;"`)


  # Deletes files
  if prompt
   print("Delete vtk files? ([y]/n) ")
   y = readline()
  else
   y = "y"
  end
  if y=="y"; run(`rm -f $(file_name).vtk`); end;
end
