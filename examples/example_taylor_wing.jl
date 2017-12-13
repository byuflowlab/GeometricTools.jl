#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

function taylor_wing(; prompt=true, file_name="temp_taylor00")
  airfoil_folder = "sg6043/"
  airfoil_files = airfoil_folder.*["sg6043-0.140.dat",
                                  "sg6043-0.150.dat",
                                  "sg6043-0.120.dat",
                                  "sg6043-0.120.dat",
                                  "sg6043-0.150.dat",
                                  "sg6043-0.050.dat"]

  println("Defining variables...")
  prev_time = time()

  # GEOMETRY DEFINITION
  b = 60.15/4          # Span
  sweep = pi/180*vcat([23.20 for i in 1:6], []) # Sweep of each section
  dihedral = pi/180*vcat([0 for i in 1:5], [23.20]) # Dihedral of each section

  y_pos = b*[0, 1/42, 1/10, 3/10, 1/2-1/487/60.15, 1/2]   # span position of each section
  x_pos = [y*tan(sweep[i]) for (i,y) in enumerate(y_pos)]
  z_pos = [y*tan(dihedral[i]) for (i,y) in enumerate(y_pos)]

  # PARAMETERS
  n_upper = 20             # Number of sections in upper surface of blade
  n_lower = 20          # Number of sections in lower surface of blade
  r = [1.0, 1.0, 8.0, 8.0, 8.0]           # Expansion ratio in both surfaces of each airfoil
  sections = [[(1.0, 10, 1.0, false)],      # Discretization between each airfoil
              [(1.0, 5, 1.0, false)],
              [(1.0, 15, 1.0, false)],
              [(1.0, 10, 1.0, false)],
              [(1.0, 10, 1.0, false)]
      ]

  chords = [1.619, 1.473, 1.456, 0.789, 0.789, 0.789]    # Chord length of each airfoil

  # Leading edge position of each airfoil
  Os = [ [x_pos[i], y_pos[i], z_pos[i]] for i in 1:size(airfoil_files)[1]]
  # Orientation of chord of each airfoil (yaw, pitch, roll)
  orien = [ [0.0, 0.0, 270.0],
              [-3.32, 0.0, 270.0],
              [-3.90, 0.0, 270.0],
              [-4.35, 0.0, 270.0],
              [-7.32, 0.0, 270.0],
              [-6.28, 0.0, 270.0]
  ]

  crosssections = []        # It will store here the cross sections for lofting
  point_datas = []          # Dummy point data for good looking visuals

  println("\t Runtime: $(round(time()-prev_time,1)) (s)")
  println("Processing airfoils...")
  prev_time = time()

  # Processes each airfoil geometry
  styles = ["--k", "--r", "--g", "--b", "--y", "--c"]
  org_points = []
  for (i,airfoil_file) in enumerate(airfoil_files)

      # Read airfoil file
      x,y = readcontour(vtk.data_path*airfoil_file; header_len=0)
      push!(org_points, [x,y])


      # Separate upper and lower sides to make the contour injective in x
      upper, lower = splitcontour(x, y)

      # Parameterize both sides independently
      fun_upper = vtk.parameterize(upper[1], upper[2], zeros(upper[1]); inj_var=1)
      fun_lower = vtk.parameterize(lower[1], lower[2], zeros(lower[1]); inj_var=1)

      # New discretization for both surfaces
      upper_points = vtk.discretize(fun_upper, 0, 1, n_upper, r[1]; central=true)
      lower_points = vtk.discretize(fun_lower, 0, 1, n_lower, r[1]; central=true)

      # Put both surfaces back together from TE over the top and from LE over the bottom.
      reverse!(upper_points)                           # Trailing edge over the top
      new_x = [point[1] for point in upper_points]
      new_y = [point[2] for point in upper_points]      # Leading edge over the bottom
      new_x = vcat(new_x, [point[1] for point in lower_points])
      new_y = vcat(new_y, [point[2] for point in lower_points])

      plot_airfoil(new_x, new_y; style=styles[i], label=airfoil_file)
  #     println("Close figure and press ENTER")
  #     readline()

      # Scales the airfoil acording to its chord length
      new_x = chords[i]*new_x
      new_y = chords[i]*new_y

      # Reformats into points
      npoints = size(new_x)[1]
      airfoil = Array{Float64, 1}[[new_x[j], new_y[j], 0] for j in 1:npoints]

      # Positions the airfoil along the blade in the right orientation
      Oaxis = vtk.rotation_matrix(orien[i][1], orien[i][2], orien[i][3])
      invOaxis = inv(Oaxis)
      airfoil = vtk.countertransform(airfoil, invOaxis, Os[i])

      push!(crosssections, airfoil)
      push!(point_datas, [j for j in npoints*(i-1)+1:npoints*i])
  end

  println("\t Runtime: $(round(time()-prev_time,1)) (s)")
  println("Lofting and generating VTK files...")
  prev_time = time()

  # Generates cells in VTK Legacy format
  out = vtk.multilines2vtkmulticells(crosssections, sections;
                                        point_datas=point_datas)
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

  println("\t Runtime: $(round(time()-prev_time,1)) (s)")

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
