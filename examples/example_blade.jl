#=##############################################################################
# # DESCRIPTION
#     Examples of VTK formatting and usage of VTKtools
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
=###############################################################################

function blade_example(; prompt=true, file_name="temp_blade00")
  airfoil_files = ["Cyl1.txt", "Cyl1.txt", "S815.txt", "S809.txt", "S826.txt"]
  # file_name = "temp_blade00"    # Name for vtk outpout file


  # PARAMETERS
  n_upper = 20             # Number of sections in upper surface of blade
  n_lower = 20          # Number of sections in lower surface of blade
  r = [1.0, 1.0, 8.0, 8.0, 8.0]           # Expansion ratio in both surfaces of each airfoil
  sections = [[(1.0, 10, 1.0, false)],      # Discretization between each airfoil
              [(1.0, 5, 1.0, false)],
              [(1.0, 15, 1.0, false)],
              [(1.0, 10, 1.0, false)],
      ]

  Rtip = 25.0           # Radius at blade tip
  Rhub = 1.0            # Radius of the hub
  pos = [0, 0.15, 0.2, 0.5, 1.0]           # Position along blade of each airfoil
  chords = [1.0, 0.6, 1.75, 3.0, 0.85]    # Chord length of each airfoil

  # Leading edge position of each airfoil
  Os = [[-chords[1]/2, 0, Rhub+pos[1]*(Rtip-Rhub)],
          [-chords[2]*7/8, 0, Rhub+pos[2]*(Rtip-Rhub)],
          [-chords[3]*4/6, 0, Rhub+pos[3]*(Rtip-Rhub)],
          [-chords[4]/4, 0, Rhub+pos[4]*(Rtip-Rhub)],
          [-chords[5]/8, 0, Rhub+pos[5]*(Rtip-Rhub)]
      ]
  # Orientation of chord of each airfoil (yaw, pitch, roll)
  orien = [ [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [10.0, 0.0, 0.0],
              [7.5, 0.0, 0.0],
              [0.0, 30.0, 0.0]
  ]

  crosssections = []        # It will store here the cross sections for lofting
  point_datas = []          # Dummy point data for good looking visuals

  # Processes each airfoil geometry
  styles = ["--k", "--r", "--g", "--b", "--y", "--c"]
  org_points = []
  for (i,airfoil_file) in enumerate(airfoil_files)

      # Read airfoil file
      x,y = readcontour(vtk.data_path*airfoil_file; header_len=2)
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
