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

using PyPlot

function airfoil_example()
  x,y = readcontour(vtk.data_path*"S809.txt"; header_len=2)
  plot_airfoil(x,y; style="--.k")

  println("Close figure and press ENTER")
  readline()

  # Separate upper and lower surfaces to make the contour injective in x
  upper, lower = splitcontour(x,y)

  # Parameterize both surfaces independently
  fun_upper = vtk.parameterize(upper[1], upper[2], zeros(upper[1]); inj_var=1)
  fun_lower = vtk.parameterize(lower[1], lower[2], zeros(lower[1]); inj_var=1)

  # New discretization for both surfaces
  upper_points = vtk.discretize(fun_upper, 0, 1, 50, 1.0)
  lower_points = vtk.discretize(fun_lower, 0, 1, 50, 1.0)

  # Put both surfaces back together from TE over the top and from LE over the bottom.
  reverse!(upper_points)                           # Trailing edge over the top
  new_x = [point[1] for point in upper_points]
  new_y = [point[2] for point in upper_points]     # Leading edge over the bottom
  new_x = vcat(new_x, [point[1] for point in lower_points])
  new_y = vcat(new_y, [point[2] for point in lower_points])

  plot_airfoil(new_x, new_y; style="^r", label="Parameterized")
  plot_airfoil(x,y; style="--.k", label="Original")
  legend(loc="best")

  println("Close figure and press ENTER")
  readline()

  # Upper surface sections
  sec1 = (0.35, 20, 3.0, true) # 35% of the line has 20 sections in ratio 3.0 around center
  sec2 = (0.65, 33, 3.0, true) # 65% of the line has 33 sections in ratio 3.0 around center

  # New discretization for both surfaces
  upper_points = vtk.multidiscretize(fun_upper, 0, 1, [sec1,sec2])
  lower_points = vtk.discretize(fun_lower, 0, 1, 50, 8.0; central=true)

  # Put both surfaces back together from TE over the top and from LE over the bottom.
  reverse!(upper_points)                           # Trailing edge over the top
  new_x = [point[1] for point in upper_points]
  new_y = [point[2] for point in upper_points]     # Leading edge over the bottom
  new_x = vcat(new_x, [point[1] for point in lower_points])
  new_y = vcat(new_y, [point[2] for point in lower_points])

  plot_airfoil(new_x, new_y; style="--^r", label="Parameterized")
end

### INTERNAL FUNCTIONS #########################################################
"Plots the contour of an airfoil given in x,y"
function plot_airfoil(x::Array{Float64,1}, y::Array{Float64,1};
                      label="Airfoil", style="-k", figfactor=1.0)
  # Sizes the figure
  figsize = [7*1.5,5*0.5]*figfactor
  xmin, xmax = -0.05, 1.05
  yrange = (xmax-xmin)/figsize[1] * figsize[2]
  ymin, ymax = -yrange/2, yrange/2
  fig1 = figure("airfoil_geometry", figsize=(figsize[1], figsize[2]))
  xlim([xmin, xmax])
  ylim([ymin, ymax])

  PyPlot.plot(x,y, style, label=label)
  xlabel("x")
  ylabel("y")
  grid(true, color="0.8", linestyle="--")
  title("Airfoil geometry")
end

"Receives a .dat file as pulled from airfoiltools.com containing the x and y
contour coordinates of an airfoil, and returns arrays x and y."
function readcontour(file_name; header_len=1)
  x, y = Float64[], Float64[]

  open(file_name) do f
    for (i,line) in enumerate(eachline(f))

      # Ignores header
      if i<=header_len
        nothing
      # Parses each line
      else
        this_x, this_y = split(line)
        push!(x, parse(Float64, this_x))
        push!(y, parse(Float64, this_y))
      end

    end
  end
  return x,y
end


"""
  Receives an airfoil contour and splits it up in upper and lower surfaces as
  divided by the chord line. It returns `(upper, lower)` with `upper=(x,y)` the
  points of the upper surface, ditto for `lower`. Both `upper` and `lower` are
  given in increasing order in x (i.e., from leading to trailing edge).
"""
function splitcontour(x,y)
  # ERROR CASES
  if !(x[1] in [0.0, 1.0])
    error("Invalid contour. x[1] must be either 0.0 or 1.0, got $(x[1]).")
  end

  # Flag indicating whether the contour start at the trailing or leading edge
  start_TE = x[1]==1.0

  # Find the opposite end of the contour
  end_i = -1
  for (i, xi) in enumerate(x)
    if i==1
      nothing
    # Case of starting from the trailing edge
    elseif start_TE && xi > x[i-1]
      end_i = i-1
      break
    # Case of leading edge
    elseif !start_TE  && xi < x[i-1]
      end_i = i-1
      break
    end
  end

  # ERROR CASE
  if end_i==-1
    error("Logic error! End of airfoil not found!")
  end

  # Splits them up
  x_sec1, y_sec1 = x[1:end_i], y[1:end_i]
  x_sec2, y_sec2 = x[end_i:end], y[end_i:end]

  # Sorts them from LE to TE
  if x_sec1[1] > 0.5; reverse!(x_sec1); reverse!(y_sec1); end;
  if x_sec2[1] > 0.5; reverse!(x_sec2); reverse!(y_sec2); end;

  # Determines upper and lower surfaces
  if mean(y_sec1) > mean(y_sec2)
    upper = [x_sec1, y_sec1]
    lower = [x_sec2, y_sec2]
  else
    upper = [x_sec2, y_sec2]
    lower = [x_sec1, y_sec1]
  end

  return upper, lower
end
################################################################################
