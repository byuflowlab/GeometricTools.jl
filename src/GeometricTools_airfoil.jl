#=##############################################################################
# DESCRIPTION
    Methods for manipulating airfoil geometries.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : MIT License
=###############################################################################



"Receives a .dat file as pulled from airfoiltools.com containing the x and y
contour coordinates of an airfoil, and returns arrays x and y."
function readcontour(file_name; header_len=1, delim=" ", path="")
  x, y = Float64[], Float64[]

  open(joinpath(path,file_name)) do f
    for (i,line) in enumerate(eachline(f))

      # Ignores header
      if i<=header_len
        nothing
      # Parses each line
      else
        this_x, this_y = split(line, delim)
        push!(x, parse(Float64, this_x))
        push!(y, parse(Float64, this_y))
      end

    end
  end
  return x,y
end

"Plots the contour of an airfoil given in x,y"
function plot_airfoil(x::Array{Float64,1}, y::Array{Float64,1};
                      label="Airfoil", style="-k", figfactor=1.0,
                      title_str="Airfoil geometry", legend_loc="best",
                      side_legend=(1.15, 1.1), zoom_factor=1.0, alpha=1.0)
  # Sizes the figure
  figsize = [7*1.5,5*0.5]*figfactor
  xmin, xmax = [-0.05, 1.05]*zoom_factor
  yrange = (xmax-xmin)/figsize[1] * figsize[2]
  ymin, ymax = -yrange/2, yrange/2
  fig1 = plt.figure(title_str, figsize=(figsize[1], figsize[2]))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])

  if label!="Airfoil"
    plt.plot(x,y, style, label=label, alpha=alpha)
    if side_legend!=nothing
      plt.legend(bbox_to_anchor=side_legend)
    else
      plt.legend(loc="best")
    end
  else
    plt.plot(x,y, style)
  end
  plt.xlabel("x")
  plt.ylabel("y")
  plt.grid(true, color="0.8", linestyle="--")
  plt.title(title_str)

end


# """
#   Receives an airfoil contour and splits it up in upper and lower surfaces as
#   divided by the chord line. It returns `(upper, lower)` with `upper=(x,y)` the
#   points of the upper surface, ditto for `lower`. Both `upper` and `lower` are
#   given in increasing order in x (i.e., from leading to trailing edge).
# """
# function splitcontour(x,y)
#   # ERROR CASES
#   if !(x[1] in [0.0, 1.0])
#     error("Invalid contour. x[1] must be either 0.0 or 1.0, got $(x[1]).")
#   end
#
#   # Flag indicating whether the contour start at the trailing or leading edge
#   start_TE = x[1]==1.0
#
#   # Find the opposite end of the contour
#   end_i = -1
#   for (i, xi) in enumerate(x)
#     if i==1
#       nothing
#     # Case of starting from the trailing edge
#     elseif start_TE && xi > x[i-1]
#       end_i = i-1
#       break
#     # Case of leading edge
#     elseif !start_TE  && xi < x[i-1]
#       end_i = i-1
#       break
#     end
#   end
#
#   # ERROR CASE
#   if end_i==-1
#     error("Logic error! End of airfoil not found!")
#   end
#
#   # Splits them up
#   x_sec1, y_sec1 = x[1:end_i], y[1:end_i]
#   x_sec2, y_sec2 = x[end_i:end], y[end_i:end]
#
#   # Sorts them from LE to TE
#   if x_sec1[1] > 0.5; reverse!(x_sec1); reverse!(y_sec1); end;
#   if x_sec2[1] > 0.5; reverse!(x_sec2); reverse!(y_sec2); end;
#
#   # Determines upper and lower surfaces
#   if mean(y_sec1) > mean(y_sec2)
#     upper = [x_sec1, y_sec1]
#     lower = [x_sec2, y_sec2]
#   else
#     upper = [x_sec2, y_sec2]
#     lower = [x_sec1, y_sec1]
#   end
#
#   return upper, lower
# end
