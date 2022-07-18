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

"Plots loft splines"
function plot_loft_splines(b_pos, LE_x, LE_z, new_LE_x, new_LE_z, twists, new_twists,
  tilt_z, new_tilt_z, chords, new_chords)

     fig = plt.figure("spl_verif", figsize=(7*2,5*1))

     plt.subplot(121)
     plt.plot(LE_x[:,1], LE_x[:,2], "og", label="Org LE x", alpha=0.5)
     plt.plot(LE_z[:,1], LE_z[:,2], "ob", label="Org LE z", alpha=0.5)
     plt.plot(b_pos, new_LE_x, "--g", label="Spline LE x")
     plt.plot(b_pos, new_LE_z, "--b", label="Spline LE z")
     plt.xlabel(plt.L"y/b_{scale}")
     plt.ylabel(plt.L"x/b_{scale}, z/b_{scale}")
     plt.grid(true, color="0.8", linestyle="--")
     plt.legend(loc="best")

     plt.subplot(122)
     p1 = plt.plot(twists[:,1], twists[:,2], "og", label="Org Twist", alpha=0.5)
     p2 = plt.plot(b_pos, new_twists, "--g", label="Spline twist")
     pextra = []
     if tilt_z!=nothing
       pextra1 = plt.plot(tilt_z[:,1], tilt_z[:,2], "or", label="Org tilt z",
                                                                       alpha=0.5)
       pextra2 = plt.plot(b_pos, new_tlt_z, "--r", label="Spline tilt z")
       pextra = vcat(pextra, [pextra1[1], pextra2[1]])
     end
     plt.ylabel("Twist (deg)")

     plt.grid(true, color="0.8", linestyle="--")
     plt.xlabel(plt.L"y/b_{scale}")

     plt.twinx()
     p3 = plt.plot(chords[:,1], chords[:,2], "ob", label="Org Chord", alpha=0.5)
     p4 = plt.plot(b_pos, new_chords, "--b", label="Spline chord")
     plt.ylabel(plt.L"c/b_{scale}")

     ps = vcat([p1[1], p2[1], p3[1], p4[1]], pextra)
     plt.legend(ps, [p.get_label() for p in ps], loc="best")
     return nothing
end
