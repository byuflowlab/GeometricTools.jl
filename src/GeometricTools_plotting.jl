const plt = PyPlot

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


"Plots the grid on PyPlot"
function plot(grid::Grid; fig_name="gridplot", fontsize=15,
                          xlims=nothing, ylims=nothing, zlims=nothing,
                          labelcells=true, labelnodes=false, labelndivs=true,
                          title_str=nothing,
                          alpha=1.0)

  if grid.dims>3
    error("There is no plotting method for $(grid.dims)-dimensional grids")
  end

  fig = PyPlot.figure(fig_name)
  ax = fig.gca(projection="3d")
  vectors_to_plot = []

  nc = grid.ncells

  # Iterates over every child plotting them
  for i in 1:nc
    nodes = [vcat(get_node(grid, node), zeros(Float64, 3-grid.dims))
                                                  for node in get_cell(grid, i)]

    if labelcells
      center = vcat(get_cellcenter(grid, i), zeros(Float64, 3-grid.dims))
      PyPlot.text3D(center[1], center[2], center[3], "$i", fontsize=fontsize, color="g")
    end


    # Connects the nodes for visualization
    dims = grid.dims - length( [1 for ndiv in grid._ndivscells if ndiv==0] )
    if dims==3
      for j in 1:3
        push!(vectors_to_plot, [nodes[j], -nodes[j]+nodes[j+1]])
        push!(vectors_to_plot, [nodes[4+j], -nodes[4+j]+nodes[4+j+1]])
      end
      push!(vectors_to_plot, [nodes[4], -nodes[4]+nodes[1]])
      push!(vectors_to_plot, [nodes[8], -nodes[8]+nodes[5]])
      for j in 1:4
        push!(vectors_to_plot, [nodes[0+j], -nodes[0+j]+nodes[4+j]])
      end
    elseif dims==2
      for j in 1:3
        push!(vectors_to_plot, [nodes[j], -nodes[j]+nodes[j+1]])
      end
      push!(vectors_to_plot, [nodes[4], -nodes[4]+nodes[1]])
    else
      push!(vectors_to_plot, [nodes[1], -nodes[1]+nodes[2]])
    end
  end

  # Iterates over every dimension labeling the subdivision
  if labelndivs
    for i in 1:grid.dims
      for j in 1:grid._ndivscells[i]
        coor = ones(Int64, grid.dims)
        coor[i] = j
        p1 = get_node(grid, coor)
        if grid.loop_dim==i && j==grid._ndivscells[i]
          coor[i] = 1
        else
          coor[i] += 1
        end
        p2 = get_node(grid, coor)
        center = vcat((p1+p2)/2, zeros(Float64, 3-grid.dims))
        PyPlot.text3D(center[1], center[2], center[3], "$j", fontsize=fontsize,
                                                                      color="r")
      end
    end
  end

  # Plots
  xyzuvw = [[] for _ in 1:6]
  for vector in vectors_to_plot
      _vector = []
      append!(_vector, vector[1])
      append!(_vector, vector[2])
      for (i,val) in enumerate(_vector)
          push!(xyzuvw[i], val)
      end
  end
  x,y,z,u,v,w = [coor for coor in xyzuvw]
  ax.quiver(x,y,z, u,v,w, arrow_length_ratio=0.0, alpha=alpha);

  # Labels nodes
  if labelnodes
    for i in 1:grid.nnodes
      pos = vcat(get_node(grid, i), zeros(Float64, 3-grid.dims))
      PyPlot.text3D(pos[1], pos[2], pos[3], "$i", fontsize=fontsize, color="k")
    end
  end


  # Format axes
  for (lim, lims, label, lbl) in [(PyPlot.xlim, xlims, PyPlot.xlabel, "x"),
                                  (PyPlot.ylim, ylims, PyPlot.ylabel, "y"),
                                  (PyPlot.zlim, zlims, PyPlot.zlabel, "z")]
    if lims!=nothing; lim(lims); end;
    label(lbl)

  end

  # Legend
  handles = []
  for (flag, clr, lbl) in [(labelnodes, "k", "Node"),
                            (labelcells, "g", "Cell"),
                            (labelndivs, "r", "Coordinate")]
    if flag; push!(handles, patch.Patch(color=clr, label=lbl)); end;
  end
  if size(handles,1)>0; PyPlot.legend(handles=handles); end;

  if title_str!=nothing; PyPlot.title(title_str); end;

end
