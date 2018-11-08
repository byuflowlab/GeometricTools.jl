#=##############################################################################
# DESCRIPTION
    Methods for generation of surface geometries.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2018
  * License   : MIT License
=###############################################################################



# LOFTING ######################################################################
"""
  `generate_loft(crosssections, bscale, b_pos, chords, twists, LE_x, LE_z)`

  Generates a lofted surface geometry. The name of the arguments are taken from
  its initial application to the lofting of an aircraft wing, so bear with it.

  **Arguments**
  * `bscale::Float64`         : Semi-span scale. It will scale the entire
                                geometry by this factor.
  * `b_pos::Array{Float64, 1}`: Normalized span positions y/bscale of the
                                following distributions.
  * `chords::Array{Float64,1}`: Chord c/bscale distribution.
  * `twists::Array{Float64,1}`: Twist (deg) distribution.
  * `LE_x::Array{Float64,1}`  : x-position (chordwise) x/bscale of leading edge
                                distribution.
  * `LE_z::Array{Float64,1}`  : z-position (dihedral-wise) z/bscale of leading
                                edge distribution.
  * `crosssections::Array{Tuple{T,Array{T,2}}, 1}`    : cross sections along the
                                span in the form [(y/bscale, crosssection)],
                                where `crosssection` is a matrix that contains
                                all points of the airfoil contour indexed by
                                row. In order to define the resulting normals
                                pointing out of the geometry, the points in
                                `crosssection` must start from the trailing
                                edge, go around the top side towards the leading
                                edge, and back to the trailing edge around the
                                bottom side. ALL CROSS SECTIONS MUST HAVE THE
                                SAME NUMBER OF POINTS.


  **Optional Arguments**
  * `tilt_z::Array{Float64,1}`: Tilting (deg) about the z-axis of every span
                                cross section. This is also a distribution.
  * `symmetric::Bool`         : Whether to consider the `crosssections` to be
                                symmetric about the semi-span. If true, only
                                positive y/bscale are neeeded.

  NOTE: The resulting geometry will be a structured surface mesh with the first
  index going in chordwise direction starting from the TE around the bottom
  surface and around back to the TE. The second index is the in the spanwise
  direction starting at left end (lowest position along y axis) to right end.
"""
function generate_loft(crosssections::Array{Tuple{T,Array{T,2}}, 1},
                        bscale::Real, b_pos::Array{T,1}, chords::Array{T,1},
                        twists::Array{T,1},
                        LE_x::Array{T,1}, LE_z::Array{T,1};
                        # MORE GEOMETRIC OPTIONS
                        tilt_z=nothing, symmetric::Bool=false,
                        loop_dim::Int64=1,
                        # OUTPUT OPTIONS
                        save_path=nothing, paraview::Bool=true,
                        file_name::String="myloft"
                       ) where{T<:Real}

  sec_NDIVS = size(crosssections[1][2], 1)-1    # Cross-section cells
  b_NDIVS = size(b_pos, 1)-1                    # Span cells

  # ERROR CASES
  for (pos, sec) in crosssections
    if sec_NDIVS!=size(sec, 1)-1
      error("All airfoil sections must have the same number of points.")
    end
  end

  for (nam,val) in [("chords",chords), ("twists",twists), ("LE_x",LE_x),
                    ("LE_z",LE_z), ("tilt_z",tilt_z)]
    if nam=="tilt_z" && tilt_z==nothing
      nothing
    else
      if size(b_pos,1)!=size(val,1)
          error("The size of `b_pos` doesn't match the size of `$nam`."*
            " Expected size $(size(b_pos,1)), got $(size(val,1)).")
      end
    end
  end

  # ----------------- PARAMETRIC GRID ------------------------------------------
  P_min = [0, 0, 0]               # Lower boundary arclength, span, dummy
  P_max = [1, 1, 0]               # Upper boundary arclength, span, dummy
  # loop_dim = 1                  # Loops the arclength dimension

  # Adds dummy division
  NDIVS = [sec_NDIVS, b_NDIVS, 0]

  # Generates parametric grid
  grid = Grid(P_min, P_max, NDIVS, loop_dim)


  # ----------------- SURFACE GRID ---------------------------------------------
  # Auxiliary function for weighting values across span
  function calc_vals(span, array)

    # Finds bounding airfoil position
    val_in, val_out = nothing, array[1]
    for val in array[2:end]
        val_in = val_out
        val_out = val
        if val[1]>=span; break; end
    end
    pos_in = val_in[1]
    val_in = val_in[2]
    pos_out = val_out[1]
    val_out = val_out[2]

    weight = (span-pos_in)/(pos_out-pos_in)

    return weight, val_in, val_out
  end

  # Space transformation function
  function my_space_transform(inds)
    span = b_pos[inds[2]]           # y/bscale span position
    chord = chords[inds[2]]         # c/bscale chord length
    twist = twists[inds[2]]         # deg twist
    le_x = LE_x[inds[2]]            # x/bscale LE position
    le_z = LE_z[inds[2]]            # z/bscale LE position

    # Merges airfoil contours at this span position
    weight, sec_in, sec_out = calc_vals(span*sign(span)^symmetric, crosssections)

    # Point over airfoil contour
    point = weight*sec_out[inds[1], :]+(1-weight)*sec_in[inds[1], :]
    point = vcat(point, 0)

    # Scales the airfoil contour by the normalized chord length
    point = chord*point

    # Applies twist to the airfoil point
    tlt_z = tilt_z!=nothing ?  tilt_z[inds[2]] : 0.0
    point = rotation_matrix(-twist, -tlt_z, 0)*point

    # Places the point relative to LE and scales by span scale
    point = [point[1]+le_x, span+point[3], point[2]+le_z]*bscale

    return point
  end

  # Transforms the quasi-two dimensional parametric grid into the wing surface
  transform2!(grid, my_space_transform)

  if save_path!=nothing
    # Outputs a vtk file
    save(grid, file_name; path=save_path)

    if paraview
      # Calls paraview
      run(`paraview --data=$(joinpath(save_path,file_name)).vtk`)
    end
  end

  return grid::Grid
end






"""
`generate_loft(crosssections, bscale, b_low, b_up, b_NDIVS, chords, twists,
LE_x, LE_z; optargs...)`

  Generates a lofted surface geometry. The name of the arguments are taken from
  its initial application to the lofting of an aircraft wing, so bear with it.

  **Arguments**
  * `bscale::Float64`         : Semi-span scale. It will scale the entire
                                geometry by this factor. All y/bscale values in
                                the following arguments must go from 0 to 1.
  * `b_low::Float64`          : Scaled lower bound of the span.
  * `b_up::Float64`           : Scaled upper bound of the span. To generate
                                a symmetric wing, give it b_low=-1, b_up=1 and .
                                symmetric=true; for a semi-span, give it
                                b_low=0, b_up=1. If generating a prop blade,
                                give it b_low=Rhub/Rtip, b_up=1.
  * `b_NDIVS`                 : Number of divisions (cells) along span. This
                                if either an Int or an array of sections in the
                                format of `multidiscretize()`.
  * `chords::Array{Float64,2}`: Chord distribution along the span in the form
                                [y/bscale c/bscale].
  * `twists::Array{Float64,2}`: Twist distribution along the span in the form
                                [y/bscale deg].
  * `LE_x::Array{Float64,2}`  : x-position (chordwise) of leading edge along the
                                span in the form [y/bscale x/bscale].
  * `LE_z::Array{Float64,2}`  : z-position (dihedral-wise) of leading edge along
                                the span in the form [y/bscale z/bscale].
  * `crosssections::Array{Tuple{T,Array{T,2}}, 1}`    : cross sections along the
                                span in the form [(y/bscale, crosssection)],
                                where `crosssection` is a matrix that contains
                                all points of the airfoil contour indexed by
                                row. In order to define the resulting normals
                                pointing out of the geometry, the points in
                                `crosssection` must start from the trailing
                                edge, go around the top side towards the leading
                                edge, and back to the trailing edge around the
                                bottom side. ALL CROSS SECTIONS MUST HAVE THE
                                SAME NUMBER OF POINTS.


  **Optional Arguments**
  * `tilt_z::Array{Float64,2}`: Tilting about the z-axis of every span cross
                                section in the form [(y/bscale, deg)].
  * `spl_k`, `spl_bc`, `spl_s`: Spline parameters with k the degree of the
                                spline, bc the boundary condition, and s the
                                smoothing or error of the spline.

  NOTE: The resulting geometry will be a structured surface mesh with the first
  index going in chordwise direction starting from the TE around the bottom
  surface and around back to the TE. The second index is the in the spanwise
  direction starting at left end (lowest position along y axis) to right end.
"""
function generate_loft(crosssections::Array{Tuple{T,Array{T,2}}, 1},
                        bscale::Real, b_low::Real, b_up::Real, b_NDIVS,
                        chords::Array{T,2}, twists::Array{T,2},
                        LE_x::Array{T,2}, LE_z::Array{T,2};
                        # MORE GEOMETRIC OPTIONS
                        tilt_z=nothing, symmetric::Bool=false,
                        # SPLINE OPTIONS
                        spl_k::Int64=5, spl_bc::String="extrapolate",
                        spl_s::Real=0.001, verify_spline::Bool=true,
                        optargs...
                       ) where{T<:Real}

 # Normalized span positions
 if typeof(b_NDIVS)==Int64
   b_pos = discretize(x->x, b_low, b_up, b_NDIVS, 1.0)
 elseif typeof(b_NDIVS)==multidisctype
   b_pos = multidiscretize(x->x, b_low, b_up, b_NDIVS)
 else
   error("Invalid b_NDIVS=$b_NDIVS."*
          " Expected type $(Int64) or $(multidiscrtype), got $(typeof(b_NDIVS))")
 end
 b_pos = T[val for val in b_pos]


  # ----------------- GEOMETRY SPLINES -----------------------------------------
  # Splines all distributions for a smooth geometry
  _spl_chord = Dierckx.Spline1D(chords[:, 1], chords[:, 2]; s=spl_s, bc=spl_bc,
                          k=size(chords,1)>=spl_k ? spl_k : size(chords,1)-1 )
  _spl_twist = Dierckx.Spline1D(twists[:, 1], twists[:, 2];s=spl_s, bc=spl_bc,
                          k=size(twists,1)>=spl_k ? spl_k : size(twists,1)-1)
  _spl_LE_x = Dierckx.Spline1D(LE_x[:, 1], LE_x[:, 2]; s=spl_s, bc=spl_bc,
                          k=size(LE_x,1)>=spl_k ? spl_k : size(LE_x,1)-1)
  _spl_LE_z = Dierckx.Spline1D(LE_z[:, 1], LE_z[:, 2]; s=spl_s, bc=spl_bc,
                          k=size(LE_z,1)>=spl_k ? spl_k : size(LE_z,1)-1)
  if tilt_z!=nothing
   _spl_tlt_z = Dierckx.Spline1D(tilt_z[:, 1], tilt_z[:, 2]; s=spl_s, bc=spl_bc,
                          k=size(tilt_z,1)>=spl_k ? spl_k : size(tilt_z,1)-1)
    new_tlt_z = [_spl_tlt_z(symmetric ? abs(b) : b) for b in b_pos]
  else
    new_tlt_z = nothing
  end

  new_chords = [_spl_chord(symmetric ? abs(b) : b) for b in b_pos]
  new_twists = [_spl_twist(symmetric ? abs(b) : b) for b in b_pos]
  new_LE_x = [_spl_LE_x(symmetric ? abs(b) : b) for b in b_pos]
  new_LE_z = [_spl_LE_z(symmetric ? abs(b) : b) for b in b_pos]


  # # ----------------- SPLINE VERIFICATION --------------------------------------
  # if verify_spline
  #  fig = plt.figure("spl_verif", figsize=(7*2,5*1))
  #
  #  plt.subplot(121)
  #  plt.plot(LE_x[:,1], LE_x[:,2], "og", label="Org LE x", alpha=0.5)
  #  plt.plot(LE_z[:,1], LE_z[:,2], "ob", label="Org LE z", alpha=0.5)
  #  plt.plot(b_pos, new_LE_x, "--g", label="Spline LE x")
  #  plt.plot(b_pos, new_LE_z, "--b", label="Spline LE z")
  #  plt.xlabel(plt.L"y/b_{scale}")
  #  plt.ylabel(plt.L"x/b_{scale}, z/b_{scale}")
  #  plt.grid(true, color="0.8", linestyle="--")
  #  plt.legend(loc="best")
  #
  #  plt.subplot(122)
  #  p1 = plt.plot(twists[:,1], twists[:,2], "og", label="Org Twist", alpha=0.5)
  #  p2 = plt.plot(b_pos, new_twists, "--g", label="Spline twist")
  #  pextra = []
  #  if tilt_z!=nothing
  #    pextra1 = plt.plot(tilt_z[:,1], tilt_z[:,2], "or", label="Org tilt z",
  #                                                                    alpha=0.5)
  #    pextra2 = plt.plot(b_pos, new_tlt_z, "--r", label="Spline tilt z")
  #    pextra = vcat(pextra, [pextra1[1], pextra2[1]])
  #  end
  #  plt.ylabel("Twist (deg)")
  #
  #  plt.grid(true, color="0.8", linestyle="--")
  #  plt.xlabel(plt.L"y/b_{scale}")
  #
  #  plt.twinx()
  #  p3 = plt.plot(chords[:,1], chords[:,2], "ob", label="Org Chord", alpha=0.5)
  #  p4 = plt.plot(b_pos, new_chords, "--b", label="Spline chord")
  #  plt.ylabel(plt.L"c/b_{scale}")
  #
  #  ps = vcat([p1[1], p2[1], p3[1], p4[1]], pextra)
  #  plt.legend(ps, [p[:get_label]() for p in ps], loc="best")
  # end

  # ----------------- LOFTING --------------------------------------------------
  return generate_loft(crosssections, bscale, b_pos, new_chords, new_twists,
                        new_LE_x, new_LE_z; tilt_z=new_tlt_z,
                        symmetric=symmetric, optargs...)

end





"""
`generate_loft(crosssections, upper_rfl_NDIVS, lower_rfl_NDIVS, args...;
rflspl_k::Int64=5, rflspl_s::Real=0.001, verify_rflspline::Bool=true,
rfloptargs=[], optargs...)`

  This function also rediscretizes the cross sections as indicated by
  `upper_rfl_NDIVS`, and `lower_rfl_NDIVS` (upper and lower surface sections,
  respectively). Hence, the original cross sections need not to have the same
  number of points.
"""
function generate_loft(crosssections::Array{Tuple{T,Array{T,2}}, 1},
                            upper_rfl_NDIVS, lower_rfl_NDIVS,
                            bscale::Real, b_low::Real, b_up::Real, b_NDIVS,
                            chords::Array{T,2}, twists::Array{T,2},
                            LE_x::Array{T,2}, LE_z::Array{T,2};
                            # AIRFOIL SPLINE OPTIONS
                            rflspl_k::Int64=5, rflspl_s::Real=0.001,
                            verify_rflspline::Bool=true, rfloptargs=[],
                            optargs...
                           ) where{T<:Real}

  # Formats the NDIVS
  if typeof(upper_rfl_NDIVS)==Int64
    urfl_NDIVS = [(1.0, upper_rfl_NDIVS, 1.0, false)]
  elseif typeof(upper_rfl_NDIVS)==multidisctype
    urfl_NDIVS = upper_rfl_NDIVS
  else
    error("Invalid upper_rfl_NDIVS=$upper_rfl_NDIVS. Expected type $(Int64)"*
            " or $(multidiscrtype), got $(typeof(upper_rfl_NDIVS))")
  end
  if typeof(lower_rfl_NDIVS)==Int64
    lrfl_NDIVS = [(1.0, lower_rfl_NDIVS, 1.0, false)]
  elseif typeof(lower_rfl_NDIVS)==multidisctype
    lrfl_NDIVS = lower_rfl_NDIVS
  else
    error("Invalid lower_rfl_NDIVS=$lower_rfl_NDIVS. Expected type $(Int64)"*
            " or $(multidiscrtype), got $(typeof(lower_rfl_NDIVS))")
  end

  # Rediscretizes cross sections
  rediscr = []

  for (pos, M) in crosssections
    new_xy = rediscretize_airfoil(M[:, 1], M[:, 2], urfl_NDIVS, lrfl_NDIVS;
                                    spl_k=rflspl_k, spl_s=rflspl_s,
                                    verify_spline=verify_rflspline,
                                    rfloptargs...)
    push!(rediscr, [new_xy[j][i] for i in 1:size(new_xy[1],1), j in 1:2])
  end

  # Formates them back to the original format
  new_crosssections = [ (pos, rediscr[i])
                                  for (i,(pos, M)) in enumerate(crosssections)]

  return generate_loft(new_crosssections, bscale, b_low, b_up, b_NDIVS,
                          chords, twists, LE_x, LE_z; optargs...)
end


"""
`generate_loft(crosssections::Array{Tuple{T,String}, 1},
                        data_path::String, args...; header_len::Int64=1,
                        delim::String=" ", optargs...) where{T<:Real}`

Loft a geometry where the cross sections are read from the files indicated by
`crosssections` found in `data_path`.
"""
function generate_loft(crosssections::Array{Tuple{T,String}, 1},
                        data_path::String, args...; header_len::Int64=1,
                        delim::String=" ", optargs...) where{T<:Real}

 secs = [(pos, readcontour(f_name; path=data_path, header_len=header_len,
                              delim=delim, output="matrix"))
                              for (pos, f_name) in crosssections]

 return generate_loft(secs, args...; optargs...)
end





# REVOLTUION ###################################################################
"""
  `surface_revolution(profile, thetaNDIVS; loop_dim=0, axis_angle=0, low_a=0,
up_a=360, save_path=nothing, paraview=true, file_name="myrev")`

Receives a contour to revolve around an axis generating a surface of revolution.

**ARGUMENTS**
  * `profile::Array{Float64,2}`   : Contour to revolve. These are two-
                                    dimensional points in the YZ-plane that will
                                    get revolve around the Z-axis.
  * `thetaNDIVS::Int64`           : Number of angle-sections (cells) in the
                                    revolution.

**OPTIONAL ARGUMENTS**
  * `loop_dim::Int64=0`           : Whether to loop any dimension of the
                                    parametric grid.
  * `axis_angle::Float64=0`       : Tilting angle (deg) about the Z-axis to
                                    revolve the contour.
  * `low_a::Float64=0`            : Lower bound of angle (deg) of revolution.
  * `up_a::Float64=0`             : Upper bound of angle (deg) of revolution.
"""
function surface_revolution(profile::Array{T,2}, thetaNDIVS::Int64;
                              loop_dim::Int64=0, axis_angle::Real=0,
                              low_a::Real=0, up_a::Real=360,
                              # OUTPUT OPTIONS
                              save_path=nothing, paraview::Bool=true,
                              file_name::String="myrev"
                              ) where{T<:Real}
  #ERROR CASES
  if size(profile,2)!=2
    error("Invalid point dimensions in `profile`."*
          " Expected 2 dimensions, got $(size(profile,2))")
  elseif profile[1,:]==profile[end,:] && loop_dim==0
    warn("Received a closed contour but parametric grid wasn't declared to"*
          " loop, resulting in overlaping start/end points. Give it"*
          " `loop_dim=1` to fix that.")
  end

  # First coordinate is a map to the points in the countor of revolution
  # Seciond coordinate is the angle of revolution
  # Third is a dummy
  NDIVS = [size(profile,1)-1, thetaNDIVS, 0]# Divisions in every coordinate
  P_min = [0, low_a, 0]                     # Lower boundaries
  P_max = [1, up_a, 0 ]                     # Upper boundaries
  # loop_dim = 1                            # Loops the countour of revolution

  # Parametric grid
  grid = Grid(P_min, P_max, NDIVS, loop_dim)

  # Rotates the profile according the revolution tilting
  if axis_angle!=0
    M = rotation_matrix(0,0,-axis_angle)
    M2D = M[2:3,2:3]
    M3D = M'
    points = hcat([M2D*profile[i,:] for i in 1:size(profile,1)]...)'
  else
    M3D = eye(3)
    points = profile
  end

  # Space transformation function
  function my_space_transform(X, ind)
      p_ind = ind[1]                  # Point index
      angle = X[2]                    # (deg) angle

      # Places the contour in the YZ-plane
      point = [0, points[p_ind,1], points[p_ind,2]]

      # Rotates the contour around the Z-axis
      return M3D*axis_rotation(Float64[0,0,1], Float64(angle))*point
  end

  # Applies the space transformation to the parametric grid
  transform3!(grid, my_space_transform)


  if save_path!=nothing
    # Outputs a vtk file
    save(grid, file_name; path=save_path)

    if paraview
      # Calls paraview
      run(`paraview --data=$(joinpath(save_path,file_name)).vtk`)
    end
  end

  return grid::Grid
end



















################################################################################
