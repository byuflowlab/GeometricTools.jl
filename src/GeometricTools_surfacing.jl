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
  its initial application to the lofting of an aircraft wing.

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
                        tilt_z=nothing, tilt_y=nothing,
                        symmetric::Bool=false,
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
                    ("LE_z",LE_z), ("tilt_z",tilt_z), ("tilt_y",tilt_y)]
    if (nam=="tilt_z" && tilt_z==nothing) || (nam=="tilt_y" && tilt_y==nothing)
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

    inds1 = (sec_NDIVS+1) - (inds[1]-1) # Here it flips the airfoil contour

    # Merges airfoil contours at this span position
    weight, sec_in, sec_out = calc_vals(span*sign(span)^symmetric, crosssections)

    # Point over airfoil contour
    point = weight*sec_out[inds1, :]+(1-weight)*sec_in[inds1, :]
    point = vcat(point, 0)

    # Scales the airfoil contour by the normalized chord length
    point = chord*point

    # Applies twist to the airfoil point
    tlt_z = tilt_z!=nothing ?  tilt_z[inds[2]] : 0.0
    tlt_y = tilt_y!=nothing ?  tilt_y[inds[2]] : 0.0
    point = rotation_matrix(-twist, -tlt_z, -tlt_y)*point

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
  its initial application to the lofting of an aircraft wing.

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


  # ----------------- SPLINE VERIFICATION --------------------------------------
  if verify_spline
    # if isdefined(Main, :PyPlot)
      plot_loft_splines(b_pos, LE_x, LE_z, new_LE_x, new_LE_z, twists, new_twists,
        tilt_z, new_tlt_z, chords, new_chords)
    # else
    #   @eval Main import PyPlot
    #   Base.invokelatest(plot_loft_splines, b_pos, LE_x, LE_z, new_LE_x, new_LE_z, twists, new_twists,
    #     tilt_z, new_tlt_z, chords, new_chords)
    # end

  end

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









# PATH LOFTING #################################################################
"""
`surface_pathloft(sections, section_NDIVS, path_NDIVS, path; injective_coordinate=1,
loop_dim=0, nperiodic=2, redisc_optargs=[], verify_spline=true,
save_path=nothing, paraview=true, file_pref="pathloft")`

Generation of a surface grid by lofting a set of sections along a path.

`sections` is an array where `sections[i] = (spos, points)` is the i-th
contour made out of the 2xn matrix `points` that goes along the path in the
nondimensional arc-length position `spos` (a number between 0 and 1, where 0
is the start of the path and 1 is the end of the path). `path` is an array where
`path[i] = (X, n, twist)` is the i-th point along the
path placed at `X` in the global coordinate system and the cross section passing
through this point will be oriented according to the normal vector `n` and
twisted by `twist` (in degrees). The number of points in `path` does not need to
match the number of sections in `sections`

This function rediscretizes the sections as indicated by `section_NDIVS`,
which could be an array of multidiscretization parameters to apply to each
section, a single multidiscretization parameter set for all the sections, or an
integer with the number of points used to discretize all sections. Similarly,
`path_NDIVS` indicates the desired discretization of the x-position along the path.

Contours can either be closed or open, and the discretization will automatically
recognize and preserve that. Contours must be concave.

`sections` has the format `sections[i] = (spos, points)` where `points`
is the nx2 matrix of n contours points (y, z) of each section, and `spos`
is the nondimensional position along the path of each section.

`path` has the format `path[i] = (X, normal, twist)`, where `X` is the position
of the i-th point (a vector of length 3), `normal` is a the normal of cross
sections at that position (a vector of length 3), and `twist` is the twist of
the cross section at that position (in degrees). `injective_coordinate` is the
coordinate that makes the path injective.
"""
function surface_pathloft(sections::AbstractVector,
                            path::AbstractVector{T},
                            sec_discretizations::AbstractVector{multidisctype},
                            xpositions::AbstractVector{<:Number};
                            # PATH PARAMETERS
                            injective_coordinate=1,
                            loop_dim=0,
                            # SECTION PARAMETERS
                            nperiodic=2,
                            sort_rediscretization=false,   # Sort rediscretized sections to align starting points
                            redisc_optargs=[],
                            # OUTPUT PARAMETERS
                            verify_spline=true,
                            save_path=nothing, paraview=true,
                            output_path_normals=false,
                            output_path_Xs=false,
                            output_path_xpos=false,
                            file_pref="pathloft"
                            ) where {T<:Tuple{<:AbstractVector, <:AbstractVector, <:Number}}

    str = ""
    nsecs = length(sections)
    npoints = sum(disc[2] for disc in sec_discretizations[1])

    # Error cases
    @assert length(sec_discretizations)==nsecs ""*
        "Length of sec_discretizations ($(length(sec_discretizations)))"*
        " different than the number of sections $(nsecs)"

    @assert prod(sum(disc[2] for disc in sec_discretization)==npoints
                            for sec_discretization in sec_discretizations) ""*
        "Requested discretization of sections with varying number of points,"*
        " but implementation requires the same number of points for all"*
        " sections.\nsec_discretizations = $(sec_discretizations)"

    @assert prod(length(X)==3 for (X, n, twist) in path) ""*
        "Received an `X` position vector along path that is not length 3"

    @assert prod(length(n)==3 for (X, n, twist) in path) ""*
        "Received an `n` normal vector along path that is not length 3"

    @assert prod(isapprox(norm(n), 1) for (X, n, twist) in path) ""*
        "Received normal in path that is not unitary;"*
        " n_norm = $([ norm(n) for (X, n, twist) in path ])"


    # ----------------- DISCRETIZE SECTIONS ------------------------------------

    # Rediscretize each section unto a uniform grid of arc lengths wider than
    # the physical arc length since later they will be used to spline and
    # discretize again in 2D
    npoints_uniform = 3*maximum(size(points, size(points, 1)==2 ? 2 : 1) for (xpos, points) in sections)

    s_lo = -(nperiodic - 1)
    s_up = 1 + (nperiodic - 1)

    discretization = [(1.0, npoints_uniform, 1.0, true)]
    pre_discretization = [(1.0, 10*npoints_uniform, 1.0, true)]

    outs = []
    for (xpos, points) in sections

        new_points = rediscretize_concavecontour(points,
                                                    sort_rediscretization ? pre_discretization : discretization;
                                                    s_lo=sort_rediscretization ? 0 : s_lo,
                                                    s_up=sort_rediscretization ? 1 : s_up,
                                                    nperiodic=nperiodic,
                                                    verify_spline=verify_spline,
                                                    plot_title=plt.L"$x=$"*"$(xpos)",
                                                    out=sort_rediscretization ? nothing : outs,
                                                    redisc_optargs...)

        # Sort rediscretized section to (hopefully) align starting points
        if sort_rediscretization

            # Force the rediscretized contour to be close if the original one was closed
            if size(points, 1)==2
                if points[1, 1]==points[1, end] && points[2, 1]==points[2, end]
                    new_points[1, 1] = new_points[1, end]
                    new_points[2, 1] = new_points[2, end]
                end
            else
                if points[1, 1]==points[end, 1] && points[1, 2]==points[end, 2]
                    new_points[1, 1] = new_points[end, 1]
                    new_points[1, 2] = new_points[end, 2]
                end
            end

            new_points = rediscretize_concavecontour(new_points, discretization;
                                                        s_lo=s_lo, s_up=s_up,
                                                        nperiodic=nperiodic,
                                                        verify_spline=false,
                                                        out=outs,
                                                        redisc_optargs...)


        end

    end

    # Fetch the x-position of each original section
    sec_xpos = [xpos for (xpos, points) in sections]

    # Fetch the uniform arc length used to re-discretize the original sections
    # sec_arclength = outs[1].new_arclengths
    # sec_arclength = Float64.(outs[1].new_arclengths)
    sec_arclength = [val for val in outs[1].new_arclengths]

    # Evenly-spaced range between 0 and 1 that maps the arc-length
    # discretization of each section
    sec_range = collect(range(0, 1, length=npoints+1))

    # Spline and discretize the arclengths of each original section,
    # interpolating the arclengths in between sections
    all_new_arclengths = hcat([multidiscretize(identity, 0, 1, disc) for disc in sec_discretizations]...)
    grid_arclengths = FLOWMath.interp2d(FLOWMath.akima,
                                        sec_range, sec_xpos, all_new_arclengths,
                                        sec_range, xpositions)

    # Spline and discretize rho and theta of the uniform sections,
    # interpolating the arclengths in between sections
    all_new_rhos = hcat([out.new_rhos for out in outs]...)
    all_new_thetas = hcat([out.new_thetas for out in outs]...)
    grid_rhos = []
    grid_thetas = []
    for (xpos, arclengths) in zip(xpositions, eachcol(grid_arclengths))

        this_rhos = FLOWMath.interp2d(FLOWMath.akima,
                                            sec_arclength, sec_xpos, all_new_rhos,
                                            arclengths, xpos)

        this_thetas = FLOWMath.interp2d(FLOWMath.akima,
                                            sec_arclength, sec_xpos, all_new_thetas,
                                            arclengths, xpos)

        push!(grid_rhos, this_rhos)
        push!(grid_thetas, this_thetas)
    end
    grid_rhos = hcat(grid_rhos...)
    grid_thetas = hcat(grid_thetas...)

    # Convert cylindrical coordinates back to Cartesian coordinates
    grid_ys = grid_rhos .* sin.(grid_thetas)
    grid_zs = grid_rhos .* cos.(grid_thetas)

    # ----------------- REDISCRETIZE PATH --------------------------------------
    # Spline and rediscretize the points along path to match the grid x-positions
    path_Xs = [X for (X, n, twist) in path]

    rediscretize_out = []
    new_path_Xs = rediscretize_line(path_Xs, xpositions;
                                      parameterization = X->X[injective_coordinate],
                                      out=rediscretize_out)

    # Fetch the arc lengths corresponding the each of the original points from
    # the intermediate calculations
    path_arclengths = rediscretize_out[1].org_arclengths
    new_path_arclengths = rediscretize_out[1].new_arclengths
    total_arclength = rediscretize_out[1].total_arclength

    # Spline normal and twist and probe it at the new x-positions
    path_n1 = [n[1] for (X, n, twist) in path]
    path_n2 = [n[2] for (X, n, twist) in path]
    path_n3 = [n[3] for (X, n, twist) in path]
    path_twist = [twist for (X, n, twist) in path]

    new_path_n1 = FLOWMath.akima(path_arclengths, path_n1, xpositions)
    new_path_n2 = FLOWMath.akima(path_arclengths, path_n2, xpositions)
    new_path_n3 = FLOWMath.akima(path_arclengths, path_n3, xpositions)
    new_path_nnorm = sqrt.(new_path_n1.^2 + new_path_n2.^2 + new_path_n3.^2)
    new_path_twist = FLOWMath.akima(path_arclengths, path_twist, xpositions)

    # Re-pack the new path
    new_path = [(X, [n1/nnorm, n2/nnorm, n3/nnorm], twist)
                    for (X, n1, n2, n3, nnorm, twist) in
                    zip(new_path_Xs, new_path_n1, new_path_n2, new_path_n3, new_path_nnorm, new_path_twist)]

    # Verification plots
    if verify_spline
        fig = plt.figure("Path discretization", figsize=[7*2, 5]*7/9)
        axs = fig.subplots(1, 2)

        fig2 = plt.figure("Path discretization - scalar", figsize=[7, 5]*1/2)
        ax2 = fig2.gca()

        for (ax, dimi) in zip(axs, 2:3)

            # Define scaling for normal vector
            scale = 0.1*total_arclength

            # Plot raw path
            ax.plot([X[1] for X in path_Xs], [X[dimi] for X in path_Xs],
                                    ".k", label="Raw", alpha=0.8, markersize=8,
                                    clip_on=false)

            for (X, n, twist) in path
                ax.plot([X[1], X[1]+scale*n[1]],
                        [X[dimi], X[dimi]+scale*n[dimi]],
                        "-k", alpha=0.5, linewidth=1, clip_on=false)
            end

            # Plot splined path
            ax.plot([X[1] for (X, n, twist) in new_path],
                    [X[dimi] for (X, n, twist) in new_path],
                                    "-r", label="Spline",
                                    alpha=1.0, markersize=4, linewidth=1,
                                    clip_on=false)
            for (X, n, twist) in new_path
                ax.plot([X[1], X[1]+scale*n[1]],
                        [X[dimi], X[dimi]+scale*n[dimi]],
                        "-g", alpha=0.25, linewidth=1, clip_on=false)
            end

            ax.set_xlabel(L"x")
            ax.set_ylabel([L"x", L"y", L"z"][dimi])

        end

        # Plot twist
        ax2.plot(path_arclengths*total_arclength,
                [twist for (X, n, twist) in path],
                                ".k", label="Raw", alpha=0.8, markersize=8,
                                clip_on=false)

        ax2.plot(new_path_arclengths*total_arclength,
                [twist for (X, n, twist) in new_path],
                                        "-r", label="Spline",
                                        alpha=1.0, markersize=4, linewidth=1,
                                        clip_on=false)

        ax2.set_xlabel(L"Arc length $s$")
        ax2.set_ylabel(L"Twist ($^\circ$)")

        for (axi, ax) in enumerate(vcat(axs, ax2))
            if axi != 3; ax.set_aspect("equal"); end
            ax.spines["right"].set_visible(false)
            ax.spines["top"].set_visible(false)
            ax.legend(loc="best", frameon=true, fontsize=8)
        end
    end

    # Output VTK of path
    if !isnothing(save_path)

        point_data = [Dict(
                            "field_name" => "normal",
                            "field_type" => "vector",
                            "field_data" => [n for (X, n, twist) in path]
                            ),
                      Dict(
                            "field_name" => "twist",
                            "field_type" => "scalar",
                            "field_data" => [twist for (X, n, twist) in path]
                            )
                    ]
        str *= generateVTK(file_pref*"_path_raw", path_Xs;
                                lines=[collect(0:length(path_Xs)-1)],
                                point_data=point_data,
                                path=save_path)

        point_data = [Dict(
                            "field_name" => "normal",
                            "field_type" => "vector",
                            "field_data" => [n for (X, n, twist) in new_path]
                            ),
                      Dict(
                            "field_name" => "twist",
                            "field_type" => "scalar",
                            "field_data" => [twist for (X, n, twist) in new_path]
                            )
                    ]
        str *= generateVTK(file_pref*"_path_spline", new_path_Xs;
                                lines=[collect(0:length(new_path_Xs)-1)],
                                point_data=point_data,
                                path=save_path)
    end


    # ----------------- PLACE SECTIONS ALONG PATH ------------------------------
    # TODO

    # ----------------- PARAMETRIC GRID ----------------------------------------
    # Create quasi-3D grid in non-dimensional space (path x, section arclength)
    Pmin = [0, 0, 0.0]                        # Lower boundary x, sec arc, dummy
    Pmax = [1, 1, 0.0]                        # Upper boundary x, sec arc, dummy
    NDIVS = [length(xpositions)-1, npoints, 0]  # Divisions

    # Generate parametric grid
    grid = Grid(Pmin, Pmax, NDIVS, loop_dim)

    nodes = grid.nodes

    # Rotation matrix of each normal
    angles = [ ( 0, 180/pi*acos(n[3]) - 90, 180/pi*sign(n[2] + 3*eps())*acos( n[1]/sqrt(n[1]^2+n[2]^2) ) ) for (X, n, twist) in new_path]
    rotationmatrix = [rotation_matrix2( angle... ) for angle in angles]

    # println([round.(n; digits=2) for (X, n, twist) in new_path])
    # println([round.(angle; digits=2) for angle in angles])

    # Convert nodes into pre-calculated grid
    function my_space_transform(inds)

        i = inds[1]                             # Index along path coordinate
        j = inds[2]                             # Index along contour coordinate

        spos = xpositions[i]                    # Position along arc length
        y = grid_ys[j, i]                       # y-coordinate of contour point
        z = grid_zs[j, i]                       # z-coordinate of contour point

        O = new_path[i][1]                      # Contour center in global coordinates
        n = new_path[i][2]                      # Normal of the control in global coordinates
        twist = new_path[i][3]                  # (deg) twist of the contour

        M = rotationmatrix[i]                   # Rotation matrix according to normal vector

        # Twist contour
        y, z = y*cosd(twist) + z*sind(twist), -y*sind(twist) + z*cosd(twist)

        # Rotate contour according to normal
        X = M*[0, y, z]

        # Translate contour to the position along the path
        point = O + X

        return point
    end

    # Transforms the quasi-two dimensional parametric grid into the wing surface
    transform2!(grid, my_space_transform)

    # Store the path normal used for each point
    if output_path_normals
        pathnormals = [new_path[CartesianIndices(grid._ndivsnodes)[i][1]][2] for i in 1:grid.nnodes]
        add_field(grid, "pathnormal", "vector", pathnormals, "node")
    end

    # Store pointer to path from each point
    if output_path_Xs
        pathXpointers = [new_path[CartesianIndices(grid._ndivsnodes)[i][1]][1] - get_node(grid, i) for i in 1:grid.nnodes]
        add_field(grid, "pathX", "vector", pathXpointers, "node")
    end

    # Store non-dimensional position along path for each point
    if output_path_xpos
        pathxpos = [xpositions[CartesianIndices(grid._ndivsnodes)[i][1]] for i in 1:grid.nnodes]
        add_field(grid, "pathposition", "scalar", pathxpos, "node")
    end

    # Output VTK of the loft
    if !isnothing(save_path)

        str *= save(grid, file_pref*"_loft"; path=save_path, format="vtk")

        if paraview
            run(`paraview --data=$(joinpath(save_path, str))`)
        end
    end

    return grid

end


function surface_pathloft(sections::AbstractVector, path,
                            sec_discretizations::AbstractVector{multidisctype},
                            spos_NDIVS::multidisctype,
                            args...;
                            spos_lo=0.0, spos_up=1.0,
                            optargs...
                            )

    xpositions = Float64.(multidiscretize(identity, spos_lo, spos_up, spos_NDIVS))

    return surface_pathloft(sections, path, sec_discretizations, xpositions,
                                                            args...; optargs...)

end
function surface_pathloft(sections::AbstractVector, path,
                            sec_discretizations::AbstractVector{multidisctype},
                            spos_NDIVS::Int,
                            args...; optargs...
                            )

    _spos_NDIVS = [(1.0, spos_NDIVS, 1.0, true)]

    return surface_pathloft(sections, path, sec_discretizations, _spos_NDIVS,
                                                            args...; optargs...)

end
function surface_pathloft(sections::AbstractVector, path,
                            sec_discretization::multidisctype,
                            spos_NDIVS::Union{Int, multidisctype},
                            args...; optargs...
                            )

    sec_discretizations = [sec_discretization for i in 1:length(sections)]

    return surface_pathloft(sections, path, sec_discretizations, spos_NDIVS,
                                                            args...; optargs...)
end
function surface_pathloft(sections::AbstractVector, path,
                            sec_NDIVS::Int,
                            spos_NDIVS::Union{Int, multidisctype},
                            args...; optargs...
                            )

    _sec_NDIVS = [(1.0, sec_NDIVS, 1.0, true)]

    return surface_pathloft(sections, path, _sec_NDIVS, spos_NDIVS,
                                                            args...; optargs...)
end

"""
`surface_pathloft(sections::Vector{Tuple{<:Real, String}}, data_path::String,
args...; header_len=1, delim=",", optargs...)`

Loft a geometry where the sections are read from the files indicated by
`sections` found in `data_path`.

`sections` has the format `sections[i] = (xpos, filename)` where `filename`
is the CSV file containing the (y, z) contours of each section, and `xpos`
is the nondimensional position along the path of each section.
"""
function surface_pathloft(sections::AbstractVector{Tuple{<:Number, String}},
                            data_path::String, args...; header_len::Int64=1,
                            delim::String=",", optargs...)

    secs = [
            (xpos, readcontour(f_name; path=data_path, header_len=header_len,
                                                    delim=delim, output="matrix"))
            for (xpos, f_name) in sections
            ]

    return surface_pathloft(secs, args...; optargs...)
end








# REVOLUTION ###################################################################
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
function surface_revolution(profile::AbstractMatrix{T}, thetaNDIVS::TN;
                              loop_dim::Integer=0, axis_angle::Number=0,
                              low_a::Number=0, up_a::Number=360,
                              # OUTPUT OPTIONS
                              save_path=nothing, paraview::Bool=true,
                              file_name::AbstractString="myrev"
                              ) where{T<:Real, TN}
  #ERROR CASES
  if size(profile,2)!=2
    error("Invalid point dimensions in `profile`."*
          " Expected 2 dimensions, got $(size(profile,2))")
  elseif profile[1,:]==profile[end,:] && loop_dim==0
    @warn("Received a closed contour but parametric grid wasn't declared to"*
          " loop, resulting in overlaping start/end points. Give it"*
          " `loop_dim=1` to fix that.")
  end

  # First coordinate is a map to the points in the countor of revolution
  # Seciond coordinate is the angle of revolution
  # Third is a dummy
  if TN <: Number
      NDIVS = [size(profile,1)-1, thetaNDIVS, 0] # Divisions in every coordinate
  else
      NDIVS = [[(1.0, size(profile,1)-1, 1.0, false)], thetaNDIVS, [(1.0, 0, 1.0, false)]]
  end

  P_min = [0, low_a, 0]                     # Lower boundaries
  P_max = [1, up_a, 0 ]                     # Upper boundaries
  # loop_dim = 1                            # Loops the countour of revolution

  # Parametric grid
  grid = Grid(P_min, P_max, NDIVS, loop_dim)

  # Rotates the profile according the revolution tilting
  if axis_angle!=0
    M = rotation_matrix(0,0,-axis_angle)
    M2D = M[2:3,2:3]
    M3D = collect(M)'
    points = collect(hcat([M2D*profile[i,:] for i in 1:size(profile,1)]...))'
  else
    M3D = I
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
