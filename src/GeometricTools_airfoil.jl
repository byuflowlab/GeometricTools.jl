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
function readcontour(file_name; header_len=1, delim=" ", path="", output="arrays")

    x, y = Float64[], Float64[]

    open(joinpath(path,file_name)) do f
        for (i,line) in enumerate(eachline(f))

            # Ignores header
            if i<=header_len
                nothing
                # Parses each line
            else
                this_x, this_y = split(line, delim; keepempty=false)
                push!(x, Base.parse(Float64, this_x))
                push!(y, Base.parse(Float64, this_y))
            end

        end
    end

    if output=="arrays"
        return x,y
    elseif output=="matrix"
        return hcat(x,y)
    else
        error("Invalid `output` argument $(output).")
    end
end
#TODO: This is type unstable, since it returns either a vector or a matrix.  Should
# we make this type stable?

function rediscretize_airfoil(x::AbstractVector{T}, y::AbstractVector{T},
    upperNDIVS::multidisctype, lowerNDIVS::multidisctype; spl_s::Real=0.00001,
    spl_k::Integer=4, verify_spline::Bool=true, pltargs...) where{T<:Real}

    # Separate upper and lower sides to make the contour injective in x
    upper, lower = splitcontour(x, y)

    # Parameterize both sides independently
    fun_upper = parameterize(upper[1], upper[2], zeros(eltype(upper[1]), size(upper[1])); inj_var=1,
                                                      s=spl_s, kspl=spl_k)
    fun_lower = parameterize(lower[1], lower[2], zeros(eltype(lower[1]), size(lower[1])); inj_var=1,
                                                      s=spl_s, kspl=spl_k)

    # Discretizes both sides
    new_upper = multidiscretize(fun_upper, 0, 1, upperNDIVS)
    new_lower = multidiscretize(fun_lower, 0, 1, lowerNDIVS)

    # Merges sides back together
    points = vcat(reverse(new_upper), new_lower[2:end])
    new_x = [p[1] for p in points]
    new_y = [p[2] for p in points]

    # Plots
    if verify_spline
        # if isdefined(Main, :PyPlot)
            plot_airfoil(x, y; label="Original", style="--^k", alpha=0.5, pltargs...)
            plot_airfoil(new_x, new_y; label="Parameterized", style=":.b", pltargs...)
        # else
        #     @eval Main import PyPlot
        #     Base.invokelatest(plot_airfoil, x, y; label="Original", style="--^k", alpha=0.5, pltargs...)
        #     Base.invokelatest(plot_airfoil, new_x, new_y; label="Parameterized", style=":.b", pltargs...)
        # end
    end

    return new_x, new_y
end

"""
`rediscretize_concavecontour(ys::Vector, zs::Vector; closed_contour=false,
nperiodic=1)`

Receives a concave contour (collection of points) and rediscretizes it based on
arc length.
Use `closed_contour=true` if the contour is closed or to close the contour.
`nperiodic` is the number of periods used to resemble periodic boundary
conditions in the internal splining process.
"""
function rediscretize_concavecontour(ys::AbstractVector, zs::AbstractVector,
                                        discretization::multidisctype;
                                        closed_contour::Bool=false,
                                        nperiodic::Int=1,
                                        verify_spline::Bool=true)

    # Catch case that it is a closed contour but not declared so
    if !closed_contour && ys[1]==ys[end] && zs[1]==zs[end]
        return rediscretize_concavecontour(ys, zs, discretization;
                                            closed_contour=true,
                                            nperiodic=nperiodic)
    end

    # Error cases
    @assert length(ys)==length(zs) ""*
        "`ys` and `zs` must have the same length;"*
        " got lengths $(length(ys)) and $(length(zs))"

    if closed_contour
        @assert ys[1]==ys[end] && zs[1]==zs[end] ""*
        "Expected closed contour, but end points do not match;"*
        " $((ys[1], zs[1])) != $((ys[end], zs[end]))"
    end

    # Calculate cylindrical coordinates (theta in range -+180deg)
    rhos = @. sqrt(ys^2 + zs^2)                 # radius
    thetas = @. atan(ys, zs)                    # (rad) azimuthal angle

    # # Check that the section is monotonic and injective in theta
    # for i in 2:length(thetas)
    #     tht2 = thetas[i]
    #     tht1 = thetas[i-1]
    #     if tht2 <= tht1
    #         if i==length(thetas) && closed_contour
    #             nothing
    #         else
    #             error("Contour is not azimuthally injective in cylindrical"*
    #                     " coordinates. Only concave contours are supported.")
    #         end
    #     end
    # end

    # Sort them by theta and remove repeated (force monotonic and injective)
    rhothetas = sort([rhotht for rhotht in zip(rhos, thetas)]; by = x->x[2])
    unique!(x->x[2], rhothetas)

    # Repeat points beyond -+180deg to resemble periodic boundary conditions in the spline
    org_rhothetas = rhothetas
    for n in 1:nperiodic
        rhothetas = vcat(
                            [[rho, theta - n*(2*pi)] for (rho, theta) in org_rhothetas],
                            rhothetas,
                            [[rho, theta + n*(2*pi)] for (rho, theta) in org_rhothetas]
                        )
    end

    # Split rho and theta
    rhos = [rho for (rho, tht) in rhothetas]
    thetas = [tht for (rho, tht) in rhothetas]

    # Generate Akima spline in cylindrical coordinates (rho as a function of theta)
    rho_spl = FLOWMath.Akima(thetas, rhos)

    # Derivative of rho(theta)
    drhodtht(tht) = FLOWMath.derivative(rho_spl, tht)

    # Arc length function
    integrand(tht) = sqrt(rho_spl(tht)^2 + drhodtht(tht)^2)
    arclength(thti, thtf) = QuadGK.quadgk(integrand, thti, thtf)[1]

    # Compute the arc length at each theta position
    arclengths = zeros(length(thetas))

    for si in 2:length(arclengths)
        arclengths[si] = arclength(thetas[si-1], thetas[si]) + arclengths[si-1]
    end

    # Normalize arc lengths to go from 0 to 1 in the original range of thetas
    index1 = nperiodic * length(org_rhothetas) + 1
    if closed_contour
        index2 = (nperiodic+1) * length(org_rhothetas) + 1  # <-- Add 1 to close the contour
    else
        index2 = (nperiodic+1) * length(org_rhothetas)
    end
    arclengths .-= arclengths[index1]
    arclengths ./= arclengths[index2]

    # Generate arc length positions to redescritize the section into
    new_arclengths = multidiscretize(identity, 0, 1, discretization)

    # Rediscretize section based on arc length positions
    new_rhos = FLOWMath.akima(arclengths, rhos, new_arclengths)
    new_thetas = FLOWMath.akima(arclengths, thetas, new_arclengths)

    # Error case: Akima spline returned a NaN point
    if !prod(isnan.(new_rhos) .== false) || !prod(isnan.(new_thetas) .== false)
        @show minimum(thetas)
        @show maximum(thetas)
        @show round.(thetas*180/pi, digits=2)
        @show hcat(round.(thetas*180/pi, digits=2), round.(rhos, digits=2))
        @show hcat(round.(new_thetas*180/pi, digits=2), round.(new_rhos, digits=2))
       error("Got NaN in Akima spline!")
    end

    new_ys = new_rhos.*sin.(new_thetas)
    new_zs = new_rhos.*cos.(new_thetas)

    # Verification plots
    if verify_spline
        fig = plt.figure(figsize=[7, 5]*1/2)
        ax = fig.gca()
        ax.set_aspect("equal")

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.plot(ys, zs, "o-k", label="Raw", alpha=0.8, markersize=8)
        ax.plot(new_ys, new_zs, "s-r", label="Spline", alpha=1.0, markersize=4, linewidth=1)

        ax.legend(loc="best", frameon=true, fontsize=8)
    end

    return new_ys, new_zs

end

function rediscretize_concavecontour(ys::AbstractVector, zs::AbstractVector,
                                        discretization::Int; optargs...)

    multidiscretization = [(1.0, discretization, 1.0, true)]

    return rediscretize_concavecontour(ys, zs, multidiscretization; optargs...)

end

function rediscretize_concavecontour(points::AbstractMatrix, args...; optargs...)

    ys = view(points, :, 1)
    zs = view(points, :, 2)

    new_ys, new_zs = rediscretize_concavecontour(ys, zs, args...; optargs...)

    return hcat(new_ys, new_zs)

end

"""
    Returns the y-coordinates of upper and lower surfaces of a 4-digit NACA
    airfoil. x input and y outputs are normalized by the chord.

    * First digit describing maximum camber as percentage of the chord.
    * Second digit describing the distance of maximum camber from the airfoil
            leading edge in tens of percents of the chord.
    * Last two digits describing maximum thickness of the airfoil as percent of
            the chord.

    https://en.wikipedia.org/wiki/NACA_airfoil#Four-digit_series
"""
function naca4digits(d1, d2, d34, xs)

    m = d1/100
    p = d2/10
    t = d34/100

    yt = 5*t * (0.2969*sqrt.(xs) .- 0.126*xs .- 0.3516*xs.^2 .+ 0.2843*xs.^3 .- 0.1015*xs.^4)

    if d1==0 && d2==0
        return xs, xs, yt, -yt

    elseif d1==0 || d2==0
        error("Invalid airfoil definition d1=$d1, d2=$d2")

    else

        yc = similar(xs)
        thetas = similar(xs)

        for (i,x) in enumerate(xs)
            if x<0 || x>1
                error("Got invalid x-coordinate: $x")
            elseif x <= p
                yc[i] = m/p^2 * (2*p*x - x^2)
                thetas[i] = atan(2*m/p^2 * (p-x))
            else
                yc[i] = m/(1-p)^2 * ((1-2*p) + 2*p*x - x^2)
                thetas[i] = atan(2*m/(1-p)^2 * (p-x))
            end
        end

        xU = xs .- yt.*sin.(thetas)
        xL = xs .+ yt.*sin.(thetas)
        yU = yc .+ yt.*cos.(thetas)
        yL = yc .- yt.*cos.(thetas)

        return xU, xL, yU, yL
    end
end


function naca4digits(d1::Int64, d2::Int64, d34::Int64, n::Int64, r::Float64;
                                                                    off::Real=0)
    xs = [float(x) for x in discretize(x->x, off, 1, n, r; central=true)]
    xs[end] = 1.0
    xU, xL, yU, yL = naca4digits(d1, d2, d34, xs)
    xU[end] = 1.0
    xL[end] = 1.0

    if off==0
        return vcat(reverse(xU), xL[2:end]), vcat(reverse(yU), yL[2:end])
    else
        return vcat(reverse(xU), xL), vcat(reverse(yU), yL)
    end
end
