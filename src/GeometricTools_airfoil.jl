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
