#=##############################################################################
# DESCRIPTION
    Functions for discretizating 1D, 2D, and 3D entities.

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jan 2024
  * License   : MIT License
=###############################################################################

"""
`rediscretize_concavecontour(ys::Vector, zs::Vector, discretization;
closed_contour=false, nperiodic=1, verify_spline=true)`

Receives a concave contour (collection of points) and rediscretizes it based on
arc length, where `discretization` can be a "multidiscretize" parameter set or
and integer with the number of evenly-spaced point to discretize the contour
into (use "multidiscretize" parameters for a more general discretization).

Use `closed_contour=true` if the contour is closed or to close the contour.
`nperiodic` is the number of periods used to resemble periodic boundary
conditions in the internal splining process.
"""
function rediscretize_concavecontour(ys::AbstractVector, zs::AbstractVector,
                                        discretization::multidisctype;
                                        closed_contour::Bool=false,
                                        nperiodic::Int=1,
                                        verify_spline::Bool=true,
                                        plot_title=nothing)

    # Catch case that it is a closed contour but not declared so
    if !closed_contour && ys[1]==ys[end] && zs[1]==zs[end]
        return rediscretize_concavecontour(ys, zs, discretization;
                                            closed_contour=true,
                                            nperiodic=nperiodic,
                                            plot_title=plot_title)
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
        fig = isnothing(plot_title) ?   plt.figure(figsize=[7, 5]*1/2) :
                                        plt.figure(plot_title, figsize=[7, 5]*1/2)
        ax = fig.gca()
        if !isnothing(plot_title); ax.set_title(plot_title); end;
        ax.set_aspect("equal")

        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)

        ax.plot(ys, zs, "o-k", label="Raw", alpha=0.8, markersize=8)
        ax.plot(new_ys, new_zs, "s-r", label="Spline",
                                alpha=1.0, markersize=4, linewidth=1)

        ax.legend(loc="best", frameon=true, fontsize=8)
    end

    return new_ys, new_zs

end

function rediscretize_concavecontour(ys::AbstractVector, zs::AbstractVector,
                                     discretization::Int, args...; optargs...)

    multidiscretization = [(1.0, discretization, 1.0, true)]

    return rediscretize_concavecontour(ys, zs, multidiscretization, args...;
                                                                    optargs...)

end

function rediscretize_concavecontour(points::AbstractMatrix, args...;
                                                                    optargs...)

    ys = view(points, :, 1)
    zs = view(points, :, 2)

    new_ys, new_zs = rediscretize_concavecontour(ys, zs, args...; optargs...)

    return hcat(new_ys, new_zs)

end



"""

`points` is an nxm matrix, where n is the dimensionality and m is the number of
points.

`parameterization` is a differentiable function `X -> s` for parameterizing the
line with a single parameter `s`. This function must be injective
(meaning that for each `X` there is a unique value of `s` and for each value of
`s` there is a unique `X`) and monotonically increasing along the line.

NOTE: The current implementation only supports using one of the coordinates
as the parameterization, otherwise don't expect to obtain sensical results.
"""
function rediscretize_line(points::AbstractMatrix, discretization::multidisctype;
                                parameterization::Function = X->X[1])

    # Obtain the parameter s corresponding to each point
    params = parameterization.(eachcol(points))

    # Check that the parameterization is strictly monotonic
    for i in 2:length(params)

        s2 = params[i]
        s1 = params[i-1]

        @assert s2 > s1 "Parameterization is not monotically increasing; s = $(params)"
    end

    # Spline each coordinate with respect to s
    X_spl = [FLOWMath.Akima(params, row) for row in eachrow(points)]


    # Derivate of X(s)
    dXds(s) = FLOWMath.derivative.(X_spl, s)

    # Arc length function
    integrand(s) = norm(dXds(s))
    arclength(si, sf) = QuadGK.quadgk(integrand, si, sf)[1]


    # Compute the arc length at each s
    arclengths = zeros(length(params))

    for si in 2:length(arclengths)
        arclengths[si] = arclength(params[si-1], params[si]) + arclengths[si-1]
    end

    # Normalize arc lengths to go from 0 to 1
    arclengths .-= arclengths[1]
    arclengths ./= arclengths[end]


    # Generate arc length positions to redescritize the section into
    new_arclengths = multidiscretize(identity, 0, 1, discretization)

    # Rediscretize line based on arc length positions
    new_points = [FLOWMath.akima(arclengths, row, new_arclengths) for row in eachrow(points)]

    # Convert coordinates of new points back to a matrix
    new_points = permutedims(hcat(new_points...))

    return new_points

end

function rediscretize_line(points::AbstractMatrix, discretization::Int, args...;
                                                                     optargs...)

    multidiscretization = [(1.0, discretization, 1.0, true)]

    return rediscretize_line(points, multidiscretization, args...; optargs...)

end

function rediscretize_line(points::AbstractVector{<:AbstractVector{<:Real}}, args...; optargs...)

    points_matrix = hcat(points...)

    points_matrix = rediscretize_line(points_matrix, args...; optargs...)

    return collect.(eachcol(points_matrix))

end
