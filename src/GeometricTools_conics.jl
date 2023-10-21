#=##############################################################################
# DESCRIPTION
    Methods for generation of conic sections
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Sep 2019
  * License   : MIT License
=###############################################################################


"""
    `line_intersect(x1, y1, x2, y2, x3, y3, x4, y4)`

    Returns the intersect between two lines, with the first line
passing through the points `(x1, y1), (x2, y2)`, and the second line passing
through the points `(x3, y3), (x4, y4)`.
"""
function line_intersect(x1, y1, x2, y2, x3, y3, x4, y4)
    Px = (
            (x1*y2-y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)
         ) / (
            (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
         )
    Py = (
            (x1*y2-y1*x2)*(y3-y4) - (y1-y2)*(x3*y4-y3*x4)
         ) / (
            (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
         )

    return Px, Py
end

"""
    `generate_conic_fun(Ax::Real, Ay::Real, Bx::Real, By::Real, Cx::Real, Cy::Real, Sx::Real, Sy::Real)`

     Returns a conic curve function with control point A, B, and C, and shoulder point S.
        See Raymer's Aircraft Design, lofting chapter (p. 128, 2nd ed.; p. 178, 5th ed.).
"""
function generate_conic_fun(Ax::Real, Ay::Real, Bx::Real, By::Real,
                            Cx::Real, Cy::Real, Sx::Real, Sy::Real)

    """
        Receives a value between 0 and 1, and returns the
    corresponding point along the conic curve, with 0 being A
    and 1 B.
    """
    function conic_fun(s)
        if s<0 || s>1
            error("Invalid conic parameter $s. Value between 0 and 1 expected.")
        end

        # D: Point in between A and B
        Dx = Ax + s*(Bx-Ax)
        Dy = Ay + s*(By-Ay)

        # E: Intersect between AS and CD
        Ex, Ey = line_intersect(Ax, Ay, Sx, Sy, Cx, Cy, Dx, Dy)

        # F: Intersect between BS and CD
        Fx, Fy = line_intersect(Bx, By, Sx, Sy, Cx, Cy, Dx, Dy)

        # P: Intersect between AF and BE
        Px, Py = line_intersect(Ax, Ay, Fx, Fy, Bx, By, Ex, Ey)

        return Px, Py
    end

    return conic_fun
end


"""
    `generate_conic_fun(Ax::Real, Ay::Real, Bx::Real, By::Real, Cx::Real, Cy::Real, rho::Real)`

    Returns a conic curve function with control point A, B, and C, and shape parameter rho.
        See Raymer's Aircraft Design, lofting chapter (p. 132, 2nd ed.;P. 182, 5th ed.).

    Hyperbola: rho > 0.5
    Parabola:  rho = 0.5
    Ellipse:   rho < 0.5
    Circle:    rho = 0.4142 and AC = BC
"""
function generate_conic_fun(Ax::Real, Ay::Real, Bx::Real, By::Real,
                            Cx::Real, Cy::Real, rho::Real)

    # D: Midpoint between A and B
    Dx = Ax + 0.5*(Bx-Ax)
    Dy = Ay + 0.5*(By-Ay)

    # S: Shoulder point in between D and C
    Sx = Dx + rho*(Cx-Dx)
    Sy = Dy + rho*(Cy-Dy)

    return generate_conic_fun(Ax, Ay, Bx, By, Cx, Cy, Sx, Sy)
end


"""
    `conic_cross_section(Ps::Array{Array{T1, 1},1}, CPs::Array{Array{T2, 1},1},
rhos::Array{T3, 1}, ss::Array{Array{T4, 1},1})`

    Receives a collection of points `Ps` along an open contour, and
stretching control points `CPs`, shape parameters `rhos`, and probing
parameters `ss` associated to every section of the contour, and it returns
a compound-conic cross section (it is a closed loop, so the number of
sections is equal to the number of points)
"""
function conic_cross_section(Ps::Array{Array{T1, 1},1},
                             CPs::Array{Array{T2, 1},1},
                             rhos::Array{T3, 1},
                             ss::Array{Array{T4, 1},1}
                            ) where {T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    nPs = size(Ps, 1)

    if size(CPs, 1)!=nPs
        error("Invalid CPs. Expected $nPs, got $(size(CPs, 1)).")
    elseif size(rhos, 1)!=nPs
        error("Invalid rhos. Expected $nPs, got $(size(rhos, 1)).")
    elseif size(ss, 1)!=nPs
        error("Invalid ss. Expected $nPs, got $(size(ss, 1)).")
    end

    points = []

    for i in 1:nPs
        A = Ps[i]
        B = Ps[i%nPs + 1]

        conic_fun = generate_conic_fun(Ps[i]..., Ps[i%nPs + 1]...,
                                                CPs[i]..., rhos[i])
        for s in ss[i]
            push!(points, [conic_fun(s)...])
        end
    end

    return points
end
