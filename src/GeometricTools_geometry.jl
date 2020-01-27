#=##############################################################################
# DESCRIPTION
    Methods for geometric manipulation.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2017
  * License   : MIT License
=###############################################################################

"""
  `AbstractDiscretization`

Super type for discretization strategies used with `discretize` and `multidiscretize`.

**Available discretization strategies**
- `Uniform`
- `Cosine`
- `FirstOverLast`
- `FirstOverCentral`
"""
abstract type AbstractDiscretization end
Base.broadcastable(m::AbstractDiscretization) = Ref(m)

"""
  `Uniform()`

Uses uniform intervals over the discretization interval
"""
struct Uniform <: AbstractDiscretization end

"""
  `Cosine()`

Uses smaller intervals at the beginning and end of the range using `cos`
"""
struct Cosine <: AbstractDiscretization end

"""
  `FirstOverLast(ratio)`

Uses the specified `ratio` of the first over the last interval length to set
interval sizes
"""
struct FirstOverLast{T} <: AbstractDiscretization
    ratio::T
end

"""
  `FirstOverSpecified(ratio, transition)`

Uses the specified `ratio` of the first/last interval length over the interval length
at the normalized location `transition` (0 < `transition` < 1) to set interval sizes.
"""
struct FirstOverSpecified{T} <: AbstractDiscretization
    ratio::T
    location::T
end

################################################################################
# MESHING AND DISCRETIZATION
################################################################################
"""
  `discretize(start, stop, n, discretization::AbstractDiscretization)`

Discretizes the range from `start` to `stop` using n intervals with the
discretization strategy specified by `discretization`.

  **Examples**

  The following lines show the discretization of a semi-circumference perimeter
  into 100 intervals of uniform length:

  ```julia
    julia> f(x) = (x, sqrt(1-round(x,digits=8)^2), 0)  # Semi-circunference of radius 1
    julia> f.(discretize(-1, 1, 100, Uniform()))
  ```
"""
function discretize(start, stop, n, discretization::AbstractDiscretization)
    error("`discretize` not defined for discretization method $discretization")
end

discretize(start, stop, n, ::Uniform) = range(start, stop, length=n+1)

function discretize(start, stop, n, ::Cosine)
   s = (1 .- cos.(0:pi/n:pi))/2
   x = start .+ (stop-start)*s
   return x
end

function discretize(start, stop, n, discretization::FirstOverLast)
    if n <= 1
        return [start, stop]
    else
        r = discretization.ratio
        len = 2/(n*(r+1))*(1 .+ (r-1)/(n-1)*(0:n-1))
        s = vcat(0, cumsum(len))
        return start .+ (stop-start)*s
    end
end

function discretize(start, stop, n, discretization::FirstOverSpecified)
    r1 = discretization.ratio
    t = discretization.location
    n1 = max(round(Int, n*t), 1)
    n2 = max(round(Int, n*(1-t)), 1)
    if n1 + n2 == n
        s1 = discretize(0, t, n1, FirstOverLast(r1))
        r2 = 2/n2*(1-t)/(s1[end]-s1[end-1])-1
        s2 = discretize(t, 1, n2, FirstOverLast(r2))
        s = vcat(s1, s2[2:end])
    else
        n1 = ceil(n*t)
        n2 = ceil(n*(1-t))
        # first part
        t1 = t/(1 - r1/(n1*(r1+1)))
        s1 = discretize(0, t1, n1, FirstOverLast(r1))
        # second part
        t2 = t - (t1-t)
        r2 = 2/n2*(1-t2)/(s1[end]-s1[end-1])-1
        s2 = discretize(t2, 1, n2, FirstOverLast(r2))
        if n1 <= 0
          s = s2
        elseif n2 <= 0
          s = s1
        else
          s = vcat(s1[1:end-1], s2[2:end])
        end
    end
    return start .+ s*(stop-start)
end

"""
  `multidiscretize(start, stop, c, n, discretization)`

Discretizes the range from `start` to `stop` into multiple sections of refinement.

  ** Arguments **
  * `start`: Lower bound.
  * `stop`: Upper bound.
  * `c`: Vector of normalized lengths of each section (must sum to one)
  * `n`: Number of intervals in each section
  * `discretization`: Discretization strategy for each section

  **Examples**

  The following lines show the discretization of a semi-circumference perimeter
  into 90 intervals done in three sections of discretization:

  ```julia
    julia> f(theta) = (cos(theta), sin(theta), 0)
    julia> c = [1/3, 1/3, 1/3]
    julia> n = 30 # or [30, 30, 30]
    julia> discretization = FirstOverLast(1/8) # or fill(FirstOverLast(1/8), 3)
    julia> points = f.(multidiscretize(0, pi, n, discretization, ))
  ```
"""
function multidiscretize(start, stop, c, n, discretization)

    # ensure the sum of `c` will add to one
    c ./= sum(c)

    # get edges
    edges = cumsum(vcat(0, c))

    # discretize each section
    ss = discretize.(edges[1:end-1], edges[2:end], n, discretization)

    # stack each section
    s = vcat(0, [s[2:end] for s in ss]...)

    return start .+ s*(stop-start)
end

"""
  Receives a closed 2D contour and splits it up in upper and lower surfaces as
  divided by lowest and highest x values. It returns `(upper, lower)` with
  `upper=(x,y)` the points of the upper surface, ditto for `lower`. Both `upper`
  and `lower` are given in increasing order in x. This function is useful for
  splitting up a closed contour into two sections that are injective in x, and
  that can be received by `parameterize()`.
"""
function splitcontour(x, y)

    le = argmin(x)
    te = argmax(x)

    # Split and sort from leading edge to trailing edge
    if le > te
        x_sec1, y_sec1 = x[le:-1:1], y[le:-1:1]
        x_sec2, y_sec2 = x[le:end], y[le:end]
    else
        x_sec1, y_sec1 = x[1:te], y[1:te]
        x_sec2, y_sec2 = x[end:-1:te], y[end:-1:te]
    end

    # Determine upper and lower surfaces
    if sum(y_sec1)/length(y_sec1) > sum(y_sec2)/length(y_sec2)
        upper = (x_sec1, y_sec1)
        lower = (x_sec2, y_sec2)
    else
        upper = (x_sec2, y_sec2)
        lower = (x_sec1, y_sec1)
    end

    return upper, lower
end


##### END OF MESHING ###########################################################

################################################################################
# PARAMETRIZATION
################################################################################
"""
  `parameterize(x, y)`

Receives a contour (line) and returns a parametrization function f(s) of the
contour. The parametrization is done on the path of the splined
contour such that f(0.5) returns the point (x,y,z) where half of the entire
contour has been walked, and f(1.0) returns the last point in the contour.

  **Arguments**
  * `x`    : x-coordinates of the contour.
  * `y`    : y-coordinates of the contour.

"""
function parameterize(x, y)

    # get arc length parameterization
    s = cumsum(vcat(0, sqrt.((x[2:end].-x[1:end]-1).^2 + (y[2:end].-y[1:end-1]).^2))

    # normalize to vary from 0 to 1
    s .= s./s[end]

    # now spline the function wrt x, y, and z
    xspl = Akima(s, x)
    yspl = Akima(s, y)

    # return splined function
    fstar = let xspl=xspl, yspl=yspl
        (t) -> xspl(t), yspl(t)
    end

    return fstar
end

"""
  `parameterize(x, y, z)`

Receives a contour (line) and returns a parametrization function f(s) of the
contour. The parametrization is done on the path of the splined
contour such that f(0.5) returns the point (x,y,z) where half of the entire
contour has been walked, and f(1.0) returns the last point in the contour.

  **Arguments**
  * `x`    : x-coordinates of the contour.
  * `y`    : y-coordinates of the contour.
  * `z`    : z-coordinates of the contour.

"""
function parameterize(x, y, z)

    # get arc length parameterization
    s = cumsum(vcat(0, sqrt.((x[2:end].-x[1:end]-1).^2 + (y[2:end].-y[1:end-1]).^2 + (z[2:end].-z[1:end-1]).^2)))

    # normalize to vary from 0 to 1
    s .= s./s[end]

    # now spline the function wrt x, y, and z
    xspl = Akima(s, x)
    yspl = Akima(s, y)
    zspl =  Akima(s, z)

    # return splined function
    fstar = let xspl=xspl, yspl=yspl, zspl=zspl
        (t) -> xspl(t), yspl(t), zspl(t)
    end

    return fstar
end

##### END OF PARAMETRIZATION ###################################################


################################################################################
# LINEAR ALGEBRA TRANSFORMATIONS
################################################################################
"""
  `transform(V::Array{Float64,1}, M::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates the vector V. Receives the i', j', k' unit vectors of an
euclidean system with origin T, and returns ``V'=M(V-T)`` (In this version, the
unit vectors have been organized as a matrix M=[i'; j'; k']).
"""
transform(V::AbstractVector, M::AbstractMatrix, T::AbstractVector) = M*(V-T)

"""
  `countertransform(Vp::Array{Float64,1}, invM::Array{Float64,2},
T::Array{Float64,1})`

Rotates and translates back a vector V' that had been rotated and translated
into the system (i', j', k') with origin T, and returns the original
``V=M^{-1}V' + T``. To ease repetitive computation, instead of giving the unit
vectors, give the inverse of their matrix.
"""
countertransform(Vp, invM, T) = invM*Vp+T

"""
  `check_coord_sys(M::Array{Float64,2}; raise_error::Bool=true)`

Checks that the unit vectors given as the matrix `M=[i;j;k]` define a coordinate
system
"""
function check_coord_sys(M::AbsractMatrix; raise_error=true)
    # Checks normalization
    for i in 1:size(M)[1]
        if abs(norm(M[i,:])-1) > 0.00000001
            println(M)
            if raise_error
                error("Not unitary axis: $(M[i,:])")
            else
                return false
            end
        end
    end

    # Checks ortogonality
    for i in size(M)[1]
        xi = M[i, :]
        xip1 = M[(i%size(M)[1])+1, :]
        proj = abs(dot(xi, xip1))
        if proj>0.00000001
            if raise_error
                error("Non-orthogonal system $M")
            else
                return false
            end
        end
    end
    return true
end

function check_coord_sys(M::AbstractVector{<:AbstractVector}; raise_error=true)
    dims = 3
    newM = zeros(Float64, dims,dims)
    for i in 1:dims
        newM[i, :] = M[i]
    end
    return check_coord_sys(newM; raise_error=raise_error)
end

"""
  `axis_rotation(r, angle)`

Returns the transformation matrix of rotation around an arbitrary axis of unit vector `r`
"""
function axis_rotation(r, angle)
  ux, uy, uz = r
  c = cos(angle)
  s = sin(angle)
  t = 1 - c
  M =  @SMatrix [t*ux^2+c     t*ux*uy-s*uz t*ux*uz+s*uy;
                 t*ux*uy+s*uz     t*uy^2+c t*uy*uz-s*ux;
                 t*ux*uz-s*uy t*uy*uz+s*ux     t*uz^2+c]
  return M
end

"""
  `rotation_matrix(yaw, pitch, roll)`

Receives yaw, pitch, and roll angles (in radians) and returns the rotation
matrix corresponding to this rotation.
(see http://planning.cs.uiuc.edu/node102.html)

NOTE: Naming follows aircraft convention, with
* roll:   rotation about x-axis.
* pitch:  rotation about y-axis.
* yaw:    rotation about z-axis.
"""
function rotation_matrix(roll, pitch, yaw)
  cr = cos(roll)
  sr = sin(roll)
  cp = cos(pitch)
  sp = sin(pitch)
  cy = cos(yaw)
  sy = sin(yaw)
  Rx = @SMatrix [1 0 0; 0 cr -sr; 0 sr cr]
  Ry = @SMatrix [cp 0 sp; 0 1 0; -sp 0 cp]
  Rz = @SMatrix [cy -sy 0; sy cy 0; 0 0 1]
  return Rz*Ry*Rx
end

################################################################################
# ORTHOGONAL SPACE TRANSFORMATIONS
################################################################################

# ORTHOGONAL SYSTEMS (See https://en.wikipedia.org/wiki/Orthogonal_coordinates)
cylindrical2D(r, theta) = @SVector [ r*cos(theta),  r*sin(theta)]
cylindrical3D(r, theta, z) = @SVector [ r*cos(theta), r*sin(theta), z ]
function spherical3D(r, theta, phi)
    st = sin(theta)
    ct = cos(theta)
    sp = sin(phi)
    cp = cos(phi)
    return @SVector [r*st*cp, r*st*sp, r*ct]
end
parabolic3D(u,v,z) = @SVector [ 1/2*(u^2-v^2), u*v, z ]
paraboloidal3D(u,v,phi) = @SVector [u*v*cos(phi), u*v*sin(phi), 1/2*(u^2-v^2)]
elliptic3D(u, v, z; a=1) = @SVector [a*cosh(u)*cos(v), a*sinh(u)*sin(v), z]
function prolate3D(xi, etha, phi; a=1)
    shx = sinh(xi)
    chx = cosh(xi)
    se = sin(etha)
    ce = cos(etha)
    sp = sin(phi)
    cp = cos(phi)
    return @SVector [a*shx*se*cp, a*shx*se*sp, a*chx*ce]
end
function oblate3D(xi, etha, phi; a=1)
    shx = sinh(xi)
    chx = cosh(xi)
    se = sin(etha)
    ce = cos(etha)
    sp = sin(phi)
    cp = cos(phi)
    return @SVector [a*chx*ce*cp, a*chx*ce*sp, a*shx*se]
end
function bipolar3D(u, v, z; a=1)
    su = sin(u)
    cu = cos(u)
    shv = sinh(v)
    chv = cosh(v)
    return @SVector [a*shv/(chv-cu), a*su/(chv-cu), z]
end
function toroidal3D(u, v, p; a=1)
    su = sin(u)
    cu = cos(u)
    shv = sinh(v)
    chv = cosh(v)
    sp = sin(phi)
    cp = cos(phi)
    return @SVector [a*shv*cp/( chv-cu ), a*shv*sp/( chv-cu ), a*su/( chv-cu )]
end
function conical3D(r, mu, nu; b=2, c=1)
    return @SVector [r*mu*nu/(b*c),
                     r/b*sqrt( (mu^2-b^2)*(nu^2-b^2) / (b^2-c^2) ),
                     r/c*sqrt( (mu^2-c^2)*(nu^2-c^2) / (c^2-b^2) )]
end
##### END OF SPACE TRANSFORMATIONS #############################################
