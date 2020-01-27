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
  `parameterize(x, y, z; inj_var::Int64=1, s=0.0001, debug=false)`

Receives a contour (line) and returns a parametrization function f(s) of the
contour. The parametrization is done on the path of the splined
contour such that f(0.5) returns the point (x,y,z) where half of the entire
contour has been walked, and f(1.0) returns the last point in the contour. To
perform this parametrization, the contour must be injective in least one of the
variables (x, y, or z) (i.e., all values of the variable are unique), and that
variable must be given in increasing order.

  **Arguments**
  * `x::Array{Float64,1}`    : x-coordinates of the contour.
  * `y::Array{Float64,1}`    : y-coordinates of the contour.
  * `z::Array{Float64,1}`    : z-coordinates of the contour.

  **Optional Arguments**
  * `inj_var::Int64`         : Indicates the variable that is injective, with
                                  1=x, 2=y, and 3=z.
  * `s::Float64`             : Spline smoothness.

See Parametrization section in documentation for more details.
"""
function parameterize(x, y, z; inj_var::Int64=1, s=0.0001, debug=false, atol=0,
                                            kspl="automatic", bc="extrapolate")
  # ERROR CASES
  if size(x)!=size(y)!=size(z)  # Invalid contour
    error("Invalid contour. "*
              "$(size(x))!=$(size(y))!=$(size(z)) (size(x)!=size(y)!=size(z))")
  elseif !(inj_var in [1,2,3])  # Invalid input
    error("Invalid `inj_var`=$inj_var (!(inj_var in [1,2,3]))")
  end

  # Identifies injective and dependant variables
  inj = nothing   # Injective variable
  dep = []        # Dependant variables
  for (i,var) in enumerate([x, y, z])
    if i==inj_var
      inj = var
    else
      push!(dep, var)
    end
  end

  # Splines
  # bc="extrapolate"           # Out of boundary case
  # s=0.0001                 # Spline smoothness
  k = kspl=="automatic" ? min(size(x)[1]-1, 5) : kspl  # Spline order
  spl = []                   # Spline of each variable respect the injective

  for var in dep
    this_spl = Dierckx.Spline1D(inj, var; k=k, bc=bc, s=s)
    push!(spl, this_spl)
  end

  # Defines the path function
  dfdx1(x) = Dierckx.derivative(spl[1], x)    # Derivative of f respect x1
  dfdx2(x) = Dierckx.derivative(spl[2], x)    # Derivative of f respect x2
  fun(x) = sqrt.(1 .+ (dfdx1(x)).^2 .+ (dfdx2(x)).^2)   # Integrand
                                                  # Integral between xmin and x
  fun_s(this_x) = QuadGK.quadgk(fun, inj[1], this_x; atol=atol)[1]

  # Defines the normalized path function
  stot = fun_s(inj[end])      # Total length of the path
  norm_s(x) = fun_s(x)/stot         # Normalized path function

  # Defines the inverse normalized function
  function inv_norm_s(sstar; debug=false)
      # ERROR CASES
      if sstar<-0.001 || sstar>1.001 # Outside of domain
          error("Invalid normalized path length $(sstar) (sstar<0 || sstar>1)")
      end

      # Finds the x that matches the target (sstar)
      this_fun(x) = sstar - norm_s(x)   # Zeroed function
      bracket = [ inj[1]*(1-0.01*sign(inj[1])),# Bracket of the zero
                  inj[end]*(1+0.01*sign(inj[end]))]

      if debug
          println("sstar=$sstar\tbracket=$bracket")
          println("flow=$(this_fun(bracket[1]))")
          println("fup=$(this_fun(bracket[2]))")
      end

      # Zero
      this_x = Roots.find_zero(this_fun, bracket, Roots.Bisection())

      return this_x
  end

  # Redefine the original function in terms of the new parametrization s*
  function fstar(sstar)
      x = inv_norm_s(sstar; debug=debug)       # Finds the x that matches the path
      out = []                    # Output
      dep_i = 1                   # Next spline to evaluate
      for i in 1:3
        if i==inj_var
          push!(out, x)           # Injective variable
        else
          push!(out, spl[dep_i](x))   # Dependant variables
          dep_i += 1
        end
      end
      return out
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
function transform(V::Array{T1,1}, M::Array{T2,2}, T::Array{T3,1}
                                        ) where{T1<:Real, T2<:Real, T3<:Real}
    return [
                M[1,1]*(V[1]-T[1]) + M[1,2]*(V[2]-T[2]) + M[1,3]*(V[3]-T[3]),
                M[2,1]*(V[1]-T[1]) + M[2,2]*(V[2]-T[2]) + M[2,3]*(V[3]-T[3]),
                M[3,1]*(V[1]-T[1]) + M[3,2]*(V[2]-T[2]) + M[3,3]*(V[3]-T[3])
            ]
end

function transform!(out::Array{T1, 1}, V::Array{T2,1},
                            M::Array{T3,2}, T::Array{T4,1}
                           ) where{T1<:Real, T2<:Real, T3<:Real, T4<:Real}
    out[1] = M[1,1]*(V[1]-T[1]) + M[1,2]*(V[2]-T[2]) + M[1,3]*(V[3]-T[3])
    out[2] = M[2,1]*(V[1]-T[1]) + M[2,2]*(V[2]-T[2]) + M[2,3]*(V[3]-T[3])
    out[3] = M[3,1]*(V[1]-T[1]) + M[3,2]*(V[2]-T[2]) + M[3,3]*(V[3]-T[3])
end

function transform(Vs::AbstractArray{Array{P,1},1}, M::AbstractArray{P,2}, T::AbstractArray{P,1}
                                                                ) where{P<:Real}
  out = Array{Float64,1}[]
  for V in Vs
    push!(out, transform(V, M, T))
  end
  return out
end

"""
  `countertransform(Vp::Array{Float64,1}, invM::Array{Float64,2},
T::Array{Float64,1})`

Rotates and translates back a vector V' that had been rotated and translated
into the system (i', j', k') with origin T, and returns the original
``V=M^{-1}V' + T``. To ease repetitive computation, instead of giving the unit
vectors, give the inverse of their matrix.
"""
function countertransform(Vp::Array{T1,1}, invM::Array{T2,2}, T::Array{T3,1}
                                        ) where{T1<:Real, T2<:Real, T3<:Real}
    return [
                invM[1,1]*Vp[1] + invM[1,2]*Vp[2] + invM[1,3]*Vp[3] + T[1],
                invM[2,1]*Vp[1] + invM[2,2]*Vp[2] + invM[2,3]*Vp[3] + T[2],
                invM[3,1]*Vp[1] + invM[3,2]*Vp[2] + invM[3,3]*Vp[3] + T[3]
            ]
end

function countertransform!(out::Array{T1, 1}, Vp::Array{T2,1},
                            invM::Array{T3,2}, T::Array{T4,1}
                           ) where{T1<:Real, T2<:Real, T3<:Real, T4<:Real}
   out[1] = invM[1,1]*Vp[1] + invM[1,2]*Vp[2] + invM[1,3]*Vp[3] + T[1]
   out[2] = invM[2,1]*Vp[1] + invM[2,2]*Vp[2] + invM[2,3]*Vp[3] + T[2]
   out[3] = invM[3,1]*Vp[1] + invM[3,2]*Vp[2] + invM[3,3]*Vp[3] + T[3]
end

function countertransform(Vps::AbstractArray{Array{P,1},1}, invM::AbstractArray{P,2},
                                                  T::AbstractArray{P,1}) where{P<:Real}
  out = Array{Float64,1}[]
  for Vp in Vps
    push!(out, countertransform(Vp, invM, T))
  end
  return out
end

"""
  `check_coord_sys(M::Array{Float64,2}; raise_error::Bool=true)`

Checks that the unit vectors given as the matrix `M=[i;j;k]` define a coordinate
system
"""
function check_coord_sys(M::Array{T,2}; raise_error::Bool=true) where{T<:Real}
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
        error("Non-ortogonal system $M")
      else
        return false
      end
    end
  end
  return true
end

function check_coord_sys(M::Array{Array{T,1},1}; raise_error::Bool=true
                                                                ) where{T<:Real}
  dims = 3
  newM = zeros(Float64, dims,dims)
  for i in 1:dims
    newM[i, :] = M[i]
  end
  return check_coord_sys(newM; raise_error=raise_error)
end

"""
  `axis_rotation(r::Array{Float64, 1}, angle_deg::Float64)`

Returns the transformation matrix of rotation around an arbitrary axis of unit
vector `r`
"""
function axis_rotation(r::Array{T, 1}, angle_deg::Real) where{T<:Real}
  ux, uy, uz = r
  C = cos(angle_deg*pi/180)
  S = sin(angle_deg*pi/180)
  t = 1 - C
  M = [t*ux^2+C t*ux*uy-S*uz t*ux*uz+S*uy;
        t*ux*uy+S*uz t*uy^2+C t*uy*uz-S*ux;
        t*ux*uz-S*uy t*uy*uz+S*ux t*uz^2+C]
  return M
end

"""
  `rotation_matrix(yaw::Real, pitch::Real, roll::Real)`

Receives yaw, pitch, and roll angles (in degrees) and returns the rotation
matrix corresponding to this rotation.
(see http://planning.cs.uiuc.edu/node102.html)

NOTE: Naming follows aircraft convention, with
* yaw:    rotation about z-axis.
* pitch:  rotation about y-axis.
* roll:   rotation about x-axis.
"""
function rotation_matrix(yaw::Real, pitch::Real, roll::Real)
  a, b, g = yaw*pi/180, pitch*pi/180, roll*pi/180
  Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
  Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
  Rx = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)]
  return Rz*Ry*Rx
end

"""
  `rotation_matrix2(roll::Real, pitch::Real, yaw::Real)`

Receives yaw, pitch, and roll angles (in degrees) and returns the rotation
matrix corresponding to this rotation.
(see http://planning.cs.uiuc.edu/node102.html)

NOTE: Naming follows aircraft convention, with
* roll:   rotation about x-axis.
* pitch:  rotation about y-axis.
* yaw:    rotation about z-axis.

**Examples**
```jldoctest
julia> M = gt.rotation_matrix(0, 0, 0)
3×3 Array{Float64,2}:
  1.0  0.0  0.0
  0.0  1.0  0.0
  0.0  0.0  1.0

julia> M = gt.rotation_matrix(90, 0, 0)
3×3 Array{Float64,2}:
  1.0  0.0   0.0
  0.0  0.0  -1.0
  0.0  1.0   0.0

julia> X = [0, 1, 0];
julia> Xp = M*X
3-element Array{Float64,1}:
  0.0
  0.0
  1.0

julia> M = gt.rotation_matrix(0, 90, 0)
3×3 Array{Float64,2}:
  0.0  0.0  1.0
  0.0  1.0  0.0
 -1.0  0.0  0.0

julia> X = [1, 0, 0];
julia> Xp = M*X
3-element Array{Float64,1}:
  0.0
  0.0
 -1.0


julia> M = gt.rotation_matrix(0, 0, 90)
3×3 Array{Float64,2}:
  0.0  -1.0  0.0
  1.0   0.0  0.0
  0.0   0.0  1.0

julia> X = [1, 0, 0];
julia> Xp = M*X
3-element Array{Float64,1}:
  0.0
  1.0
  0.0


julia> M = gt.rotation_matrix(0, 45, 0)
3×3 Array{Float64,2}:
  0.707  0.0  0.707
  0.0    1.0  0.0
 -0.707  0.0  0.707

julia> X = [1, 0, 0];
julia> Xp = M*X
3-element Array{Float64,1}:
  0.707
  0.0
 -0.707


julia> M = gt.rotation_matrix(0, 45, 45)
3×3 Array{Float64,2}:
  0.5    -0.707  0.5
  0.5     0.707  0.5
 -0.707   0.0    0.707

julia> X = [1, 0, 0];
julia> Xp = M*X
3-element Array{Float64,1}:
  0.5
  0.5
 -0.707
```
"""
function rotation_matrix2(roll::Real, pitch::Real, yaw::Real)
  return rotation_matrix(yaw, pitch, roll)
end
##### END OF ALGEBRA ###########################################################





################################################################################
# ORTHOGONAL SPACE TRANSFORMATIONS
################################################################################

# ORTHOGONAL SYSTEMS (See https://en.wikipedia.org/wiki/Orthogonal_coordinates)
function cylindrical3D(X)
    r, theta, z = X
    return [ r*cos(theta), r*sin(theta), z ]
end
function cylindrical2D(X)
    r, theta = X
    return [ r*cos(theta),  r*sin(theta)]
end
function spherical3D(X)
    r, theta, phi = X
    return [  r*sin(theta)*cos(phi),
              r*sin(theta)*sin(phi),
              r*cos(theta) ]
end
function parabolic3D(X)
    u,v,z = X
    return [ 1/2*(u^2-v^2), u*v, z ]
end
function paraboloidal3D(X)
    u,v,phi = X
    return [u*v*cos(phi), u*v*sin(phi), 1/2*(u^2-v^2)]
end
function elliptic3D(X; a=1)
    u,v,z = X
    return [a*cosh(u)*cos(v), a*sinh(u)*sin(v), z]
end
function prolate3D(X; a=1)
    xi,etha,phi = X
    return [a*sinh(xi)*sin(etha)*cos(phi),
            a*sinh(xi)*sin(etha)*sin(phi),
            a*cosh(xi)*cos(etha)]
end
function oblate3D(X; a=1)
    xi,etha,phi = X
    return [a*cosh(xi)*cos(etha)*cos(phi),
            a*cosh(xi)*cos(etha)*sin(phi),
            a*sinh(xi)*sin(etha)]
end
function bipolar3D(X; a=1)
    u,v,z = X
    return [a*sinh(v)/(cosh(v)-cos(u)),
            a*sin(u)/(cosh(v)-cos(u)),
            z]
end
function toroidal3D(X; a=1)
    u,v,phi = X
    return [a*sinh(v)*cos(phi)/( cosh(v)-cos(u) ),
            a*sinh(v)*sin(phi)/( cosh(v)-cos(u) ),
            a*sin(u)/( cosh(v)-cos(u) )]
end
function conical3D(X; b=2, c=1)
    r,mu,nu = X
    return [r*mu*nu/(b*c),
            r/b*sqrt( (mu^2-b^2)*(nu^2-b^2) / (b^2-c^2) ),
            r/c*sqrt( (mu^2-c^2)*(nu^2-c^2) / (c^2-b^2) )]
end
##### END OF SPACE TRANSFORMATIONS #############################################
