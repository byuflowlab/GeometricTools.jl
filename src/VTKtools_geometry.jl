#=##############################################################################
# DESCRIPTION
    Methods for geometric manipulation.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2017
  * License   : MIT License
=###############################################################################




################################################################################
# MESHING AND DISCRETIZATION
################################################################################
"""
  `discretize(f, xlow, xhigh, n::Int64, r::Float64; central::Bool=false)`

Discretizes the continuous function `f` between the range `xlow` and `xhigh`
into `n` intervals, with `r` the ratio between first and last interval if
`central=false` or between first and central interval if `central=true`.

  **Examples**

  The following lines show the discretization of a semi-circumference perimeter
  into 100 intervals of uniform length:

  ```julia
    julia> f(x) = (x, sqrt(1-round(x,8)^2), 0)  # Semi-circunference of radius 1
    julia> discretize(f, -1, 1, 100, 1.0)
  ```
"""
function discretize(f, xlow, xhigh, n::Int64, r::Float64; central::Bool=false,
                      check::Bool=true)

  # ERROR CASES
  if n <= 0
    error("Invalid number of intervals (n <= 0)")
  elseif r <= 0
    error("Invalid expansion ratio (r <= 0)")
  end

  out = Any[f(xlow)]

  l = xhigh - xlow        # Length of domain to be discretized
  cumlen = 0              # Cumulative length already walked

  for i in 1:n

    # Case of uniform discretization
    if r==1.0
      len = l/n

    # Linear increment discretization (see notebook entry 20170519)
    # Case of no central expansion
    elseif !central
      p = l/( (n*(n-1)/2)*(r+1)/(r-1) )
      d1 = p*(n-1)/(r-1)
      len = d1 + p*(i-1)
      # println("i=$i\tp=$p\td1=$d1\tlen=$len")

    # Case of central expansion
    else
      _central = 0.5
      # Left of the center
      if i<=floor(n*_central)
        _l = l*_central
        _n = floor(n*_central)
        _r = r
        _i = i
      # Right of the center
      else
        _l = l*(1-_central)
        _n = n-floor(n*_central)
        _r = 1/r
        _i = i-floor(n*_central)
      end
      p = _l/( (_n*(_n-1)/2)*(_r+1)/(_r-1) )
      d1 = p*(_n-1)/(_r-1)
      len = d1 + p*(_i-1)
      # println("$i\t$_r\t$cumlen\t_n=$_n\t_i=$_i")
    end

    cumlen += len
    this_x = xlow + (cumlen/l)*(xhigh-xlow)
    # println("i=$i\tx=$this_x\tlen=$len")
    this_f = f(this_x)
    push!(out, this_f)
  end

  # Verifies correct discretization
  if check && abs((cumlen-l)/l)>0.0001
    error("Critical logic error! cumlen!=l ($cumlen!=$l)")
  end

  return out
end

"""
  `multidiscretize(f, xlow, xhigh, sections)`

Discretizes the continuous function `f` between the range `xlow` and `xhigh`
into multiple sections of refinement as specified in `sections`.

  ** Arguments **
  * `f`         : Continuous function of the form `f(x)` to be discretized
                  between `xlow` and `xhigh`,
  * `xlow`      : Lower bound.
  * `xhigh`     : Upper bound.
  * `sections`  : Array `[sec1, sec2, ...]`specifying the
                  sections of discretization in the format
                  `sec = (c::Float64, n::Int64, r::Float64, central::Bool)`,
                  with `c` the normalized length of this section (the sum of all
                  c must equal one), `n` the number of intervals in this section
                  , `r` the increment ratio between first and last interval if
                  `central=false` or between first and central interval if
                  `central=true`.

  **Examples**

  The following lines show the discretization of a semi-circumference perimeter
  into 90 intervals done in three sections of discretization:

  ```julia
    julia> f(theta) = (cos(theta), sin(theta), 0)
    julia> sec = (1/3, 30, 1/8, true)
    julia> points = multidiscretize(f, 0, pi, [sec, sec, sec])
  ```
"""
function multidiscretize(f, xlow, xhigh,
                          sections::Array{Tuple{Float64,Int64,Float64,Bool},1};
                          check::Bool=true)
  out = Any[f(xlow)]
  ctot = sum([sec[1] for sec in sections]) # Sum of all `c`s

  # Iterates over sections
  prev_xhigh = xlow
  for (c,n,r,central) in sections

    this_c = c/ctot # Makes sure that the sum of `c`s will add 1
    this_xlow = prev_xhigh
    this_xhigh = prev_xhigh + this_c*(xhigh-xlow)

    this_out = discretize(f, this_xlow, this_xhigh, n, r; central=central,
                                  check=check)
    out = vcat(out, this_out[2:end])

    prev_xhigh = this_xhigh
  end

  return out
end



"""
  Receives a closed 2D contour and splits it up in upper and lower surfaces as
  divided by lowest and highest x values. It returns `(upper, lower)` with
  `upper=(x,y)` the points of the upper surface, ditto for `lower`. Both `upper`
  and `lower` are given in increasing order in x. This function is useful for
  splitting up a closed contour into two sections that are injective in x, and
  that can be received by `parameterize()`.
"""
function splitcontour(x,y)

  # Flag indicating whether the contour start at the trailing or leading edge
  start_TE = x[1]==max(x)

  # Find the opposite end of the contour
  end_i = -1
  for (i, xi) in enumerate(x)
    if i==1
      nothing
    # Case of starting from the trailing edge
    elseif start_TE && xi > x[i-1]
      end_i = i-1
      break
    # Case of leading edge
    elseif !start_TE  && xi < x[i-1]
      end_i = i-1
      break
    end
  end

  # ERROR CASE
  if end_i==-1
    error("Logic error! End of contour not found!")
  end

  # Splits them up
  x_sec1, y_sec1 = x[1:end_i], y[1:end_i]
  x_sec2, y_sec2 = x[end_i:end], y[end_i:end]

  # Sorts them from LE to TE
  if x_sec1[1] > min(x); reverse!(x_sec1); reverse!(y_sec1); end;
  if x_sec2[1] > min(x); reverse!(x_sec2); reverse!(y_sec2); end;

  # Determines upper and lower surfaces
  if mean(y_sec1) > mean(y_sec2)
    upper = [x_sec1, y_sec1]
    lower = [x_sec2, y_sec2]
  else
    upper = [x_sec2, y_sec2]
    lower = [x_sec1, y_sec1]
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
function parameterize(x, y, z; inj_var::Int64=1, s=0.0001, debug=false,
                                                              kspl="automatic")
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
  bc="extrapolate"           # Out of boundary case
  # s=0.0001                 # Spline smoothness
  k = kspl=="automatic" ? min(size(x)[1]-1, 3) : kspl  # Spline order
  spl = []                   # Spline of each variable respect the injective
  for var in dep
    this_spl = Dierckx.Spline1D(inj, var; k=k, bc=bc, s=s)
    push!(spl, this_spl)
  end

  # Defines the path function
  dfdx1(x) = Dierckx.derivative(spl[1], x)    # Derivative of f respect x1
  dfdx2(x) = Dierckx.derivative(spl[2], x)    # Derivative of f respect x2
  fun(x) = sqrt.(1+(dfdx1(x)).^2+(dfdx2(x)).^2)   # Integrand
                                                  # Integral between xmin and x
  fun_s(this_x) = QuadGK.quadgk(fun, inj[1], this_x)[1]

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
  `transform(V::Array{Float64,1},
                      M::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates the vector V. Receives the i', j', k' unit vectors of an
euclidean system with origin T, and returns V'. (In this version, the unit
vectors have been organized as a matrix M)
"""
function transform(V::Array{Float64,1},
                    M::Array{Float64,2}, T::Array{Float64,1})
  return M*(V-T)
end

function transform(Vs::Array{Array{Float64,1},1},
                    M::Array{Float64,2}, T::Array{Float64,1})
  out = Array{Float64,1}[]
  for V in Vs
    push!(out, transform(V, M, T))
  end
  return out
end

"""
  `countertransform(Vp::Array{Float64,1},
                            invM::Array{Float64,2}, T::Array{Float64,1})`

Rotates and translates back a vector V' that had been rotated and translated
into the system (i', j', k') with origin T, and returns the original V.
To ease repetitive computation, instead of giving the unit vectors, give the
inverse of their matrix.
"""
function countertransform(Vp::Array{Float64,1},
                          invM::Array{Float64,2}, T::Array{Float64,1})
  return invM*Vp + T
end

function countertransform(Vps::Array{Array{Float64,1},1},
                          invM::Array{Float64,2}, T::Array{Float64,1})
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
function check_coord_sys(M::Array{Float64,2}; raise_error::Bool=true)
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

function check_coord_sys(M::Array{Array{Float64,1},1}; raise_error::Bool=true)
  dims = 3
  newM = zeros(dims,dims)
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
function axis_rotation(r::Array{Float64, 1}, angle_deg::Float64)
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
  Receives yaw, pitch, and roll angles (in degrees) and return the rotation
  matrix corresponding to this rotation.
  (see http://planning.cs.uiuc.edu/node102.html)
"""
function rotation_matrix(yaw::Real, pitch::Real, roll::Real)
  a, b, g = yaw*pi/180, pitch*pi/180, roll*pi/180
  Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]
  Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)]
  Rx = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)]
  return Rz*Ry*Rx
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
