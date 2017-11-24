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

  # Examples
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

  # Arguments
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
##### END OF MESHING ###########################################################






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
##### END OF ALGEBRA ###########################################################
