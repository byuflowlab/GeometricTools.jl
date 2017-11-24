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
