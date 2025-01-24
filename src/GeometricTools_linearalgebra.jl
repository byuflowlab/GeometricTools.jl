#=##############################################################################
# DESCRIPTION
    Linear algebra tools

# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Oct 2021
  * License   : MIT License
=###############################################################################

dot(A, B) = sum(a*b for (a,b) in zip(A, B))
norm(A) = sqrt(mapreduce(x->x^2, +, A))

cross1(A, B) = A[2]*B[3] - A[3]*B[2]
cross2(A, B) = A[3]*B[1] - A[1]*B[3]
cross3(A, B) = A[1]*B[2] - A[2]*B[1]

function cross!(out, A, B)
    out[1] = cross1(A, B)
    out[2] = cross2(A, B)
    out[3] = cross3(A, B)
    return out
end

function cross(A::AbstractVector{T1}, B::AbstractVector{T2}) where {T1, T2}
    out = zeros(promote_type(T1, T2), 3)
    return cross!(out, A, B)
end
