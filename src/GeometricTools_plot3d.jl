#=##############################################################################
# DESCRIPTION
    Functions for exporting geometries into PLOT3D xyz format
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Jul 2019
  * License   : MIT License
=###############################################################################


"""
Generate a structured PLOT3D mesh
"""
function generatePLOT3D(grid::GridExtentions, name::String; path="",
                        num=nothing, ext="xyz", trickwopwop=false)

    fname = joinpath(path, name*(num!=nothing ? ".$num" : "")*"."*ext)
    f = open(fname, "w")

    ndivs = get_ndivsnodes(grid)
    dims = grid.dims

    # imax jmax kmax
    if trickwopwop
        print(f, ndivs[1], " ", dims==3  ? ndivs[2]*ndivs[3] :  ndivs[2], " 1")
    else
        for (i,ndiv) in enumerate(ndivs)
            print(f, i!=1 ? " " : "", ndiv)
        end
        for i in 1:3-dims
            print(f, i!=1 ? " " : "", 1)
        end
    end
    print(f, "\n")

    # imax × jmax ×kmax floating point x coordinates
    for n in 1:grid.nnodes
        s = @sprintf "%8.6f" get_node(grid, n)[1]
        print(f, n!=1 ? (n%10==0 ? "\n" : " ") : "", s)
    end
    print(f, "\n")
    # imax × jmax ×kmax floating point y coordinates
    for n in 1:grid.nnodes
        s = @sprintf "%8.6f" get_node(grid, n)[2]
        print(f, n!=1 ? (n%10==0 ? "\n" : " ") : "", s)
    end
    print(f, "\n")
    # imax × jmax ×kmax floating point z coordinates
    for n in 1:grid.nnodes
        s = @sprintf "%8.6f" get_node(grid, n)[3]
        print(f, n!=1 ? (n%10==0 ? "\n" : " ") : "", s)
    end
    close(f)
end
