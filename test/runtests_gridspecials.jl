using Test
import GeometricTools as gt

try
    verbose
catch
    global verbose = true
end

@testset verbose=verbose "isedge Tests" begin

    # --- isedge() ---
    if verbose
        println("Testing isedge()...")
        println("Generating 4 x 3 unit grid...")
    end

    Pmin = [0.0, 0.0, 0.0]
    Pmax = [1.0, 1.0, 0.0]
    n = [4, 3, 0]

    loop_dims = [0, 1, 2]
    dim_splits = [1, 2]

    # Truth table for isEdge
    edgeCells = Array{Int, 3}(undef, 
                              length(loop_dims), 
                              length(dim_splits), 
                              2*n[1]*n[2])
    edgeCellsEval = Array{Int, 3}(undef,
                                  length(loop_dims), 
                                  length(dim_splits), 
                                  2*n[1]*n[2])

    # Case dim_split = 1
    edgeCells[1,1,:] = [1,1,1,0,1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,1]
    edgeCells[2,1,:] = [1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1]
    edgeCells[3,1,:] = [0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0]

    # Case dim_split = 2
    edgeCells[1,2,:] = [1,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1]
    edgeCells[2,2,:] = [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
    edgeCells[3,2,:] = [0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0]

    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)

            for ci = 1:trigrid.ncells
                edgeCellsEval[loop_dim+1, dim_split, ci] = gt.isedge(trigrid, ci)
            end

            @test all(edgeCellsEval[loop_dim+1, dim_split, :] == 
                      edgeCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end
end
