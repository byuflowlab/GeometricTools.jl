using Test
import GeometricTools as gt

try
    verbose
catch
    global verbose = true
end

@testset verbose=verbose "Gridspecials Tests" begin

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

    # --- isedge() ---
    if verbose
        println("Testing neighbor()...")
        println("Generating 4 x 3 unit grid...")
    end

    Pmin = [0.0, 0.0, 0.0]
    Pmax = [1.0, 1.0, 0.0]
    n = [4, 3, 0]

    loop_dims = [0, 1, 2]
    dim_splits = [1]

    # Case: inner cell
    ci = 11
    # Correct values for neigbor
    neighCells = Array{Int, 3}(undef, length(loop_dims), length(dim_splits), 3)
    neighCellsEval = Array{Int, 3}(undef, length(loop_dims), length(dim_splits), 3)

    neighCellsEval .= 0
    # dim_split = 1
    neighCells[1, 1, :] = [4, 14, 12]
    neighCells[2, 1, :] = [4, 14, 12]
    neighCells[3, 1, :] = [4, 14, 12]

    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # Case: cell on Xmin edge
    ci = 10

    # dim_split = 1
    neighCells[1, 1, :] = [17, 0, 9]
    neighCells[2, 1, :] = [17, 15, 9]
    neighCells[3, 1, :] = [17, 0, 9]

    neighCellsEval .= 0
    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                if !all(cidx == [0,0,0])
                    neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
                end
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # Case: cell on Xmax edge
    ci = 23

    # dim_split = 1
    neighCells[1, 1, :] = [16, 0, 24]
    neighCells[2, 1, :] = [16, 18, 24]
    neighCells[3, 1, :] = [16, 0, 24]

    neighCellsEval .= 0
    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                if !all(cidx == [0,0,0])
                    neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
                end
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # Case: cell on Ymin edge
    ci = 3

    # dim_split = 1
    neighCells[1, 1, :] = [0, 6, 4]
    neighCells[2, 1, :] = [0, 6, 4]
    neighCells[3, 1, :] = [20, 6, 4]

    neighCellsEval .= 0
    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                if !all(cidx == [0,0,0])
                    neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
                end
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # Case: cell on Ymax edge
    ci = 24

    # dim_split = 1
    neighCells[1, 1, :] = [0, 21, 23]
    neighCells[2, 1, :] = [0, 21, 23]
    neighCells[3, 1, :] = [7, 21, 23]

    neighCellsEval .= 0
    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                if !all(cidx == [0,0,0])
                    neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
                end
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # Case: cell on Xmin, Ymax edge
    ci = 18

    # dim_split = 1
    neighCells[1, 1, :] = [0, 0, 17]
    neighCells[2, 1, :] = [0, 23, 17]
    neighCells[3, 1, :] = [1, 0, 17]

    neighCellsEval .= 0
    for dim_split in dim_splits
        if verbose; print("Testing dim_split=$dim_split loop_dim="); end

        for loop_dim in loop_dims
            if verbose; print("$loop_dim, "); end

            grid = gt.Grid(Pmin, Pmax, n, loop_dim)
            trigrid = gt.GridTriangleSurface(grid, dim_split)
            ndivs = Tuple(collect(1:(d !=0 ? d : 1) for d in trigrid._ndivscells))
            linc = LinearIndices(ndivs)

            for ni = 1:3
                cidx = gt.neighbor(trigrid, ni, ci; preserveEdge=true)
                if !all(cidx == [0,0,0])
                    neighCellsEval[loop_dim+1, dim_split, ni] = linc[cidx...]
                end
            end

            @test all(neighCellsEval[loop_dim+1, dim_split, :] == 
                      neighCells[loop_dim+1, dim_split, :])
        end

        if verbose; println(""); end
    end

    # --- get_num_cells_around_node() ---
    if verbose
        println("Testing gt.get_num_cells_around_node()...")
        println("Generating 4 x 3 unit grid...")
    end

    P_min = [0, 0, 0]
    P_max = [2, 1, 0]
    n = [2, 3, 0]

    # loop_dim = 0 case
    loop_dim = 0
    rect_grid = gt.Grid(P_min, P_max, n, loop_dim)

    for dim_split in [1, 2]
        if verbose; println("Testing loop_dim=$loop_dim dim_split = $dim_split"); end
        tri_grid = gt.GridTriangleSurface(rect_grid, dim_split)

        cart = CartesianIndices((1:n[1]+1, 1:n[2]+1))

        # Corners
        @test gt.get_num_cells_around_node(tri_grid, cart[1]) == 2
        @test gt.get_num_cells_around_node(tri_grid, cart[3]) == 1
        @test gt.get_num_cells_around_node(tri_grid, cart[10]) == 1
        @test gt.get_num_cells_around_node(tri_grid, cart[12]) == 2

        # Edges
        @test gt.get_num_cells_around_node(tri_grid, cart[2]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[7]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[9]) == 3

        # Interior
        @test gt.get_num_cells_around_node(tri_grid, cart[8]) == 6
        @test gt.get_num_cells_around_node(tri_grid, cart[5]) == 6
    end

    # loop_dim = 1 case
    loop_dim = 1
    rect_grid = gt.Grid(P_min, P_max, n, loop_dim)

    for dim_split in [1, 2]
        if verbose; println("Testing loop_dim=$loop_dim dim_split = $dim_split"); end
        tri_grid = gt.GridTriangleSurface(rect_grid, dim_split)

        cart = CartesianIndices((1:n[1], 1:n[2]+1))

        # Edges
        @test gt.get_num_cells_around_node(tri_grid, cart[1]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[2]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[7]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[8]) == 3

        # Interior
        @test gt.get_num_cells_around_node(tri_grid, cart[5]) == 6
        @test gt.get_num_cells_around_node(tri_grid, cart[4]) == 6
    end

    # loop_dim = 2 case
    loop_dim = 2
    rect_grid = gt.Grid(P_min, P_max, n, loop_dim)

    for dim_split in [1, 2]
        if verbose; println("Testing loop_dim=$loop_dim dim_split = $dim_split"); end
        tri_grid = gt.GridTriangleSurface(rect_grid, dim_split)

        cart = CartesianIndices((1:n[1]+1, 1:n[2]))

        # Edges
        @test gt.get_num_cells_around_node(tri_grid, cart[1]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[3]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[4]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[6]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[7]) == 3
        @test gt.get_num_cells_around_node(tri_grid, cart[9]) == 3

        # Interior
        @test gt.get_num_cells_around_node(tri_grid, cart[1]) == 6
        @test gt.get_num_cells_around_node(tri_grid, cart[7]) == 6
    end
end
