#=##############################################################################
# DESCRIPTION
    Methods for VTK Legacy formatting
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2017
  * License   : MIT License
=###############################################################################


"""
  Receives an array of lines, with lines[i] the i-th line, and lines[i][j] the
  j-th point in the i-th line, and formats them in vtk format. Give it point
  data through `values` in the same format, and it will output the data formated
  as corresponding to the vtk lines.

  ** Output **
  * (points, vtk_lines, vtk_values) where `points` is an array with all points
    in all lines, `lines` is an array of arrays containing the indices of
    every point in each line, and `vtk_values` is an array with the point data
    of each point in `points`.
"""
function lines2vtk(lines; values=nothing)
  points = Array{Float64,1}[]
  vtk_lines = Array{Int64,1}[]

  if values!=nothing
    vtk_values = []
  else
    vtk_values = nothing
  end

  n_points = 0
  for (i,line) in enumerate(lines)
    vtk_line = []
    for (j,point) in enumerate(line)
      push!(points, point)
      push!(vtk_line, n_points)
      if values!=nothing
        push!(vtk_values, values[i][j])
      end
      n_points += 1
    end
    push!(vtk_lines, vtk_line)
  end

  return (points, vtk_lines, vtk_values)
end


"""
  `generateCells(line1::Array{Array{Float64,1}},line2::Array{Array{Float64,1}};
                    point_data1, point_data2)`

Given two lines with the same amount of points, it generate cells by matching
the points between lines. For instance:

  ```julia
  julia> line1 = [p11,p12,p13,p14,p15,p16,p17]
  julia> line2 = [p21,p22,p23,p24,p25,p26,p27]
  julia> generateCells(line1, line2)
  ([cell1, cell2, ...], point_data)
  ```
where `cell1=[p11, p12, p22, p21]`, `cell2=[p12,p13,p23,p22]`, etc.

Give it point data corresponding to points in each line through `point_data1`
and `point_data2`, and it return it through `point_data` formatted for
`lines2vtk`'s `values` input.
"""
function generateCells(line1, line2; point_data1=nothing, point_data2=nothing)
  # Error case
  if size(line1)!=size(line2)
    error("Cells can't be generated from lines of unequal divisions!")
  end

  ncells = size(line1)[1]-1   # Number of cells
  cells = Array{Float64,1}[]  # Cells
                              # Point data on each cell
  point_data = (nothing in [point_data1, point_data2]) ? nothing : []

  # Builds cells
  for i in 1:ncells
    push!(cells, [line1[i], line1[i+1], line2[i+1], line2[i]])
    if point_data!=nothing
      push!(point_data, [point_data1[i], point_data1[i+1], point_data2[i+1],
                        point_data2[i]])
    end
  end

  return cells, point_data
end


"""
  `lines2vtkcells(line1::Array{Array{Float64,1}},line2::Array{Array{Float64,1}};
                    point_data1, point_data2)`

Given two lines with the same amount of points, it generate cells in VTK format
by matching the points between lines. For instance:

  ```julia
  julia> line1 = [p11,p12,p13]
  julia> line2 = [p21,p22,p23]
  julia> lines2vtkcells(line1, line2)
  (points, vtk_cells, point_data)
  ```
where `points=[p11,p12,p13,p21,p22,p23]`, and `vtk_cells=[[0,1,4,3],[1,2,5,4]]`.

Give it point data corresponding to points in each line through `point_data1`
and `point_data2`, and it return it through `point_data` formatted for
`generateVTK`'s `point_data` input.

Prefer this method over `generateCells()` since it will store the data
efficiently when generating the VTK.
"""
function lines2vtkcells(line1, line2; point_data1=nothing, point_data2=nothing)
  # Error case
  if size(line1)!=size(line2)
    error("Cells can't be generated from lines of unequal divisions!")
  end

  npoints = size(line1)[1]        # Number of points per line
  ncells = npoints-1              # Number of cells
  points = vcat(line1, line2)     # Points
  vtk_cells = Array{Int64,1}[]    # Cells
  # Point data
  if !(nothing in [point_data1, point_data2])
    if size(point_data1)!=size(point_data2)
      error("Invalid point data! "*
            "$(size(point_data1))!=$(size(point_data2))"*
            "(size(point_data1)!=size(point_data2))")
    else
      point_data = vcat(point_data1, point_data2)
    end
  else
    point_data = nothing
  end

  # Builds cells
  for i in 1:ncells
    point_index = i-1
    push!(vtk_cells, [point_index, point_index+1,
                        npoints+point_index+1, npoints+point_index])
  end

  return points, vtk_cells, point_data
end

"""
  `lines2vtkmulticells(line1, line2,
                          sections::Array{Tuple{Float64,Int64,Float64,Bool},1},
                          point_data1=nothing, point_data2=nothing)`

Like `lines2vtkcells()` it generates cells between two lines, but allows to
define multiple rows of cells in between. The rows are given by `sections`,
which is the same variable `sections` defined in `VTKtools_geometry.jl`'s
`multidiscretize()` function, which divides the space in between the lines.
For more details see docstring of `lines2vtkcells()` and `multidiscretize()`.
"""
function lines2vtkmulticells(line1, line2,
                          sections::Array{Tuple{Float64,Int64,Float64,Bool},1};
                          point_data1=nothing, point_data2=nothing)
  # ERROR CASES
  if size(line1)!=size(line2)
    error("Cells can't be generated from lines of unequal divisions!")
  end
  if !(nothing in [point_data1, point_data2])
    if size(point_data1)!=size(point_data2)
      error("Invalid point data! "*
            "$(size(point_data1))!=$(size(point_data2))"*
            "(size(point_data1)!=size(point_data2))")
    else
      point_data = []
    end
  else
    point_data = nothing
  end


  npoints = size(line1)[1]        # Number of points per line
  points = []                     # Points
  vtk_cells = Array{Int64,1}[]    # Cells

  # Discretize the space between the lines
  f(x) = x              # Function to discretize
  xlow, xhigh = 0, 1    # The space is parameterized between 0 and 1
  row_pos = multidiscretize(f, xlow, xhigh, sections) # Position of each row
                                                      # between the lines

  # Iterates over each row
  prev_line2 = line1
  prev_point_data2 = point_data1
  for (i,this_pos) in enumerate(row_pos[2:end])

    # Creates the upper and lower lines of this row
    this_line1 = prev_line2
    this_line2 = [line1[j] + this_pos*(line2[j]-line1[j]) for j in 1:npoints]

    # Linearly interpolates point data
      this_point_data1 = prev_point_data2
    if point_data!=nothing
      this_point_data2 = [point_data1[j] + this_pos*(point_data2[j]-point_data1[j]) for j in 1:npoints]
    else
      this_point_data2 = nothing
    end


    # Discretizes this row
    row = lines2vtkcells(this_line1, this_line2; point_data1=this_point_data1,
                                                  point_data2=this_point_data2)
    row_points, row_vtk_cells, row_point_data = row

    # Adds this row to the overall mesh
    if i==1 # Case of the first row of cells (adds both upper and lower bound)
      # Adds points
      points = vcat(points, row_points)
      # Adds point data
      if point_data!=nothing
        point_data = vcat(point_data, row_point_data)
      end
    else # Case of other rows (only adds upper bound)
      # Adds points
      points = vcat(points, row_points[npoints+1:end])
      # Adds point data
      if point_data!=nothing
        point_data = vcat(point_data, row_point_data[npoints+1:end])
      end
    end
    # Adds point indexing of cells
    vtk_cells = vcat(vtk_cells,
                            [cells+(i-1)*npoints for cells in row_vtk_cells])

    prev_line2 = this_line2
    prev_point_data2 = this_point_data2
  end

  return points, vtk_cells, point_data
end


"""
  `generateVTK(filename, points; lines, cells, point_data, path, num, time)`

Generates a vtk file with the given data.

  **Arguments**
  * `points::Array{Array{Float64,1},1}`     : Points to output.

  **Optional Arguments**
  * `lines::Array{Array{Int64,1},1}`  : line definitions. lines[i] contains the
                            indices of points in the i-th line.
  * `cells::Array{Array{Int64,1},1}`  : VTK polygons definiton. cells[i]
                            contains the indices of points in the i-th polygon.
  * `data::Array{Dict{String,Any},1}` : Collection of data point fields in the
                            following format:
                              data[i] = Dict(
                                "field_name" => field_name::String
                                "field_type" => "scalar" or "vector"
                                "field_data" => point_data
                              )
                            where point_data[i] is the data at the i-th point.

See `examples.jl` for an example on how to use this function.
"""
function generateVTK(filename::String, points;
                    lines::Array{Array{Int64,1},1}=Array{Int64,1}[],
                    cells::Array{Array{Int64,1},1}=Array{Int64,1}[],
                    point_data=nothing, num=nothing, time=nothing,
                    path="", comments="")
  aux = num!=nothing ? ".$num" : ""
  ext = aux*".vtk"
  if path !=""
    _path = string(path, (path[end]!="/" ? "/" : ""))
  else
    _path = path
  end
  f = open(string(_path, filename, ext), "w")

  # HEADER
  header = "# vtk DataFile Version 4.0" # File version and identifier
  header = string(header, "\n", " ", comments) # Title
  header = string(header, "\n", "ASCII") # File format
  header = string(header, "\n", "DATASET UNSTRUCTURED_GRID")
  write(f, header)

  # TIME
  if time!=nothing
    line0 = "\nFIELD FieldData 1"
    line1 = "\nSIM_TIME 1 1 double"
    line2 = "\n$(time)"
    write(f, line0*line1*line2)
  end

  np = size(points)[1]
  nl = size(lines)[1]
  nc = size(cells)[1]

  # POINTS
  write(f, string("\n", "POINTS ", np, " float"))
  for i in 1:np
    p = points[i]
    line = ""
    for (j,coor) in enumerate(p)
      line *= "$coor"
      if j!=size(p)[1]
        line *= " "
      end
    end
    write(f, "\n"*line)
  end

  # CELLS
  auxl = size(lines)[1]
  for line in lines
    auxl += size(line)[1]
  end
  auxc = size(cells)[1]
  for cell in cells
    auxc += size(cell)[1]
  end
  write(f, "\n\nCELLS $(np+nl+nc) $(2*np+auxl+auxc)")

  for i in 1:np+nl+nc
    if i<=np
      pts = [i-1]
    elseif i<=np+nl
      pts = lines[i-np]
    else
      pts = cells[i-(nl+np)]
    end
    line = "$(size(pts)[1])"
    for pt in pts
      line *=" $(pt)"
    end
    write(f, "\n"*line)
  end

  write(f, "\n\nCELL_TYPES $(np+nl+nc)")
  for i in 1:np+nl+nc
    if i<=np
      tpe = 1
    elseif i<=np+nl
      tpe = 4
    else
      tpe = 7
    end
    write(f, "\n"*"$tpe")
  end

  # POINT DATA
  if point_data!=nothing
      write(f, "\n\nPOINT_DATA $np")
  end
  _p_data = point_data!=nothing ? point_data : []
  for field in _p_data
    field_name = field["field_name"]
    field_type = field["field_type"]
    data = field["field_data"]
    if size(data)[1]!=np
      warn("Corrupted field $(field_name)! Field size != number of points.")
    end
    if field_type=="scalar"
      write(f, "\n\nSCALARS $field_name float\nLOOKUP_TABLE default")
      for entry in data
        write(f, "\n$entry")
      end
    elseif field_type=="vector"
      write(f, "\n\nVECTORS $field_name float")
      for entry in data
        line = ""
        for (j,coor) in enumerate(entry)
          line *= "$coor"
          if j!=size(entry)[1]
            line *= " "
          end
        end
        write(f, "\n"*line)
      end
    else
      error("Unknown field type $(field_type).")
    end
  end

  close(f)
end
