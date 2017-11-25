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
