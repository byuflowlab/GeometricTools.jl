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
  Receives an array of lines as lines[i] with the i-th line and lines[i][j] the
  j-th point in the i-th line, and formats them in vtk format.

  returns (points, vtk_lines)
    where `points` is an array with all points in all lines, and `lines` contains
    is an array of arrays containing the indexes of every point in each line.
"""
function lines2vtk(lines; values=nothing)
  points = []
  vtk_lines = []

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
    `generateVTK(filename, points [, lines, cells, point_data, line_n_cell_data,
                  path, comments])`
  Generates a vtk file with the given data.

  # Arguments
  *   `points`        :   ([float[]]) array of points.
  * (KEYWORDS)
  *   `lines`         :   ([float[]]) lines definition. lines[i] contains the
                            index of all the points in the i-th line.
  *   `cells`         :   ([float[]]) polygons definiton. cells[i] contains the
                            index of all the points in the i-th polygon
  *   `data`          :   ([Dict]) collection of data fields in the following
                            format:
                              data[i] = Dict(
                                "field_name" => field_name
                                "field_type" => "scalar" or "vector"
                                "field_data" => data
                              )
                            where data[i] is the data at the i-th point.
"""
function generateVTK(filename::String, points;
                    lines=nothing, cells=nothing,
                    point_data=nothing, line_n_cell_data=nothing,
                    path="", comments="", num=nothing)
  if num!=nothing
    nt = ".$(num)"
  else
    nt = ""
  end
  ext = "$(nt).vtk"
  if path !=""
    _path = string(path, (path[end]!="/" ? "/" : ""))
  else
    _path = ""
  end
  f = open(string(_path, filename, ext), "w")

  # HEADER
  header = "# vtk DataFile Version 4.0" # File version and identifier
  header = string(header, "\n", " ", comments) # Title
  header = string(header, "\n", "ASCII") # File format
  header = string(header, "\n", "DATASET UNSTRUCTURED_GRID")
  write(f, header)


  np = size(points)[1]
  nl = lines!=nothing ? size(lines)[1] : 0
  nc = cells!=nothing ? size(cells)[1] : 0

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
  if lines!=nothing
    auxl = size(lines)[1]
    for line in lines
      auxl += size(line)[1]
    end
  else
    auxl = 0
  end
  if cells!=nothing
    auxc = size(cells)[1]
    for cell in cells
      auxc += size(cell)[1]
    end
  else
    auxc = 0
  end
  write(f, "\n\nCELLS $(np+nl+nc) $(2*np+auxl+auxc)")

  # println("$np, $nl, $nc")
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

  # # CELL DATA
  # if line_n_cell_data!=nothing
  #   write(f, "\n\nCELL_DATA $(nl+nc)")
  #   _c_data = line_n_cell_data
  # else
  #   _c_data = []
  # end
  # for field in _c_data
  #   field_name = field["field_name"]
  #   field_type = field["field_type"]
  #   data = field["field_data"]
  #   if size(data)[1]!=nl+nc
  #     warn("Corrupted field $(field_name)!"*
  #           " Field size != number of lines + cells.")
  #   end
  #   if field_type=="scalar"
  #     write(f, "\n\nSCALARS $field_name float\nLOOKUP_TABLE default")
  #     for entry in data
  #       write(f, "\n$entry")
  #     end
  #   elseif field_type=="vector"
  #     write(f, "\n\nVECTORS $field_name float")
  #     for entry in data
  #       line = ""
  #       for (j,coor) in enumerate(entry)
  #         line *= "$coor"
  #         if j!=size(p)[1]
  #           line *= " "
  #         end
  #       end
  #       write(f, "\n"*line)
  #     end
  #   else
  #     error("Unknown field type $(field_type).")
  #   end
  # end
  close(f)
end
