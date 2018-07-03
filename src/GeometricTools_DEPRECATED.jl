# Place here all functions that are listed for deprication

function add_field(grid::Grid, field_name::String, field_type::String,
                                                                field_data)

  warn("`add_field` function with no entry type argument is depricated."*
          " Entry type will default to `node`")

  add_field(grid, field_name, field_type, field_data, "node")
end



function calculate_field(grid::Grid, f, field_name::String, field_type::String)
  warn("`calculate_field` function with no entry type argument is depricated."*
          " Entry type will default to `node`")

  calculate_field(grid, f, field_name, field_type, "node")
end
