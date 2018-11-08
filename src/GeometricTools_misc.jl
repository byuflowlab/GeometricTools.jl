#=##############################################################################
# DESCRIPTION
    Miscellanous methods.
# AUTHORSHIP
  * Author    : Eduardo J Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 2017
  * License   : MIT License
=###############################################################################


"""
  `create_path(save_path::String, prompt::Bool)`

Create folder `save_path`. `prompt` prompts the user if true.
"""
function create_path(save_path::String, prompt::Bool)
  # Checks if folder already exists
  if isdir(save_path)
    if prompt
      inp1 = ""
      opts1 = ["y", "n"]
      while false==(inp1 in opts1)
        print("\n\nFolder $save_path already exists. Remove? (y/n) ")
        inp1 = readline()[1:end]
      end
      if inp1=="y"
        run(`rm $save_path -rf`)
        println("\n")
      else
        return
      end
    else
      run(`rm $save_path -rf`)
    end
  end
  run(`mkdir $save_path`)
end
