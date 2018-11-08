################################################################################
# # DESCRIPTION
#     Python functions for the manipulation of Paraview.
# # AUTHORSHIP
#   * Author    : Eduardo J Alvarez
#   * Email     : Edo.AlvarezR@gmail.com
#   * Created   : Nov 2017
#   * License   : MIT License
################################################################################

import paraview.simple as prvs

def render_vtk(vtkfilename, outputfilename, processing=None):
    """
      It receives the name of a vtk file, it renders it and saves a png image of
      it. Give it a function `processing` to apply any processing to the
      paraview view.
    """
    ext = ".png"
    reader = prvs.OpenDataFile(vtkfilename)
    renderView1 = prvs.GetActiveViewOrCreate('RenderView')
    testvtkDisplay = prvs.Show(reader, renderView1)

    if processing!=None:
        processing(renderView1)

    prvs.WriteImage(outputfilename + (ext if ext not in outputfilename else ""))



################################################################################
#                   VIEWS
################################################################################
def view1(renderView1):
    renderView1.CameraPosition = [-34.15196513494814, -4.703195703640361, 19.236477879153906]
    renderView1.CameraFocalPoint = [4.453237056732174, 10.56249999999999, 1.0000000000000009]
    renderView1.CameraViewUp = [0.2962788554858759, 0.3240403764420813, 0.8984523772728602]
    renderView1.CameraParallelScale = 11.735587947447558
