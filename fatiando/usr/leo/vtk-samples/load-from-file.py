from enthought.mayavi import mlab
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader

casca = VTKFileReader()
casca.initialize('casca.vtk')

mlab.pipeline.surface(casca)
mlab.show()
