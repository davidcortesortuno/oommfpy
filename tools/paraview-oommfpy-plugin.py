import numpy as np
from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.util import numpy_support
from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid

import oommfpy

paraview_plugin_version = oommfpy.__version__
oommfpy_extensions = ["omf", "ovf", "ohf"]
oommfpy_input_filetypes = ["automatic"] + ["omf", "ovf", "ohf"]


@smproxy.reader(
    name="oommfpy reader",
    extensions=oommfpy_extensions,
    file_description="oommfpy-supported files",
    support_reload=False,
)
class OOMMFPyReader(VTKPythonAlgorithmBase):

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkRectilinearGrid"
        )
        self._filename = None
        self._file_format = None

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=oommfpy_extensions,
        file_description="oommfpy-supported files"
    )
    def SetFileName(self, filename):
        if self._filename != filename:
            self._filename = filename
            self.Modified()

    @smproperty.stringvector(name="StringInfo", information_only="1")
    def GetStrings(self):
        return oommfpy_input_filetypes

    # @smproperty.stringvector(name="FileFormat", number_of_elements="1")
    # @smdomain.xml(
    #     """
    #     <StringListDomain name="list">
    #         <RequiredProperties>
    #             <Property name="StringInfo" function="StringInfo"/>
    #         </RequiredProperties>
    #     </StringListDomain>
    #     """
    # )
    # def SetFileFormat(self, file_format):
    #     # Automatically deduce input format
    #     if file_format == "automatic":
    #         file_format = None

    #     if self._file_format != file_format:
    #         self._file_format = file_format
    #         self.Modified()

    def RequestData(self, request, inInfoVec, outInfoVec):
        output = vtkRectilinearGrid.GetData(outInfoVec)
        output.Initialize()

        # Use oommfpy to read the mesh
        if self._filename.endswith('omf'):
            mesh = oommfpy.MagnetisationData(self._filename)
        else:
            mesh = oommfpy.FieldData(self._filename)

        mesh.generate_coordinates(compute_vertex_grid=True)
        mesh.generate_field()

        # points, cells = mesh.coordinates, mesh.cells
        output.SetDimensions(mesh.nx + 1, mesh.ny + 1, mesh.nz + 1)
        output.SetExtent(0, mesh.nx, 0, mesh.ny, 0, mesh.nz)

        # Scale the coordinates in nm (somehow small values do not render
        # well in Paraview)
        output.SetXCoordinates(numpy_support.numpy_to_vtk(mesh.vertex_grid[0] * 1e9))
        output.SetYCoordinates(numpy_support.numpy_to_vtk(mesh.vertex_grid[1] * 1e9))
        output.SetZCoordinates(numpy_support.numpy_to_vtk(mesh.vertex_grid[2] * 1e9))

        # Cell data
        # Spin directions
        spinData = numpy_support.numpy_to_vtk(mesh.field.reshape(-1))
        spinData.SetNumberOfComponents(3)
        spinData.SetName('spin')
        output.GetCellData().AddArray(spinData)

        # Magnetisation
        magData = numpy_support.numpy_to_vtk(mesh.field_norm)
        magData.SetNumberOfComponents(1)
        magData.SetName('M')
        output.GetCellData().AddArray(magData)

        # TODO: - HSL colours for in plane orientation?
        #       - Skyrmion number density

        return 1
