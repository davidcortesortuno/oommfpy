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
        # output = dsa.WrapDataObject(vtkRectilinearGrid.GetData(outInfoVec))
        output = vtkRectilinearGrid.GetData(outInfoVec)
        output.Initialize()

        # Use oommfpy to read the mesh
        if self._filename.endswith('omf'):
            mesh = oommfpy.MagnetisationData(self._filename)
        else:
            mesh = oommfpy.FieldData(self._filename)

        mesh.generate_coordinates()
        mesh.generate_field()

        # points, cells = mesh.coordinates, mesh.cells
        output.SetDimensions(mesh.nx + 1, mesh.ny + 1, mesh.nz + 1)

        # output.SetXCoordinates(mesh.grid[0])
        # output.SetYCoordinates(mesh.grid[1])
        # output.SetZCoordinates(mesh.grid[2])

        output.SetXCoordinates(numpy_support.numpy_to_vtk(mesh.grid[0]))
        output.SetYCoordinates(numpy_support.numpy_to_vtk(mesh.grid[1]))
        output.SetZCoordinates(numpy_support.numpy_to_vtk(mesh.grid[2]))

        # Cell data
        # output.CellData.append(mesh.field.reshape(-1), 'M')
        spinData = numpy_support.numpy_to_vtk(mesh.field.reshape(-1))
        spinData.SetNumberOfComponents(3)
        spinData.SetName('spin')
        output.GetCellData().AddArray(spinData)

        magData = numpy_support.numpy_to_vtk(mesh.field_norm)
        magData.SetNumberOfComponents(1)
        magData.SetName('M')
        output.GetCellData().AddArray(magData)

        return 1


# @smproxy.writer(
#     name="oommfpy Writer",
#     extensions=oommfpy_extensions,
#     file_description="oommfpy-supported files",
#     support_reload=False,
# )
# @smproperty.input(name="Input", port_index=0)
# @smdomain.datatype(dataTypes=["vtkRectilinearGrid"], composite_data_supported=False)
# class MeshioWriter(VTKPythonAlgorithmBase):
#     def __init__(self):
#         VTKPythonAlgorithmBase.__init__(
#             self, nInputPorts=1, nOutputPorts=0, inputType="vtkRectilinearGrid"
#         )
#         self._filename = None
#
#     @smproperty.stringvector(name="FileName", panel_visibility="never")
#     @smdomain.filelist()
#     def SetFileName(self, filename):
#         if self._filename != filename:
#             self._filename = filename
#             self.Modified()
#
#     def RequestData(self, request, inInfoVec, outInfoVec):
#         mesh = dsa.WrapDataObject(vtkRectilinearGrid.GetData(inInfoVec[0]))
#
#         # Read points
#         points = np.asarray(mesh.GetPoints())
#
#         # Read cells
#         # Adapted from test/legacy_reader.py
#         cell_conn = mesh.GetCells()
#         cell_offsets = mesh.GetCellLocations()
#         cell_types = mesh.GetCellTypes()
#         cells_dict = {}
#         for vtk_cell_type in np.unique(cell_types):
#             offsets = cell_offsets[cell_types == vtk_cell_type]
#             ncells = len(offsets)
#             npoints = cell_conn[offsets[0]]
#             array = np.empty((ncells, npoints), dtype=int)
#             for i in range(npoints):
#                 array[:, i] = cell_conn[offsets + i + 1]
#             cells_dict[vtk_to_oommfpy_type[vtk_cell_type]] = array
#         cells = [oommfpy.CellBlock(key, cells_dict[key]) for key in cells_dict]
#
#         # Read point and field data
#         # Adapted from test/legacy_reader.py
#         def _read_data(data):
#             out = {}
#             for i in range(data.VTKObject.GetNumberOfArrays()):
#                 name = data.VTKObject.GetArrayName(i)
#                 array = np.asarray(data.GetArray(i))
#                 out[name] = array
#             return out
#
#         point_data = _read_data(mesh.GetPointData())
#         field_data = _read_data(mesh.GetFieldData())
#
#         # Read cell data
#         cell_data_flattened = _read_data(mesh.GetCellData())
#         cell_data = {}
#         for name, array in cell_data_flattened.items():
#             cell_data[name] = []
#             for cell_type in cells_dict:
#                 vtk_cell_type = oommfpy_to_vtk_type[cell_type]
#                 mask_cell_type = cell_types == vtk_cell_type
#                 cell_data[name].append(array[mask_cell_type])
#
#         # Use oommfpy to write mesh
#         oommfpy.write_points_cells(
#             self._filename,
#             points,
#             cells,
#             point_data=point_data,
#             cell_data=cell_data,
#             field_data=field_data,
#         )
#         return 1
#
#     def Write(self):
#         self.Modified()
#         self.Update()
