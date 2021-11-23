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
    support_reload=True,
)
class OOMMFPyReader(VTKPythonAlgorithmBase):

    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1,
            outputType="vtkRectilinearGrid")
        self._filename = None
        self._file_format = None

        self._fileList = []

        self._output = None

    # @smproperty.stringvector(name="FileName")
    # @smdomain.filelist()
    # @smhint.filechooser(
    #     extensions=oommfpy_extensions,
    #     file_description="oommfpy-supported files"
    # )
    # def SetFileName(self, filename):
    #     print(filename)
    #     if self._filename != filename:
    #         self._filename = filename
    #         self.Modified()

    # "StringInfo" and "String" demonstrate how one can add a selection widget
    # that lets user choose a string from the list of strings.
    @smproperty.stringvector(name="FileList", information_only="1")
    def GetStrings(self):
        return self._fileList

    @smproperty.stringvector(name="Files", number_of_elements="1")
    @smdomain.xml(
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="FileList" function="FileList"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def SetString(self, value):
        self._filename = value
        self.Modified()

    # Trying to read multiple files:
    @smproperty.xml("""
      <StringVectorProperty animateable="0"
                                clean_command="RemoveAllFileNames"
                                command="AddFileName"
                                name="FileName"
                                number_of_elements="1"
                                panel_visibility="never"
                                repeat_command="1">
            <FileListDomain name="files" />
          </StringVectorProperty>
    """)
    def AddFileName(self, filename):
        if self._filename != filename:
            self._filename = filename

            self._fileList.append(filename)

            self.Modified()

    @smproperty.xml("""
        <StringVectorProperty
                name="FileNames"
                clean_command="RemoveAllFileNames"
                command="AddFileName"
                animateable="0"
                number_of_elements="0"
                repeat_command="1">
                <FileListDomain name="files"/>
               <Documentation>
                 The list of files to be read by the reader. If more than 1
                 file is specified, the reader will switch to file series mode
                 in which it will pretend that it can support time and provide
                 1 file per time step.
               </Documentation>
             </StringVectorProperty>
    """)
    def RemoveAllFileNames(self):
        print('Cleaning file list')
        self._fileList.clear()

    # @smproperty.doublevector(name="TimestepValues",
    #                          information_only="1",
    #                          repeatable="1")
    # def GetTimestepValues(self):
    #     print(self._filename)

    def InitializeSystem(self, request, inInfoVec, outInfoVec):

        # Avoid reloading mesh if if was already loaded
        # (NOT WORKING)
        # if self._output is None:

        output = vtkRectilinearGrid.GetData(outInfoVec)
        output.Initialize()

        # Use oommfpy to read the mesh
        if self._filename.endswith('omf'):
            mesh = oommfpy.MagnetisationData(self._filename)
        else:
            mesh = oommfpy.FieldData(self._filename)

        mesh.generate_coordinates()
        print('Reloading field')
        mesh.generate_field()

        return output, mesh

    def RequestData(self, request, inInfoVec, outInfoVec):

        output, mesh = self.InitializeSystem(request, inInfoVec, outInfoVec)

        # output = self._output
        # mesh = self.mesh

        # points, cells = mesh.coordinates, mesh.cells
        output.SetDimensions(mesh.nx + 1, mesh.ny + 1, mesh.nz + 1)
        output.SetExtent(0, mesh.nx, 0, mesh.ny, 0, mesh.nz)

        # Scale the coordinates in nm (somehow small values do not render
        # well in Paraview)
        # TODO: scale spatial data from value in OVF file header
        output.SetXCoordinates(numpy_support.numpy_to_vtk(mesh.grid[0] * 1e9))
        output.SetYCoordinates(numpy_support.numpy_to_vtk(mesh.grid[1] * 1e9))
        output.SetZCoordinates(numpy_support.numpy_to_vtk(mesh.grid[2] * 1e9))

        # Cell data
        # Spin directions
        spinData = numpy_support.numpy_to_vtk(mesh.field.reshape(-1))
        spinData.SetNumberOfComponents(3)
        spinData.SetName('spin')
        # spinData._numpy_reference[:300] = 0.0

        # Add an array to the array list.
        # If an array with the same name already exists - then the added array
        # will replace it. Return the index of the added array.
        # (not working in Paraview)
        output.GetCellData().AddArray(spinData)

        # Magnetisation
        magData = numpy_support.numpy_to_vtk(mesh.field_norm)
        magData.SetNumberOfComponents(1)
        magData.SetName('M')
        output.GetCellData().AddArray(magData)

        # TODO: - HSL colours for in plane orientation?
        #       - Skyrmion number density

        return 1
