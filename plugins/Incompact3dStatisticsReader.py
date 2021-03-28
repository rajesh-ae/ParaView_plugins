#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 19:22:23 2021

@author: Rajesh Venkatesan
"""

from paraview.util.vtkAlgorithm import *
from vtk.util import numpy_support
import numpy as np
import os
import re
#------------------------------------------------------------------------------
# Incompact3D Statistics Data Reader
#------------------------------------------------------------------------------
def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified

@smproxy.reader(name="Incompact3dStatisticsReader", label="Incompact3d Simulation Statistics Data Reader",
                extensions="i3d", file_description="Incompact3d files")
class Incompact3dStatisticsReader(VTKPythonAlgorithmBase):
    """ A reader that reads a Incompact3D simulation configuration file (.i3d) 
        along with solution statistics data (umean.dat##,vmean.dat##, etc.).

        Author: Rajesh Venkatesan
        e-mail: vrajesh.ae[at]gmail.com
    """
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkRectilinearGrid')
        self._filename = None
        self._ndata = None
        self._fileID = None
        self._nFillZeros = None
        self._domExtent = None
        self._precision = 0
        self._scalarFlag = 1

    # Get the data - Called during RequestData phase
    def _get_raw_data(self):
        if self._ndata is not None:
            return self._ndata

        if self._filename is None:
            raise RuntimeError("No filename specified")

        pathToSimulation = os.path.dirname(os.path.abspath(self._filename))

        nx = 0
        ny = 0
        nz = 0
        nscalar = 0

        # Read the parameters from i3d configuration file
        i3d_config = open(self._filename,'r')        
        for num, line in enumerate(i3d_config,1):
            if ('nx' in line) and (nx == 0):
                mynums = re.findall(r"\d+", line)
                nx = int(mynums[0])
            if ('ny' in line) and (ny == 0):
                mynums = re.findall(r"\d+", line)
                ny = int(mynums[0])
            if ('nz' in line) and (nz == 0):
                mynums = re.findall(r"\d+", line)
                nz = int(mynums[0])
            if 'xlx' in line:
                mynums = re.findall(r"\d*\.\d+|\d+", line)
                xlx = float(mynums[0])
            if 'yly' in line:
                mynums = re.findall(r"\d*\.\d+|\d+", line)
                yly = float(mynums[0])
            if 'zlz' in line:
                mynums = re.findall(r"\d*\.\d+|\d+", line)
                zlz = float(mynums[0])
            if 'istret' in line:
                mynums = re.findall(r"\d+", line)
                istret = int(mynums[0])
            if 'numscalar' in line:
                mynums = re.findall(r"\d+", line)
                nscalar = int(mynums[0])
            if 'nclx1' in line:
                mynums = re.findall(r"\d+", line)
                nclx = int(mynums[0])
            if 'ncly1' in line:
                mynums = re.findall(r"\d+", line)
                ncly = int(mynums[0])
            if 'nclz1' in line:
                mynums = re.findall(r"\d+", line)
                nclz = int(mynums[0])

        i3d_config.close()

        # Construct the grid - to be improved 
        if (nclx==0):
            dx = xlx/nx
        elif ((nclx==1) or (nclx==2)):
            dx = xlx/(nx-1.)
        if (ncly==0):
            dy = yly/ny
        elif ((ncly==1) or (ncly==2)):
            dy = yly/(ny-1.)
        if (nclz==0):
            dz = zlz/nz
        elif ((nclz==1) or (nclz==2)):
            dz = zlz/(nz-1.)

        x = np.zeros(nx)
        y = np.zeros(ny)
        z = np.zeros(nz)

        for i in range(nx):
            x[i] = i*dx
        
        if (istret != 0):
            y = np.loadtxt(pathToSimulation+'/yp.dat')
        else:
            for j in range(ny):
                y[j] = j*dx

        for k in range(nz):
            z[k] = k*dz

        # Total grid size
        nsize = nx*ny*nz

        # Set the precision to read in data files
        if (self._precision):
            file_precision = np.dtype(np.float64)
        else:
            file_precision = np.dtype(np.float32)

        # Create a list of all available statistics
        list_of_arrays = ['umean','vmean','wmean','pmean','uumean','vvmean','wwmean','uvmean','uwmean','vwmean','kmean']

        # Add scalars to the list if needed
        if (nscalar) and (self._scalarFlag):
            for iscalar in range(nscalar):
                this_scalar = 'phi'+str(iscalar+1).zfill(2)+'mean'
                list_of_arrays.append(this_scalar)
                this_scalar = 'phiphi'+str(iscalar+1).zfill(2)+'mean'
                list_of_arrays.append(this_scalar)
        
        # The data will be stored as vtkRectilinearGrid
        from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid
        self._ndata = vtkRectilinearGrid()

        # Set the grid parameters
        self._ndata.SetDimensions(nx,ny,nz);
        self._ndata.SetXCoordinates(numpy_support.numpy_to_vtk(x))
        self._ndata.SetYCoordinates(numpy_support.numpy_to_vtk(y))
        self._ndata.SetZCoordinates(numpy_support.numpy_to_vtk(z))
        point_data = self._ndata.GetPointData()

        # Try to read all statistics variables
        for this_array in list_of_arrays:
            try:
                curr_array = np.fromfile(pathToSimulation+"/"+this_array+".dat"+str(self._fileID).zfill(self._nFillZeros),dtype=file_precision,count=nsize)
                vtk_array = numpy_support.numpy_to_vtk(curr_array)
                vtk_array.SetName(this_array)
                point_data.AddArray(vtk_array)          
            except FileNotFoundError:
                # Skip the unavailable ones - move this to log message
                print('Statistics file',this_array+".dat"+str(self._fileID).zfill(self._nFillZeros),'does not exist! Skipping this file ...')

        return self._get_raw_data()   

    # Get the grid size/extent - called during RequestInformation phase
    def _get_grid_size(self):
        if self._filename is None:
            raise RuntimeError("No filename specified")

        # Read the grid parameters from i3d configuration file
        i3d_config = open(self._filename,'r')
        nx = 0
        ny = 0
        nz = 0

        i3d_config = open(self._filename,'r')        
        for num, line in enumerate(i3d_config,1):
            if ('nx' in line) and (nx == 0):
                mynums = re.findall(r"\d+", line)
                nx = int(mynums[0])
            if ('ny' in line) and (ny == 0):
                mynums = re.findall(r"\d+", line)
                ny = int(mynums[0])
            if ('nz' in line) and (nz == 0):
                mynums = re.findall(r"\d+", line)
                nz = int(mynums[0])
        i3d_config.close()
        
        # Save the extent
        self._domExtent = (0,nx-1,0,ny-1,0,nz-1)

        return

# Get i3d file name 
    @smproperty.stringvector(name="Config.FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="i3d", file_description="Incompact3d config files")
    def GetConfigFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self._ndata = None
            self.Modified()

# Get data file ID
    @smproperty.xml("""
        <IntVectorProperty name="File ID"
            command="GetFileID"
            number_of_elements="1"
            default_values="1">
            <IntRangeDomain name="range" />
            <Documentation>If the files are umean.dat0100, vmean.dat0100, etc., enter 100</Documentation>
        </IntVectorProperty>""")
    def GetFileID(self, fileID):
        self._fileID = fileID
        return

# Get leading zeros in data file names
    @smproperty.xml("""
        <IntVectorProperty name="Leading zeroes in file names"
            command="GetNumberOfZeros"
            number_of_elements="1"
            default_values="7">
            <IntRangeDomain name="range" />
            <Documentation>If set to 4, the files are umean.dat0100,vmean.dat0100, etc., for File ID 100 </Documentation>
        </IntVectorProperty>""")
    def GetNumberOfZeros(self, nFillZeros):
        self._nFillZeros = nFillZeros
        return

# Get the precision of data files
    @smproperty.xml("""
        <IntVectorProperty name="Data in double precision?"
                           command="GetPrecision"
                           number_of_elements="1"
                           default_values="0">
          <BooleanDomain name="bool"/>
          <Documentation>If ticked, the data will be read in as double. Otherwise, single.</Documentation>
        </IntVectorProperty>""")
    def GetPrecision(self, precision):
        self._precision = precision
        return

# Get the flag for reading scalar data
    @smproperty.xml("""
        <IntVectorProperty name="Read scalars if available? "
                           command="GetScalars"
                           number_of_elements="1"
                           default_values="1">
          <BooleanDomain name="bool"/>
          <Documentation>When checked, scalars will be read if available</Documentation>
        </IntVectorProperty>""")
    def GetScalars(self, scalarFlag):
        self._scalarFlag = scalarFlag
        return

#    @smproperty.xml("""
#        <StringVectorProperty name="Path to data"
#                            command="SetPathToData"
#                            number_of_elements="1"
#                            default_values="./data/">
#            <Documentation>Set the path to the data files</Documentation>
#        </StringVectorProperty>""")
#    def SetPathToData(self, pathToData):
#        self._pathToData = pathToData
#        return

    # Provides meta data information to other pipeline members
    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        
        # Set the extent of the grid - as this reader creates a structured grid
        self._get_grid_size()
        outInfo.Set(executive.WHOLE_EXTENT(), self._domExtent, 6)

        return 1

    # Provides the actual data information to other pipeline members
    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid
        from vtkmodules.numpy_interface import dataset_adapter as dsa

        # Get the data
        raw_data = self._get_raw_data()
        output = dsa.WrapDataObject(vtkRectilinearGrid.GetData(outInfoVec, 0))
        output.ShallowCopy(raw_data)

        # Update the extent of the grid - not sure this is needed
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Get(executive.UPDATE_EXTENT(), raw_data.GetExtent());

        return 1


