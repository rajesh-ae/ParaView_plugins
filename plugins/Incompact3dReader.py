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
# Incompact3D Instantaneous Data Reader
#------------------------------------------------------------------------------
@smproxy.reader(name="Incompact3dReader", label="Incompact3d Simulation Data Reader",
                extensions="i3d", file_description="Incompact3d files")
class Incompact3dReader(VTKPythonAlgorithmBase):
    """ A reader that reads a Incompact3D simulation configuration file (.i3d) 
        along with instantaneous solution data (ux##, uy##, uz##, etc.).

        Author: Rajesh Venkatesan
        e-mail: vrajesh.ae[at]gmail.com
    """
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkRectilinearGrid')
        self._filename = None
        self._ndata = None
        self._firstFileID = None
        self._nFiles = None
        self._nFillZeros = None
        self._domExtent = None
        self._precision = 0
        self._collection = None
        self._timesteps = None
        self._scalarFlag = 1

    # Get the data - Called during RequestData phase
    def _get_raw_data(self, requested_time=None):
        if self._collection is not None:
            if requested_time is not None:
                self._ndata = self._collection.GetItem(int(requested_time))
                return self._ndata
            self._ndata = self._collection.GetItem(0)
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
                y[j] = j*dy

        for k in range(nz):
            z[k] = k*dz

        # Total grid size
        nsize = nx*ny*nz

        # Set the precision to read in data files
        if (self._precision):
            file_precision = np.dtype(np.float64)
        else:
            file_precision = np.dtype(np.float32) 

        # Create a list of all available instantaneous variables
        list_of_arrays = ['ux','uy','uz','pp']

        # Add scalars to the list if needed
        if (nscalar) and (self._scalarFlag):
            for iscalar in range(nscalar):
                this_scalar = 'phi'+str(iscalar+1)
                list_of_arrays.append(this_scalar)
            

        # The data will be stored as vtkDataSetCollection - with each instant being a vtkRectilinearGrid
        from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid,vtkDataSetCollection
        self._collection = vtkDataSetCollection()

        # Read arrays for time instants specified
        for ifile in range(self._firstFileID,(self._firstFileID + self._nFiles)):
            this_block = vtkRectilinearGrid()

            # Set the grid parameters - see if these can be mapped instead of setting for each instant
            this_block.SetDimensions(nx,ny,nz);
            this_block.SetXCoordinates(numpy_support.numpy_to_vtk(x))
            this_block.SetYCoordinates(numpy_support.numpy_to_vtk(y))
            this_block.SetZCoordinates(numpy_support.numpy_to_vtk(z))

            point_data = this_block.GetPointData()
            # Try to read all instantaneous variables
            for this_array in list_of_arrays:
                try:
                    curr_array = np.fromfile(pathToSimulation+"/data/"+this_array+str(ifile).zfill(self._nFillZeros),dtype=file_precision,count=nsize)
                    vtk_array = numpy_support.numpy_to_vtk(curr_array)
                    vtk_array.SetName(this_array)
                    point_data.AddArray(vtk_array)          
                except FileNotFoundError:
                    # Skip the unavailable ones - move this to log message
                    print('File',this_array+str(ifile).zfill(self._nFillZeros),'does not exist! Skipping this file ...')

            # Add this instant to the collection
            self._collection.AddItem(this_block)

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

# Get the list of timesteps - called during RequestInformation phase among others
    def _get_timesteps(self):
        self._timesteps = np.arange(self._nFiles,dtype=np.float32)
        return self._timesteps.tolist() if self._timesteps is not None else None

# Update the list of timesteps - called during RequestData phase
    def _get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self._get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

# Get i3d file name 
    @smproperty.stringvector(name="Config.FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="i3d", file_description="Incompact3d config files")
    def GetConfigFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self._collection = None
            self._ndata = None
            self.Modified()

# TimestepValues - only to enable time step feature
    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._get_timesteps()

# Get first data file ID
    @smproperty.xml("""
        <IntVectorProperty name="First file ID"
            command="GetFirstFileID"
            number_of_elements="1"
            default_values="1">
            <IntRangeDomain name="range" />
            <Documentation>First file index to read</Documentation>
        </IntVectorProperty>""")
    def GetFirstFileID(self, firstFileID):
        self._firstFileID = firstFileID
        return

# Get number of data files to read
    @smproperty.xml("""
        <IntVectorProperty name="Number of files"
            number_of_elements="1"
            default_values="1"
            command="GetNumberOfFiles">
            <IntRangeDomain name="range" />
            <Documentation>Number of files to read</Documentation>
        </IntVectorProperty>""")
    def GetNumberOfFiles(self, nFiles):
        if self._nFiles != nFiles:
            self._nFiles = nFiles
            self._collection = None
            self._ndata = None
            self.Modified()


# Get leading zeros in data file names
    @smproperty.xml("""
        <IntVectorProperty name="Leading zeroes in file names"
            command="GetNumberOfZeros"
            number_of_elements="1"
            default_values="5">
            <IntRangeDomain name="range" />
            <Documentation>If set to 5, the files are ux00001,uy00002, etc.</Documentation>
        </IntVectorProperty>""")
    def GetNumberOfZeros(self, nFillZeros):
        self._nFillZeros = nFillZeros
        return

# Get the precision of data files
    @smproperty.xml("""
        <IntVectorProperty name="Data in double precision? "
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

    # Provides meta data information to other pipeline members
    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        
        # Set the extent of the grid - as this reader creates a structured grid
        self._get_grid_size()
        outInfo.Set(executive.WHOLE_EXTENT(), self._domExtent, 6)

        # Set the time information of the data
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())
        timesteps = self._get_timesteps()
        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])

        return 1

    # Provides the actual data information to other pipeline members
    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid
        from vtkmodules.numpy_interface import dataset_adapter as dsa

        # Get the timestamp of the data requested
        data_time = self._get_update_time(outInfoVec.GetInformationObject(0))

        # Get the actual data corresponding to this timestamp - as vtkRectilinearGrid dataset
        raw_data = self._get_raw_data(data_time)
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)

        output = dsa.WrapDataObject(vtkRectilinearGrid.GetData(outInfoVec, 0))
        output.ShallowCopy(raw_data)

        # Update the extent of the grid - not sure this is needed
        outInfo.Get(executive.UPDATE_EXTENT(), raw_data.GetExtent());

        if data_time is not None:
            output.GetInformation().Set(output.DATA_TIME_STEP(), data_time)

        return 1


