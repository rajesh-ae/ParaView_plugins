#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 19:22:23 2021

@author: Rajesh Venkatesan
"""

from paraview.util.vtkAlgorithm import *
from vtk.util import numpy_support
from vtk import vtkStreamingDemandDrivenPipeline
import numpy as np
import os
import re
#------------------------------------------------------------------------------
# Incompact3D Reader
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

@smproxy.reader(name="Incompact3dReader", label="Incompact3d Simulation Data Reader",
                extensions="i3d", file_description="Incompact3d files")
class Incompact3dReader(VTKPythonAlgorithmBase):
    """ A reader that reads a Incompact3D simulation configuration file (.i3d) 
        along with solution data.

        Author: Rajesh Venkatesan
        e-mail: vrajesh.ae[at]gmail.com
    """
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkMultiBlockDataSet')
        self._filename = None
        self._ndata = None
        self._firstFileID = None
        self._nFiles = None
        self._nFillZeros = None
        self._domExtent = None

    def _get_raw_data(self):
        if self._ndata is not None:
            return self._ndata

        if self._filename is None:
            # Include more exceptions like this!
            raise RuntimeError("No filename specified")

        pathToSimulation = os.path.dirname(os.path.abspath(self._filename))

#        print(pathToSimulation)
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
            if 'nclx' in line:
                mynums = re.findall(r"\d+", line)
                nclx = int(mynums[0])
            if 'ncly' in line:
                mynums = re.findall(r"\d+", line)
                ncly = int(mynums[0])
            if 'nclz' in line:
                mynums = re.findall(r"\d+", line)
                nclz = int(mynums[0])

#        print (nx,ny,nz,xlx,yly,zlz,istret,nclx,ncly,nclz)
        i3d_config.close()

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

        from vtkmodules.vtkCommonDataModel import vtkRectilinearGrid,vtkMultiBlockDataSet
        self._ndata = vtkMultiBlockDataSet()
        self._ndata.SetNumberOfBlocks(self._nFiles);

        iblock = 0
        for ifile in range(self._firstFileID,(self._firstFileID + self._nFiles)):
            this_block = vtkRectilinearGrid()
    #        print(dir(this_block))
            this_block.SetDimensions(nx,ny,nz);
            this_block.SetXCoordinates(numpy_support.numpy_to_vtk(x))
            this_block.SetYCoordinates(numpy_support.numpy_to_vtk(y))
            this_block.SetZCoordinates(numpy_support.numpy_to_vtk(z))

    #        print(this_block.GetExtent())
    #        help(numpy_support.numpy_to_vtk)

            nsize = nx*ny*nz

            ux = np.fromfile(pathToSimulation+"/data/ux"+str(ifile).zfill(self._nFillZeros),dtype=np.float32,count=nsize)
            uy = np.fromfile(pathToSimulation+"/data/uy"+str(ifile).zfill(self._nFillZeros),dtype=np.float32,count=nsize)
            uz = np.fromfile(pathToSimulation+"/data/uz"+str(ifile).zfill(self._nFillZeros),dtype=np.float32,count=nsize)
            pressure = np.fromfile(pathToSimulation+"/data/pp"+str(ifile).zfill(self._nFillZeros),dtype=np.float32,count=nsize)
            point_data = this_block.GetPointData()
            vtk_array = numpy_support.numpy_to_vtk(ux)
            vtk_array.SetName('ux')
            point_data.AddArray(vtk_array)

            vtk_array = numpy_support.numpy_to_vtk(uy)
            vtk_array.SetName('uy')
            point_data.AddArray(vtk_array)

            vtk_array = numpy_support.numpy_to_vtk(uz)
            vtk_array.SetName('uz')
            point_data.AddArray(vtk_array)

            vtk_array = numpy_support.numpy_to_vtk(pressure)
            vtk_array.SetName('p')
            point_data.AddArray(vtk_array)

            self._ndata.SetBlock(iblock,this_block);
            iblock = iblock + 1

        return self._get_raw_data()   

    def _get_grid_size(self):
        if self._filename is None:
            # Include more exceptions like this!
            raise RuntimeError("No filename specified")

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
        self._domExtent = (0,nx-1,0,ny-1,0,nz-1)

        return

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="i3d", file_description="Incompact3d config files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self._filename != name:
            self._filename = name
            self._ndata = None
            self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Leading zeroes in file names"
            number_of_elements="1"
            default_values="5"
            command="number_of_zeros">
            <IntRangeDomain name="range" />
            <Documentation>If set to 5, the files are ux00001,uy00002, etc. </Documentation>
        </IntVectorProperty>""")
    def number_of_zeros(self, nFillZeros):
        self._nFillZeros = nFillZeros
        return

    @smproperty.xml("""
        <IntVectorProperty name="File ID"
            number_of_elements="1"
            default_values="1"
            command="get_first_fileID">
            <IntRangeDomain name="range" />
            <Documentation>First file index to read</Documentation>
        </IntVectorProperty>""")
    def get_first_fileID(self, firstFileID):
        self._firstFileID = firstFileID
        return

    @smproperty.xml("""
        <IntVectorProperty name="Number of files"
            number_of_elements="1"
            default_values="1"
            command="get_number_of_files">
            <IntRangeDomain name="range" />
            <Documentation>Number of files to read</Documentation>
        </IntVectorProperty>""")
    def get_number_of_files(self, nFiles):
        self._nFiles = nFiles
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


#        vtkDataObject::SetPointDataActiveScalarInfo(outInfo, this->DataScalarType,
#            this->NumberOfScalarComponents);

    def RequestInformation(self, request, inInfoVec, outInfoVec):
        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        
        # Set the extent of the grid - as this reader creates a structured grid
        self._get_grid_size()
        outInfo.Set(vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(), self._domExtent, 6)

        return 1

    def RequestData(self, request, inInfoVec, outInfoVec):
        from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet
        from vtkmodules.numpy_interface import dataset_adapter as dsa

        raw_data = self._get_raw_data()
        output = dsa.WrapDataObject(vtkMultiBlockDataSet.GetData(outInfoVec, 0))
        output.ShallowCopy(raw_data)
        print(raw_data.GetNumberOfBlocks())
        # Update the extent of the grid - not sure this is needed
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Get(vtkStreamingDemandDrivenPipeline.UPDATE_EXTENT(), raw_data.GetBlock(0).GetExtent());

        return 1


