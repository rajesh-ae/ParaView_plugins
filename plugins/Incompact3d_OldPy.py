#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 06:22:23 2021

@author: rajesh
"""

from vtkmodules.vtkCommonDataModel import vtkDataSet
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs
from vtkmodules.vtkCommonCore import vtkDataArraySelection
import numpy as np

# new module for ParaView-specific decorators.
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

@smproxy.filter(label="Incompact3d")
@smproperty.input(name="Input")
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class Incompact3dFilter(VTKPythonAlgorithmBase):
    """
    This filter provides post-processing support for Incompact3d simulation data
    Input dataset must have: pressure (p), ux, uy, uz

    Author: Rajesh Venkatesan
    e-mail: vrajesh.ae[at]gmail.com
    """
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,nInputPorts=1, nOutputPorts=1,outputType='vtkRectilinearGrid')

    def RequestData(self, request, inInfo, outInfo):
        inData = self.GetInputData(inInfo, 0, 0)
#        outData = self.GetOutputData(outInfo, 0)

        input0 = dsa.WrapDataObject(inData)

        vel = algs.make_vector(input0.PointData['ux'],input0.PointData['uy'],input0.PointData['uz'])
        vmag = np.linalg.norm(vel,axis=1)
        ux2 = input0.PointData['ux'] * 2.0
        output = dsa.WrapDataObject(vtkDataSet.GetData(outInfo))
        output.ShallowCopy(input0.VTKObject)
        output.PointData.append(ux2, "VelIntoTw0");
        output.PointData.append(vel, "Velocity");
        output.PointData.append(vmag, "Vmag");

        return 1