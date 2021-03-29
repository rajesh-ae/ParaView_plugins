from vtkmodules.vtkCommonDataModel import vtkDataSet
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs
from vtkmodules.vtkCommonCore import vtkDataArraySelection
import numpy as np

# new module for ParaView-specific decorators.
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain

@smproxy.filter(label="Compressible Flow")
@smproperty.input(name="Input")
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class CompressibleFlowFilter(VTKPythonAlgorithmBase):
    """
    This filter calculates the compressible flow properties for the dataset.
    Input dataset must have: pressure,density,velocity (vector or components)
    Other input: Gas properties - specific heat, specific heat ratio
    Author: Rajesh Venkatesan
    e-mail: vrajesh.ae[at]gmail.com
    """
    def __init__(self):
        super().__init__(nInputPorts=1, nOutputPorts=1, outputType="vtkDataSet")
        self._gamma = 0.0
        self._Cp = 0.0
        self._pressure_array = "whats_in_a_name"
        self._density_array = "whats_in_a_name"
        self._vel_vector = 0
        self._velocity_array = "whats_in_a_name"
        self._velocity_array1 = "whats_in_a_name"
        self._velocity_array2 = "whats_in_a_name"
        self._velocity_array3 = "whats_in_a_name"

    def RequestDataObject(self, request, inInfo, outInfo):
        inData = self.GetInputData(inInfo, 0, 0)
        outData = self.GetOutputData(outInfo, 0)
        assert inData is not None
        if outData is None or (not outData.IsA(inData.GetClassName())):
            outData = inData.NewInstance()
            outInfo.GetInformationObject(0).Set(outData.DATA_OBJECT(), outData)
        return super().RequestDataObject(request, inInfo, outInfo)

    def RequestData(self, request, inInfo, outInfo):
        input0 = dsa.WrapDataObject(vtkDataSet.GetData(inInfo[0],0))

        gamma = self._gamma
        Cp = self._Cp
        Cv = Cp/gamma
        R = Cp - Cv
        pr_name = self._pressure_array
        rho_name = self._density_array
        ga1 = gamma/(gamma - 1)

        if(self._vel_vector):
            vel = input0.PointData[self._velocity_array]
        else:
            vel = algs.make_vector(input0.PointData[self._velocity_array1],input0.PointData[self._velocity_array2],input0.PointData[self._velocity_array3])
        Ts = input0.PointData[pr_name] /(R*input0.PointData[rho_name])
        a = np.sqrt(gamma*R*Ts)
        vmag = np.linalg.norm(vel,axis=1)
        Mach = vmag/a
        T0 = Ts + (0.5/Cp)*(vmag*vmag)
        P0 = input0.PointData[pr_name]*((T0/Ts)**ga1)

        output = dsa.WrapDataObject(vtkDataSet.GetData(outInfo))
        output.ShallowCopy(input0.VTKObject)
        output.PointData.append(T0, "T0");
        output.PointData.append(P0, "P0");
        if(self._vel_vector == 0):
            output.PointData.append(vel, "vel");
        output.PointData.append(Ts, "Ts");
        output.PointData.append(Mach, "M");

        return 1

    @smproperty.xml("""
        <DoubleVectorProperty name="Specific heat, Cp"
            number_of_elements="1"
            default_values="1006.0"
            command="get_Cp">
            <DoubleRangeDomain name="range" />
            <Documentation>Set the specific heat at constant pressure (J/kg/K)</Documentation>
        </DoubleVectorProperty>""")
    def get_Cp(self, Cp):
        self._Cp = Cp
        return

    @smproperty.xml("""
        <DoubleVectorProperty name="Specific heat ratio"
            number_of_elements="1"
            default_values="1.4"
            command="get_gamma">
            <DoubleRangeDomain name="range" />
            <Documentation>Set the specific heat ratio</Documentation>
        </DoubleVectorProperty>""")
    def get_gamma(self, gamma):
        self._gamma = gamma
        return

    @smproperty.xml("""
        <StringVectorProperty name="SelectPressure"
                            label="Pressure"
                            command="SelectPressureArray"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")
    def SelectPressureArray(self, idx, port, connection, field, name):
        self._pressure_array = name
        return

    @smproperty.xml("""
        <StringVectorProperty name="SelectDensity"
                            label="Density"
                            command="SelectDensityArray"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")

    def SelectDensityArray(self, idx, port, connection, field, name):
        self._density_array = name
        return

    @smproperty.xml("""
  <IntVectorProperty name="Velocity vector available? "
                     command="check_vel_vector"
                     number_of_elements="1"
                     default_values="1">
    <BooleanDomain name="bool"/>
  </IntVectorProperty>""")
    def check_vel_vector(self, x):
        self._vel_vector = x
        return

    @smproperty.xml("""
        <StringVectorProperty name="SelectVelocity"
                            label="Velocity"
                            command="SelectVelocityArray"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Vectors"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")
    def SelectVelocityArray(self, idx, port, connection, field, name):
        self._velocity_array = name
        return

    @smproperty.xml("""
        <StringVectorProperty name="SelectVelocity1"
                            label="VelocityX"
                            command="SelectVelocityComp1"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")
    def SelectVelocityComp1(self, idx, port, connection, field, name):
        self._velocity_array1 = name
        return

    @smproperty.xml("""
        <StringVectorProperty name="SelectVelocity2"
                            label="VelocityY"
                            command="SelectVelocityComp2"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")
    def SelectVelocityComp2(self, idx, port, connection, field, name):
        self._velocity_array2 = name
        return
    @smproperty.xml("""
        <StringVectorProperty name="SelectVelocity3"
                            label="VelocityZ"
                            command="SelectVelocityComp3"
                            number_of_elements="5"
                            element_types="0 0 0 0 2"
                            animateable="0">
        <ArrayListDomain name="array_list"
                         attribute_type="Scalars"
                         input_domain_name="inputs_array">
          <RequiredProperties>
            <Property name="Input"
                      function="Input" />
          </RequiredProperties>
        </ArrayListDomain>
      </StringVectorProperty>""")
    def SelectVelocityComp3(self, idx, port, connection, field, name):
        self._velocity_array3 = name
        return
