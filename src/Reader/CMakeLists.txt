set(classes
  vtkLightConeReader)

set(private_headers
  vtkLightConeReader.h)

vtk_module_add_module(LightConeReader
  CLASSES ${classes}
  PRIVATE_HEADERS ${private_headers})

paraview_add_server_manager_xmls(
  XMLS  LightCone_Server.xml)

ADD_EXECUTABLE(TestReaderSimple TestReaderSimple.cxx)
TARGET_LINK_LIBRARIES(TestReaderSimple LightConeReader)
