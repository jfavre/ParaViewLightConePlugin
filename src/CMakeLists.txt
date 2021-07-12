cmake_minimum_required(VERSION 3.12...3.16 FATAL_ERROR)
paraview_add_plugin(pvLightConeReader
  REQUIRED_ON_SERVER
  REQUIRED_ON_CLIENT
  VERSION "1.0"
  MODULES LightConeReader
  MODULE_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Reader/vtk.module")

#if (BUILD_TESTING)
#  add_subdirectory(Testing)
#endif ()

