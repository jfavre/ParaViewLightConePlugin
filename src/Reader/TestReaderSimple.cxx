#include <iostream>
#include <stdio.h>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
#include "vtkLightConeReader.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

int main(int argc, char **argv)
{
  VTK_CREATE(vtkLightConeReader, reader);
  std::cout << "opening file " << argv[1] << std::endl;
  reader->SetFileName(argv[1]);
  reader->SetDirectoryName(argv[1]);
  
  reader->SetPointArrayStatus("id", 1);
  reader->SetPointArrayStatus("velocity", 1);
  reader->Update();

  exit(1);
}
