#include <iostream>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
#include "vtkLightConeReader.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

int main(int argc, char **argv)
{
  VTK_CREATE(vtkLightConeReader, reader);
  reader->SetFileName(argv[1]);
  reader->UpdateInformation();
  //reader->Enable("Id");
  //reader->Enable("Velocity");
  reader->Update();

  exit(1);
}
