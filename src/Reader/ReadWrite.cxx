#include <iostream>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
#include "vtkLightConeReader.h"
#include "vtkPointData.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

int main(int argc, char **argv)
{
  VTK_CREATE(vtkLightConeReader, reader);

  reader->SetFileName(argv[1]);
  reader->SetDirectoryName(argv[1]);

  reader->SetPointArrayStatus("id", 1);
  reader->SetPointArrayStatus("velocity", 1);
  reader->Update();

  std::cout << "opening file in write mode" << argv[2] << std::endl;
  FILE *fp = fopen(argv[2],"wb");

  int beg_marker=256, end_marker=256;
  double Mass[6]={0,0,0,0,0,0}, lambda=0.;
  unsigned int NpartTotal[6]={0,0,0,0,0,0}, luint, NumFiles=4, Npart[6]={0,0,0,0,0,0};
  unsigned int NpartTotalHW[6]={0,0,0,0,0,0};
  
  Npart[1] = reader->GetOutput()->GetNumberOfPoints();
  // using 4 files and must execute 4 time to create 4 distinct files
  NumFiles = 4;
  NpartTotal[1] = NumFiles * Npart[1];
  
  fwrite(&beg_marker,     sizeof(int),          1, fp);
  fwrite(Npart,          sizeof(unsigned int), 6, fp);
  fwrite(Mass,            sizeof(double),       6, fp);
  fwrite(&lambda,           sizeof(double),       1, fp);
  fwrite(&lambda,       sizeof(double),       1, fp);
  fwrite(&luint,        sizeof(unsigned int), 1, fp);
  fwrite(&luint,   sizeof(unsigned int), 1, fp);
  fwrite(NpartTotal,      sizeof(unsigned int), 6,fp);
  fwrite(&luint,    sizeof(unsigned int), 1, fp);
  fwrite(&NumFiles,       sizeof(unsigned int), 1, fp);
  fwrite(&lambda,        sizeof(double),       1, fp);
  fwrite(&lambda,         sizeof(double),       1, fp);
  fwrite(&lambda,    sizeof(double),       1, fp);
  fwrite(&lambda,    sizeof(double),       1, fp);
  fwrite(&luint,        sizeof(unsigned int), 1, fp);
  fwrite(&luint,     sizeof(unsigned int), 1, fp);
  fwrite(NpartTotalHW,    sizeof(unsigned int), 6, fp);
  fwrite(&luint,  sizeof(unsigned int), 1, fp);
  //fill up with 15 4-byte numbers to fil up
  for(auto i=0; i < 15; i++)
    fwrite(&beg_marker,     sizeof(int),          1, fp);
  fwrite(&end_marker,       sizeof(int),          1, fp);
  //////////////////////////
  // now coordinates;
  end_marker = beg_marker = 3*Npart[1]*sizeof(float);
  fwrite(&beg_marker,     sizeof(int),          1, fp);
  float *coords = (float*)reader->GetOutput()->GetPoints()->GetVoidPointer(0);
  for(auto i=0; i < Npart[1]; i++)
    coords[3*i] = coords[3*i] + 0.0;
    //coords[3*i] += 50000;
    //coords[3*i+1] += + 50000;
    //{coords[3*i] += 50000; coords[3*i+1] += 50000;}
  fwrite(coords, sizeof(float)*3, Npart[1], fp);
  fwrite(&end_marker,     sizeof(int),          1, fp);
    //////////////////////////
    // now velocity;
  end_marker = beg_marker = 3*Npart[1]*sizeof(float);
  fwrite(&beg_marker,     sizeof(int),          1, fp);
  float *vel = (float*)reader->GetOutput()->GetPointData()->GetArray("velocity")->GetVoidPointer(0);
  fwrite(vel, sizeof(float)*3, Npart[1], fp);
  fwrite(&end_marker,     sizeof(int),          1, fp);
    //////////////////////////
    // now 64-bit id;
  end_marker = beg_marker = Npart[1]*sizeof(long);
  fwrite(&beg_marker,     sizeof(int),          1, fp);
  long *id = (long*)reader->GetOutput()->GetPointData()->GetArray("id")->GetVoidPointer(0);
  fwrite(id, sizeof(long), Npart[1], fp);
  fwrite(&end_marker,     sizeof(int),          1, fp);
  fclose(fp);
  exit(1);
}
