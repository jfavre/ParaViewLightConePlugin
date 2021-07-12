/*=========================================================================

  Program:   ParaView
  Module:    vtkLightConeReader.cxx

=========================================================================*/

#include "vtkLightConeReader.h"

#include "vtkCellType.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDataArraySelection.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#ifdef ALL_TYPES
#include "vtkMultiBlockDataSet.h"
#endif
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vtksys/RegularExpression.hxx>
#include <vtksys/SystemTools.hxx>

#ifdef PARAVIEW_USE_MPI
#include "vtkMultiProcessController.h"
vtkCxxSetObjectMacro(vtkLightConeReader, Controller, vtkMultiProcessController);
#endif

#include <algorithm>
#include <functional>
#include <map>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkLightConeReader);
//----------------------------------------------------------------------------
vtkLightConeReader::vtkLightConeReader()
{
  this->SetNumberOfInputPorts(0);
  this->DistributedSnapshot      = true;
  this->FileName                 = nullptr;
  this->fp                       = nullptr;
  this->TimeStep                 = 0;
  this->ActualTimeStep           = 0;
  this->CellType                 = 0;
  this->UpdatePiece              = 0;
  this->UpdateNumPieces          = 0;
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  for (int i=0; i< 6; i++)
    {
    this->PartTypes[i] = false;
    this->NumPart_Total[i] = 0;
    }

  this->FieldArrays = {
    {"PartType0", {"Foo", "Bar"} },
    {"PartType1", {"velocity", "id"} }, 
#ifdef ALL_TYPES
    {"PartType1", {"velocity"} },
    {"PartType4", {"Density"} }
#endif
  };
#ifdef PARAVIEW_USE_MPI
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
#endif
}
//----------------------------------------------------------------------------
vtkLightConeReader::~vtkLightConeReader()
{
  this->CloseFile();

  delete [] this->FileName;
  this->PointDataArraySelection->Delete();
  this->PointDataArraySelection = nullptr;

#ifdef PARAVIEW_USE_MPI
  this->SetController(NULL);
#endif
  std::cerr << "destructor..." << std::endl;
}
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
void vtkLightConeReader::CloseFile()
{
  if (this->fp)
    {
    fclose(this->fp);
    std::cerr << "closing " << this->FileName << std::endl;
    }
  this->fp = nullptr;
}
//----------------------------------------------------------------------------
// returns 1 for success, 0 otherwise
int vtkLightConeReader::OpenFile()
{
  if (!this->FileName)
    {
    vtkErrorMacro(<<"FileName must be specified.");
    return 0;
    }
  if(!this->fp)
    {
    this->fp = fopen(this->FileName, "r");
    fseek(this->fp, 0L, SEEK_END);
    long TotalSize = ftell(fp);
    fseek(this->fp, 0L, SEEK_SET);
    
    int beg_marker, end_marker;
 
    size_t r;
    r = fread(&beg_marker,           sizeof(int),          1, this->fp);
    r = fread(&this->Npart,          sizeof(unsigned int), 6, this->fp);
    r = fread(&this->Massarr,        sizeof(double),       6, this->fp);
    r = fread(&this->Time,           sizeof(double),       1, this->fp);
    r = fread(&this->Redshift,       sizeof(double),       1, this->fp);
    r = fread(&this->FlagSfr,        sizeof(unsigned int), 1, this->fp);
    r = fread(&this->FlagFeddback,   sizeof(unsigned int), 1, this->fp);
    r = fread(&this->Nall,           sizeof(unsigned int), 6,this->fp);
    r = fread(&this->FlagCooling,    sizeof(unsigned int), 1, this->fp);
    r = fread(&this->NumFiles,       sizeof(unsigned int), 1, this->fp);
    r = fread(&this->BoxSize,        sizeof(double),       1, this->fp);
    r = fread(&this->Omega0,         sizeof(double),       1, this->fp);
    r = fread(&this->OmegaLambda,    sizeof(double),       1, this->fp);
    r = fread(&this->HubbleParam,    sizeof(double),       1, this->fp);
    r = fread(&this->FlagAge,        sizeof(unsigned int), 1, this->fp);
    r = fread(&this->FlagMetals,     sizeof(unsigned int), 1, this->fp);
    r = fread(&this->NallHW,         sizeof(unsigned int), 6, this->fp);
    r = fread(&this->Flag_entr_ics,  sizeof(unsigned int), 1, this->fp);

    fseek(this->fp, 256+4, SEEK_SET);
    r = fread(&end_marker, sizeof(int), 1, this->fp);
    if(beg_marker != end_marker)
        return 0;

    this->PrintHeader();
    //std::cout << "\nSize on disk   " << TotalSize << std::endl;
    //std::cout << "Estimated Size " << (4+256+4) + (4+this->Npart[1]*12L+4) + (4+this->Npart[1]*12L+4) + (4+this->Npart[1]*8L+4) << "\n" <<std::endl;
    }

  return 1;
}

//----------------------------------------------------------------------------
void vtkLightConeReader::PrintHeader()
{
  int i;
  std::cout << "header file: "          << this->FileName << ": ";
  std::cout << this->Npart[1] << ", " << this->Nall[1] << std::endl;
  /*
    
  std::cout << "this->Npart = ";
  for(i=0; i<6;i++)
    std::cout << this->Npart[i] << ", ";
  std::cout << std::endl;

  std::cout << "this->Massarr = ";
  for(i=0; i<6;i++)
    std::cout << this->Massarr[i] << ", ";
  std::cout << std::endl;
  
  std::cout << "this->Time = "          << this->Time << std::endl;
  std::cout << "this->Redshift = "      << this->Redshift << std::endl;
  std::cout << "this->FlagSfr = "       << this->FlagSfr << std::endl;
  std::cout << "this->FlagFeddback = "  << this->FlagFeddback << std::endl;

  std::cout << "this->Nall = ";
  for(i=0; i<6;i++)
    std::cout << this->Nall[i]  << ", ";
  std::cout << std::endl;

  std::cout << "this->FlagCooling = "   << this->FlagCooling  << std::endl;
  std::cout << "this->NumFiles = "      << this->NumFiles << std::endl;
  std::cout << "this->BoxSize = "       << this->BoxSize  << std::endl;
  std::cout << "this->Omega0 = "        << this->Omega0  << std::endl;
  std::cout << "this->OmegaLambda = "   << this->OmegaLambda << std::endl;
  std::cout << "this->HubbleParam = "   << this->HubbleParam << std::endl;
  std::cout << "this->FlagAge = "       << this->FlagAge << std::endl;
  std::cout << "this->FlagMetals = "    << this->FlagMetals << std::endl;

  std::cout << "this->NallHW = ";
  for(i=0; i<6; i++)
    std::cout << this->NallHW[i] << ", ";
  std::cout << std::endl;

  std::cout << "this->Flag_entr_ics = " << this->Flag_entr_ics << std::endl;
  */
}

//----------------------------------------------------------------------------
int vtkLightConeReader::RequestInformation(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  outInfo->Set(CAN_HANDLE_PIECE_REQUEST(), 1);

#ifdef PARAVIEW_USE_MPI
  if (this->Controller)
    {
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
    }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
#endif

  if(this->OpenFile())
    {
    for (int i=0; i< 6; i++)
      {
      if(this->Nall[i])
        this->PartTypes[i] = true;
      }
    this->PointDataArraySelection->AddArray("velocity");
    this->PointDataArraySelection->AddArray("id");
    return 1;
    }
  else
    return 0;
}

void vtkLightConeReader::ReadINT64Dataset(const char *name, vtkTypeInt64* p, long offset=0, long size=0)
{
  size_t r;
  int beg_marker, end_marker;
  // skip over header, coordinates, and velocity
  offset += 4+256+4 + 8*(1+1) + this->Npart[1] * sizeof(float) * (3 + 3);

  fseek(this->fp, offset, SEEK_SET);
  r = fread(&beg_marker, sizeof(int), 1, this->fp);
  if(!size) // by default size is 0 and we read ALL
    size = this->Npart[1];
  r = fread(p, sizeof(vtkTypeInt64), size, this->fp);
  if(r != size)
    std::cerr << "error reading bulk data array for" << name << std::endl;

  r = fread(&end_marker, sizeof(int), 1, this->fp);
  if(beg_marker != end_marker)
    std::cerr << "ReadINT64Dataset(): error reading begin(" << beg_marker <<  ") and end markers(" << end_marker << ") for " << name << std::endl;
};

void vtkLightConeReader::ReadFloatDataset(const char *name, float* p, long offset=0, size_t size=0)
{
// offset will be given as 3*number of particles in subset
  size_t r;
  //int beg_marker, end_marker;
  if(!strcmp(name, "Coordinates"))
    offset += 4+256+4;
  else // velocity
    offset += 4+256+4 + 8 + this->Npart[1]*sizeof(float) * 3;

  offset += sizeof(int); // to account for begin-of-block marker

  fseek(this->fp, offset, SEEK_SET);
  //r = fread(&beg_marker, sizeof(int), 1, this->fp);
  if(!size) // by default size is 0 and we read ALL
    size = this->Npart[1] * 3;

  std::cerr << "offset = " << offset <<  " and size = " << size << std::endl;
  //size_t ChunkSize = 268435456; //2^28;
  size_t ChunkSize = 536870912; //2^29;
  long Loops = size / ChunkSize;
  // repeat "Loops" times and then finish up with smaller chunk
  for(long i=0 ; i < Loops; i++){
    std::cerr << "offset = " << i*ChunkSize <<  " and size = " << ChunkSize << std::endl;
    r = fread(&p[i*ChunkSize], sizeof(float), ChunkSize, this->fp);
    if(r != ChunkSize)
      std::cerr << "error reading bulk data array for" << name << std::endl;
  }
  std::cerr << "offset = " << Loops*ChunkSize <<  " and size = " << size - Loops*ChunkSize << std::endl;
  r = fread(&p[Loops*ChunkSize], sizeof(float), size - Loops*ChunkSize, this->fp);
  if(r != (size - Loops*ChunkSize))
    std::cerr << "error reading bulk data array for" << name << std::endl;
}

unsigned int split_particlesSet(unsigned int N, int piece, int numPieces, long &offset)
{
// I call standard_load, the size of all the first pieces except the last one
// I call load, the size of the last piece
  unsigned int load;
  if(numPieces == 1)
    {
    load = N;
    offset=0;
    }
  else
    {
    load = N / numPieces;
    if (piece == (numPieces-1))
      {
      load = N - (numPieces-1) * load;
      }
    }
  offset = (N / numPieces) * piece;
  return load;
}

int vtkLightConeReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkDataObject* doOutput = outInfo->Get(vtkDataObject::DATA_OBJECT());
#ifdef ALL_TYPES
  vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::SafeDownCast(doOutput);
  if (!mb)
    {
    return 0;
    }
#else
#ifdef OUTPUT_UG
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(doOutput);
#else
  vtkPolyData* output = vtkPolyData::SafeDownCast(doOutput);
#endif
#endif

#ifdef PARAVIEW_USE_MPI
  if (this->Controller &&
      (this->UpdatePiece != this->Controller->GetLocalProcessId() ||
       this->UpdateNumPieces != this->Controller->GetNumberOfProcesses()))
  {
    vtkDebugMacro(<< "Parallel failure, Id's not right (ignore)");
    this->UpdatePiece = this->Controller->GetLocalProcessId();
    this->UpdateNumPieces = this->Controller->GetNumberOfProcesses();
  }
#else
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
#endif

#define PARALLEL_DEBUG 1
#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/scratch/snx3000/jfavre/out." << this->UpdatePiece << ".txt" << ends;
  std::ofstream errs;
  errs.open(fname.str().c_str(), ios::app);
/*
  errs << "piece " << this->UpdatePiece << " out of " << this->UpdateNumPieces << endl;
  errs << "int of size " <<  sizeof(int) << endl;
  errs << "size_t of size " <<  sizeof(size_t) << endl;
  errs << "long of size " <<  sizeof(long) << endl;
  errs << "vtkIdType of size " <<  sizeof(vtkIdType) << endl;
*/
#endif

// must correct for overflow.
  long temp;
  for(int i=0; i < 6; i++)
    if(this->Nall[i] < 0)
      {
      temp = this->Nall[i];
      temp += 4294967296;
#ifdef PARALLEL_DEBUG
  errs << "OVERFLOW detected temp= " << temp << ", this->Nall[i]= " << this->Nall[i] << endl;
#endif
      this->Nall[i] = temp;
      }

  long LoadPart_Total[6] ={0,0,0,0,0,0};
  long ParallelOffset[6];

  for (int i=0; i< 6; i++)
      {
      LoadPart_Total[i] = split_particlesSet(this->Npart[i],  this->UpdatePiece,  this->UpdateNumPieces, ParallelOffset[i]);
      errs << "LoadPart_Total["<< i << "] = " << LoadPart_Total[i] << ", Nall["<< i << "] = " << this->Npart[i]<< endl;
      }


#ifdef ALL_TYPES
#ifdef OUTPUT_UG
  vtkUnstructuredGrid *output;
#else
  vtkPolyData *output;
#endif
#endif

  vtkFloatArray  *data;
  vtkIdTypeArray *uidata;

  int validPart, myType;

  for(validPart=0, myType = Gas; myType<= Stars; myType++)
    {
    if(PartTypes[myType])
      {
#ifdef PARALLEL_DEBUG
      errs << " creating PolyData for PartType " << myType << " with " << LoadPart_Total[myType] << " points at offset " << ParallelOffset[myType] << std::endl;
#endif

#ifdef ALL_TYPES
#ifdef OUTPUT_UG
      output = vtkUnstructuredGrid::New();
#else
      output = vtkPolyData::New();
#endif

      mb->SetBlock(validPart, output);
      mb->GetMetaData(validPart)->Set(vtkCompositeDataSet::NAME(), ParticleTypes[myType]);
      output->Delete();
#endif
      vtkDoubleArray *cst = vtkDoubleArray::New();
      cst->SetNumberOfComponents(1);
      cst->SetNumberOfTuples(6);
      cst->SetName("Massarr");
      for(auto i=0; i < 6; i++)
        cst->SetTuple1(i, this->Massarr[i]);
      output->GetFieldData()->AddArray(cst);
      cst->Delete();

      vtkDoubleArray *cst0 = vtkDoubleArray::New();
      cst0->SetNumberOfComponents(1);
      cst0->SetNumberOfTuples(1);
      cst0->SetName("Time");
      cst0->SetTuple1(0, this->Time);
      output->GetFieldData()->AddArray(cst0);
      cst0->Delete();
      
      vtkDoubleArray *cst1 = vtkDoubleArray::New();
      cst1->SetNumberOfComponents(1);
      cst1->SetNumberOfTuples(5);
      cst1->SetName("constants");
      cst1->SetTuple1(0, this->Redshift);
      cst1->SetTuple1(1, this->BoxSize);
      cst1->SetTuple1(2, this->Omega0);
      cst1->SetTuple1(3, this->OmegaLambda);
      cst1->SetTuple1(4, this->HubbleParam);
      output->GetFieldData()->AddArray(cst1);
      cst1->Delete();
  
      vtkFloatArray *coords = vtkFloatArray::New();
      coords->SetNumberOfComponents(3);
      coords->SetNumberOfTuples(LoadPart_Total[myType]);
      std::cerr << " creating coordinates array of size " << LoadPart_Total[myType] << " points\n";
      coords->SetName("coords");

      vtkPoints *points = vtkPoints::New();
      points->SetData(coords);
      output->SetPoints(points);
      coords->Delete();
      points->Delete();

      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i)))
          {
          const char *name = &this->GetPointArrayName(i)[0];
#ifdef PARALLEL_DEBUG
      errs << " creating data array " << name << " for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif

          if(!strcmp(name, "velocity"))
            {
            data = vtkFloatArray::New();
            data->SetNumberOfComponents(3);
            data->SetNumberOfTuples(LoadPart_Total[myType]);
      std::cerr << " creating    velocity array of size " << LoadPart_Total[myType] << " points\n";
            data->SetName("velocity");
            output->GetPointData()->AddArray(data);
            data->Delete();
            }
          else if(!strcmp(name, "id"))
            {
            uidata = vtkIdTypeArray::New();
            uidata->SetNumberOfComponents(1);
            uidata->SetNumberOfTuples(LoadPart_Total[myType]);
      std::cerr << " creating id array of size " << LoadPart_Total[myType] << " points\n";
            uidata->SetName("id");
            output->GetPointData()->AddArray(uidata);
            output->GetPointData()->SetGlobalIds(uidata);
            uidata->Delete();
            }
          }
        }

      if (this->CellType == CellTypes::Vertex)
        {
        vtkCellArray *vertices =  vtkCellArray::New();
#ifdef VTK_CELL_ARRAY_V2
#else
        vtkIdType* cells = vertices->WritePointer(LoadPart_Total[myType], 2 * LoadPart_Total[myType]);
        for (vtkIdType i = 0; i < LoadPart_Total[myType]; ++i)
          {
          cells[2 * i] = 1;
          cells[2 * i + 1] = i;
          }
#endif
#ifdef OUTPUT_UG
       output->SetCells(VTK_VERTEX, vertices);
#else
       output->SetVerts(vertices);
#endif
       vertices->Delete();
       }
     else if (this->CellType == CellTypes::PolyVertex)
        {
        vtkIdList *list = vtkIdList::New();
        list->SetNumberOfIds(LoadPart_Total[myType]);
        for(vtkTypeInt64 i=0; i < LoadPart_Total[myType]; i++)
          list->SetId(i, i);
        output->Allocate(1);
        output->InsertNextCell(VTK_POLY_VERTEX, list);
        list->Delete();
        }
      validPart++;
      }
    }

  int offset, fileOffsetNodes[6] = {0,0,0,0,0,0};

    for(validPart=0, myType = Gas; myType<= Stars; myType++)
      {
      char name[16];
      sprintf(name,"PartType%1d", myType);
      if(PartTypes[myType])
        {
#ifdef ALL_TYPES
#ifdef OUTPUT_UG
        output = static_cast<vtkUnstructuredGrid*>(mb->GetBlock(validPart));
#else
        output = static_cast<vtkPolyData*>(mb->GetBlock(validPart));
#endif
#endif
        long offset;
        offset = ParallelOffset[myType]*sizeof(float)*3;
        size_t size = LoadPart_Total[myType] * 3;
// insert coordinates read
        float *dptr = static_cast<vtkFloatArray *>(output->GetPoints()->GetData())->GetPointer(fileOffsetNodes[myType]*3);
        ReadFloatDataset("Coordinates", dptr, offset, size);
// end of coordinates read

// insert PointData here
      for(int i = 0 ; i < this->GetNumberOfPointArrays(); i++)
        {
        if(this->GetPointArrayStatus(this->GetPointArrayName(i)))
          {
            const char *name = &this->GetPointArrayName(i)[0];
#ifdef PARALLEL_DEBUG
      errs << this->FileName << ": reading data array " << name << " for PartType " << myType << " with " << LoadPart_Total[myType] << " points\n";
#endif
          if(!strcmp(name, "id"))
            {
            uidata = static_cast<vtkIdTypeArray *>(output->GetPointData()->GetArray(name));
            vtkTypeInt64 *lptr = uidata->GetPointer( 0 );
            offset = ParallelOffset[myType]*sizeof(vtkTypeInt64);
            size = LoadPart_Total[myType];
            ReadINT64Dataset("id", lptr, offset, size);
            }
          else
            {
            data = static_cast<vtkFloatArray *>(output->GetPointData()->GetArray(name));
            dptr = data->GetPointer( 0 );
            offset = ParallelOffset[myType]*sizeof(float)*3;
            size = LoadPart_Total[myType] * 3;
            ReadFloatDataset(name, dptr, offset, size);
            // split by component
/*
            if(!strcmp(name, "velocity"))
              {
              double tuple[3];
              vtkIdType NbTuples = data->GetNumberOfTuples();
              vtkFloatArray *vx = vtkFloatArray::New();
              vx->SetName("vx");
              vx->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vx);
              vx->Delete();

              vtkFloatArray *vy = vtkFloatArray::New();
              vy->SetName("vy");
              vy->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vy);
              vy->Delete();

              vtkFloatArray *vz = vtkFloatArray::New();
              vz->SetName("vz");
              vz->SetNumberOfTuples(NbTuples);
              output->GetPointData()->AddArray(vz);
              vz->Delete();
              for(vtkIdType i=0; i < NbTuples; i++)
                {
                data->GetTuple(i, tuple);
                vx->SetTuple1(i, tuple[0]);
                vy->SetTuple1(i, tuple[1]);
                vz->SetTuple1(i, tuple[2]);
                }
              }
*/
            }
          }
        }
// end of PointData read

        validPart++;
        }
      fileOffsetNodes[myType] += offset;
      } // for all part types

  this->CloseFile();

#ifdef PARALLEL_DEBUG
  errs.close();
#endif
  return 1;
}


//----------------------------------------------------------------------------
const char* vtkLightConeReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkLightConeReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkLightConeReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
    {
    if (status)
      {
      this->PointDataArraySelection->EnableArray(name);
      }
    else
      {
      this->PointDataArraySelection->DisableArray(name);
      }
    this->Modified();
    }
}

void vtkLightConeReader::EnableAllPointArrays()
{
    this->PointDataArraySelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkLightConeReader::DisableAllPointArrays()
{
    this->PointDataArraySelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
int vtkLightConeReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
void vtkLightConeReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "File Name: " <<
    (this->FileName ? this->FileName : "(none)") << "\n";
}
