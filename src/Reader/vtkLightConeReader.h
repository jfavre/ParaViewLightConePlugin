// .NAME vtkLightConeReader
// .SECTION Description
// vtkLightConeReader reads collection of "LightCone" Gadget files archived by Julian Adamek
// 
// .SECTION Thanks
// Jean Favre
// CSCS - Swiss National Supercomputing Centre for creating and contributing
// this class.

#ifndef vtkLightConeReader_h
#define vtkLightConeReader_h

#include <string>
#include <vector>

#define MULTI_BLOCKS 1

#ifdef MULTI_BLOCKS
#include "vtkMultiBlockDataSetAlgorithm.h"
static std::vector<std::string> ParticleTypes = {"PartType0", "PartType1", "PartType2", "PartType3", "PartType4", "PartType5"};
#else
#include "vtkPolyDataAlgorithm.h"
static std::vector<std::string> ParticleTypes = {"PartType1"};
#endif

#include <map>
#include <sstream>

class LightConeFile;
class vtkDataArraySelection;
//class vtkStdString;
class vtkMultiProcessController;

enum  ParticleType : int {Gas=0, Halo=1, Disk=2, Bulge=3, Stars=4, Bndry=5};
enum  CellTypes : int {None=0, Vertex=1, PolyVertex=2};

#ifdef MULTI_BLOCKS
class vtkLightConeReader : public vtkMultiBlockDataSetAlgorithm
#else
class vtkLightConeReader : public vtkPolyDataAlgorithm
#endif
{
public:
  static vtkLightConeReader *New();
#ifdef MULTI_BLOCKS
  vtkTypeMacro(vtkLightConeReader,vtkMultiBlockDataSetAlgorithm);
#else
  vtkTypeMacro(vtkLightConeReader,vtkPolyDataAlgorithm);
#endif
  void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  
  vtkSetMacro(DistributedSnapshot, vtkTypeBool);
  vtkGetMacro(DistributedSnapshot, vtkTypeBool);
  vtkBooleanMacro(DistributedSnapshot, vtkTypeBool);

  void SetDirectoryName(const char* dn);
  vtkGetStringMacro(DirectoryName);
  
  // Description:
  // Set/Get the timestep to be read
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);

  // Description:
  // When set (default no), the reader will generate a vertex cell
  // for each point/particle read. When using the points directly
  // this is unnecessary and time can be saved by omitting cell generation
  // vtkPointSpriteMapper does not require them.
  // When using ParaView, cell generation is recommended, without them
  // many filter operations are unavailable

  vtkGetMacro(CellType, int);
  vtkSetClampMacro(CellType, int, CellTypes::None, CellTypes::PolyVertex);
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAllPointArrays();
  void        EnableAllPointArrays();
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }
  int         CanReadFile(const char* fname);
#ifdef PARAVIEW_USE_MPI

    // Description:
    // Set/Get the controller use in compositing (set to
    // the global controller by default)
    // If not using the default, this must be called before any
    // other methods.
    virtual void SetController(vtkMultiProcessController* controller);
    vtkGetObjectMacro(Controller, vtkMultiProcessController);

#endif

protected:
   vtkLightConeReader();
  ~vtkLightConeReader();
  //
  int   RequestInformation(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  int   OpenFile(const char*filename);
  void  PrintHeader(const char *filename);
  void  CloseFile(const char *filename);
  void ReadINT64Dataset(const char *, vtkTypeInt64*, long, long);
  void ReadFloatDataset(const char *, float*, long, size_t);
  //
  // Internal Variables

  unsigned int Npart[6];     //local to the current file
  double       Mass[6];
  double       Time;
  double       Redshift;
  int          FlagSfr;
  int          FlagFeddback;
  unsigned int NpartTotal[6];
  int          FlagCooling;
  int          NumFiles;
  double       BoxSize;
  double       Omega0;
  double       OmegaLambda;
  double       HubbleParam;
  int          FlagAge;
  int          FlagMetals;
  unsigned int NpartTotalHW[6];
  int          Flag_entr_ics;

  //
  FILE *fp;
  char          *FileName;
  char          *DirectoryName;
  int           TimeStep;
  int           ActualTimeStep;
  int           CellType;
  bool          PartTypes[6];
  long          NumPart_Total[6];
  vtkTypeBool   DistributedSnapshot;
  typedef std::vector<std::string>  stringlist;
  
  stringlist                        LightConeFileNames;
  std::map<std::string, stringlist> FieldArrays;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

  #ifdef PARAVIEW_USE_MPI
  vtkMultiProcessController* Controller;
  #endif

private:
  vtkLightConeReader(const vtkLightConeReader&) = delete;  // Not implemented.
  void operator=(const vtkLightConeReader&) = delete;  // Not implemented.
};

#endif
