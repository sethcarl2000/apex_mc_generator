#ifndef TFileWriter_HH
#define TFileWriter_HH

#include "TrackData_t.hh"

#include <TFile.h>
#include <TTree.h>

#include <mutex>
#include <memory>
#include <vector>


// meyers singleton which can be accessed by any class, and can write data to a TFile.
class TFileWriter {
private:

  TFileWriter() {};
  ~TFileWriter(); 

  std::vector<TrackData_t> fData_R, fData_L; 
  
public:

  //delete copy constructor and assignment operator
  TFileWriter(const TFileWriter&) = delete;
  TFileWriter& operator=(const TFileWriter&) = delete;   

  //file operations
  void OpenFile();
  void CloseFile(); 

  //writes tracks to the current TFile  
  void WriteTrack(); 
  
  //this is how any class or method can get acess to the single, static instance
  static TFileWriter& Instance()
  {
    static TFileWriter inst;
    return inst; 
  }
};

#endif
