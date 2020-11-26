#ifndef _CYLAKS_FILE_UTILITIES_HPP_
#define _CYLAKS_FILE_UTILITIES_HPP_
#include "definitions.hpp"

template <typename DATA_T> struct DataType {
  Str key_{"nope"};
  Str value_type_{"null"};
  DATA_T value_;
  DataType(int val) :
};

class FileManager {
private:
  struct SysFile {
    Str name_{"nope"};
    Str filename_{"simName_nope.file"};
    FILE *file_ptr_;
    SysFile(char *sim_name, char *name) : name_{name} { OpenFile(sim_name); }
    ~SysFile() { fclose(file_ptr_); }
    void OpenFile(char *sim_name) {
      filename_ = *sim_name + "_" + name_ + ".file";
      file_ptr_ = fopen(filename_.c_str(), "w");
      if (file_ptr_ == nullptr) {
        printf("Error; cannot open '%s'\n", filename_.c_str());
        exit(1);
      }
    }
  };
  template <typename DATA_T> void WriteData(DATA_T *array, size_t count) {
    int n_chars_written{fwrite(array, sizeof(DATA_T), count, file_ptr_)};
    if (n_chars_written < 0) {
      printf("Error writing to '%s'\n", filename_.c_str());
      exit(1);
    }
  }

public:
  SysFile log_;
  UMap<Str, SysFile> data_;

private:
public:
  void GenerateLog();
};
#endif