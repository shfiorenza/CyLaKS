#ifndef _CYLAKS_FILE_MANAGER_HPP_
#define _CYLAKS_FILE_MANAGER_HPP_
#include "definitions.hpp"

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
    template <typename DATA_T> void WriteData(DATA_T *array, size_t count) {
      int n_chars_written{fwrite(array, sizeof(DATA_T), count, file_ptr_)};
      if (n_chars_written < 0) {
        printf("Error writing to '%s'\n", filename_.c_str());
        exit(1);
      }
    }
  };
  template <typename DATA_T> struct DataType {
    Str key_{"nope"};
    Str value_type_{"null"};
    DATA_T value_;
    DataType(int val) :
  };

public:
  FILE *log_;
  UMap<Str, SysFile> data_;

private:
public:
  void GenerateLog(char *sim_name) {
    char log_filename[256];
    sprintf(log_filename, "%s.log", sim_name);
    // Check to see if sim files already exist
    if (std::filesystem::exists(log_filename)) {
      printf("Simulation log file with this name already exists!\n");
      printf("Do you wish to overwrite these data? y/n\n");
      std::string response;
      size_t n_responses{0};
      bool response_unacceptable{true};
      while (response_unacceptable) {
        std::getline(std::cin, response);
        if (response == "n" or response == "N") {
          printf("Simulation terminated.\n");
          exit(1);
        } else if (response == "y" or response == "Y") {
          printf("Very well. ");
          printf("Overwriting data for sim '%s'\n\n", sim_name);
          response_unacceptable = false;
        } else {
          printf("ayo I said y or n. try again plz\n");
        }
        n_responses++;
        if (n_responses > 5) {
          printf("aight if u gonna be like that,");
          printf("let's just cancel this whole thing\n");
          exit(1);
        }
      }
    }
    log_ = fopen(log_filename, "w");
    if (log_ == nullptr) {
      printf("Error; cannot open log file '%s'\n", log_filename);
    }
  }
  void AddDataFile(char *sim_name, char *name) {
    data_.emplace(name, SysFile(sim_name, name));
  }
};
#endif