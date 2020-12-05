#ifndef _CYLAKS_SYSTEM_FILES_HPP_
#define _CYLAKS_SYSTEM_FILES_HPP_
#include <filesystem>
#include <iostream>
#include <unordered_map>

struct SysFiles {
private:
  struct DataFile {
    FILE *file_ptr_;
    std::string name_{"nope"};
    std::string filename_{"simName_nope.file"};
    DataFile() {}
    DataFile(std::string sim_name, std::string name) : name_{name} {
      filename_ = sim_name + "_" + name_ + ".file";
      file_ptr_ = fopen(filename_.c_str(), "w");
      if (file_ptr_ == nullptr) {
        printf("Error; cannot open '%s'\n", filename_.c_str());
        exit(1);
      }
    }
    ~DataFile() { fclose(file_ptr_); }
    template <typename DATA_T> void WriteData(DATA_T *array, size_t count) {
      size_t n_chars_written{fwrite(array, sizeof(DATA_T), count, file_ptr_)};
      if (n_chars_written < count) {
        printf("Error writing to '%s'\n", filename_.c_str());
        exit(1);
      }
    }
  };

public:
  FILE *log_;
  std::unordered_map<Str, DataFile> data_;

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
  void AddDataFile(std::string sim_name, std::string name) {
    data_.emplace(name, DataFile(sim_name, name));
  }
};
#endif