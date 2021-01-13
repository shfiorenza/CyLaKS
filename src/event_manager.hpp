#ifndef _CYLAKS_EVENT_MANAGER_HPP_
#define _CYLAKS_EVENT_MANAGER_HPP_
#include "event.hpp"

class EventManager {
private:
  int n_events_to_exe_;
  Vec<Event *> events_to_exe_;

public:
  Vec<Event> events_;

private:
  void SampleEventStatistics();
  void GenerateExecutionSequence();

public:
  EventManager();
  ~EventManager() {
    for (auto &&event : events_) {
      printf("p_%s = %g\n", event.name_.c_str(),
             double(event.n_executed_tot_) / event.n_opportunities_tot_);
    }
  }
  void Initialize();
  void ExecuteEvents();
};

#endif