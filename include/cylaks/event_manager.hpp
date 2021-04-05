#ifndef _CYLAKS_EVENT_MANAGER_HPP_
#define _CYLAKS_EVENT_MANAGER_HPP_
#include "event.hpp"
#include "system_namespace.hpp"

// EventManager: Handles all Events(); ensures no populations go negative
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
      Sys::Log("p_%s = %g [%zu exe]\n", event.name_.c_str(),
               double(event.n_executed_tot_) / event.n_opportunities_tot_,
               event.n_executed_tot_);
    }
  }
  void Initialize();
  void ExecuteEvents();
};

#endif