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
  void Initialize();
  void ExecuteEvents();
};

#endif