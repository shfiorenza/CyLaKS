#ifndef _CYLAKS_EVENT_MANAGER_HPP_
#define _CYLAKS_EVENT_MANAGER_HPP_
#include "event.hpp"

class Curator;
struct SysRNG;

class EventManager {
private:
  int n_events_to_exe_;
  Vec<Event *> events_to_exe_;

  Curator *wally_;
  SysRNG *gsl_;

public:
  Vec<Event> events_;

private:
  void SampleEventStatistics();
  void GenerateExecutionSequence();

public:
  EventManager();
  void Initialize(Curator *wally);
  void ExecuteEvents();
};

#endif