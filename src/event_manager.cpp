#include "cylaks/event_manager.hpp"
#include "cylaks/curator.hpp"
#include "cylaks/system_namespace.hpp"
#include "cylaks/system_rng.hpp"

void EventManager::SampleEventStatistics() {

  n_events_to_exe_ = 0;
  for (auto &&event : events_) {
    // Sys::Log(1, "Sampling event %s\n", event.name_.c_str());
    n_events_to_exe_ += event.SampleStatistics();
  }
  if (n_events_to_exe_ <= 1) {
    return;
  }
  Sys::Log(1, "Validating event statistics\n");
  // Put all events with >1 target into an active_ array
  Pair<Event *, Object *> active_events[n_events_to_exe_];
  int i_active{0};
  for (auto &&event : events_) {
    // Add a ptr to the event for each target it has
    for (int i_tar{0}; i_tar < event.n_expected_; i_tar++) {
      active_events[i_active++] = std::make_pair(&event, event.targets_[i_tar]);
    }
  }
  // Scan through all active events to ensure that no 2 target the same
  for (int i_entry{0}; i_entry < n_events_to_exe_; i_entry++) {
    Event *event_i{active_events[i_entry].first};
    Object *tar_i{active_events[i_entry].second};
    Sys::Log(2, "   event_i = %s @ %zu\n", event_i->name_.c_str(),
             tar_i->GetID());
    for (int j_entry{i_entry + 1}; j_entry < n_events_to_exe_; j_entry++) {
      Event *event_j{active_events[j_entry].first};
      Object *tar_j{active_events[j_entry].second};
      Sys::Log(2, "     event_j = %s @ %zu\n", event_j->name_.c_str(),
               tar_j->GetID());
      bool correction_needed{false};
      // If event_i and event_j target the same protein, remove one
      if (tar_i->GetID() == tar_j->GetID()) {
        correction_needed = true;
      }
      // Also check teth partners since they're coupled
      if (tar_i->IsTethered()) {
        if (tar_i->GetTethPartner()->GetID() == tar_j->GetID()) {
          correction_needed = true;
        }
      }
      if (correction_needed) {
        double p_one{event_i->p_occur_};
        double p_two{event_j->p_occur_};
        double ran{SysRNG::GetRanProb()};
        if (ran < p_one / (p_one + p_two)) {
          event_i->RemoveTarget(active_events[i_entry].second);
          active_events[i_entry] = active_events[n_events_to_exe_ - 1];
          i_entry--;
          n_events_to_exe_--;
          Sys::Log(2, "   Removed event '%s' @ %zu\n", event_i->name_.c_str(),
                   tar_i->GetID());
          break;
        } else {
          event_j->RemoveTarget(active_events[j_entry].second);
          active_events[j_entry] = active_events[n_events_to_exe_ - 1];
          j_entry--;
          n_events_to_exe_--;
          Sys::Log(2, "   Removed event '%s' @ %zu\n", event_j->name_.c_str(),
                   tar_j->GetID());
        }
      }
    }
  }
  Sys::Log(1, " - DONE\n");
}

void EventManager::GenerateExecutionSequence() {

  if (n_events_to_exe_ == 0) {
    return;
  }
  Sys::Log(1, "Generating execute sequence\n");
  int i_array{0};
  Event *pre_array[n_events_to_exe_];
  for (auto &&event : events_) {
    for (int i_entry{0}; i_entry < event.n_expected_; i_entry++) {
      pre_array[i_array++] = &event;
    }
  }
  if (i_array != n_events_to_exe_) {
    Sys::ErrorExit("K_MGMT::GenerateExecutionSequence()");
  }
  if (n_events_to_exe_ > 1) {
    SysRNG::Shuffle(pre_array, n_events_to_exe_, sizeof(Event *));
  }
  if (n_events_to_exe_ > events_to_exe_.size()) {
    events_to_exe_.resize(n_events_to_exe_);
  }
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    events_to_exe_[i_event] = pre_array[i_event];
  }
  Sys::Log(1, " - DONE\n");
}

bool EventManager::ExecuteEvents() {

  SampleEventStatistics();
  GenerateExecutionSequence();
  bool event_executed{false};
  for (int i_event{0}; i_event < n_events_to_exe_; i_event++) {
    Event *event{events_to_exe_[i_event]};
    Sys::Log(1, "Executing event %s on protein #%i\n", event->name_.c_str(),
             event->targets_[event->n_expected_ - 1]->GetID());
    bool successful{event->Execute()};
    if (successful) {
      event_executed = true;
    }
    Sys::Log(1, " DONE\n");
  }
  return event_executed;
}