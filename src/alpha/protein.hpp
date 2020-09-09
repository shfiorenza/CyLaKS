#ifndef _PROTEIN
#define _PROTEIN
#include "base_object.hpp"
class Tubulin;

class Protein : public BaseObject {
private:
  // Index of this protein in its manager's list
  unsigned int index_;
  // Whether or not protein is active in simulation (i.e., bound to something)
  bool active_;
  // Points to occupied site when bound; otherwise nullptr
  Tubulin *site_;

  unsigned int n_neighbor_sites_;
  Vec<Tubulin *> neighbor_sites_;

public:
  Protein(int id, int index) : BaseObject(id) {
    index_ = index;
    active_ = false;
    site_ = nullptr;
    n_neighbor_sites_ = 0;
  };
  ~Protein() {}

  void InitializeNeighborList(int size) { neighbor_sites_.reserve(size); }

  int GetIndex() { return index_; }

  void Activate() { active_ = true; }
  void Deactivate() { active_ = false; }

  Tubulin *GetSite() { return site_; }
  void SetSite(Tubulin *new_site) { site_ = new_site; }
};

#endif