#ifndef _CYLAKS_OBJECT_HPP_
#define _CYLAKS_OBJECT_HPP_
#include "definitions.hpp"

class Object {
private:
  size_t unique_id_{0};
  size_t species_id_{0};

protected:
  bool visible_{true};
  Vec<double> pos_;
  int n_neighbors_{0};
  Vec<Object *> neighbors_;

public:
  Object(size_t sid, size_t id) {
    unique_id_ = id;
    species_id_ = sid;
    pos_.resize(_n_dims_max);
  }
  virtual ~Object();

  size_t GetID() { return unique_id_; }
  size_t GetSpeciesID() { return species_id_; }
  double GetPos(int i_dim) { return pos_[i_dim]; }

  void SetPos(int i_dim, double new_pos) { pos_[i_dim] = new_pos; }
};
#endif
