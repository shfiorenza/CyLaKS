#ifndef _CYLAKS_OBJECT_HPP_
#define _CYLAKS_OBJECT_HPP_
#include "system_definitions.hpp"

// Object: Basic building block of CyLaKS; foundation of all proteins/etc
class Object {
protected:
  size_t unique_id_{0}; // Unique ID within CyLaKS among ALL objects
  size_t species_id_{0};

public:
  bool visible_{true};
  Vec<double> pos_; // C.O.M. position in lab frame

  Object *teth_partner_{nullptr};

public:
  Object() {}
  virtual ~Object() {}
  void Initialize(size_t sid, size_t id) {
    unique_id_ = id;
    species_id_ = sid;
    pos_.resize(_n_dims_max);
  }
  size_t GetID() { return unique_id_; }
  size_t GetSpeciesID() { return species_id_; }

  virtual Object *GetTethPartner() { return teth_partner_; }

  virtual bool IsOccupied() { return true; }
  virtual bool IsTethered() {
    return teth_partner_ == nullptr ? false : true;
    // if (teth_partner_ != nullptr) {
    //   return true;
    // }
    // return false;
  }
  virtual bool HasSatellite() { return false; }

  virtual void AddForce(Vec<double> f) {}
  virtual void AddTorque(double tq) {}

  virtual int GetNumNeighborsOccupied() { return -1; }
  virtual int GetNumNeighborsOccupied_Side() { return -1; }

  virtual int GetNumHeadsActive() { return -1; }

  virtual Object *GetHeadOne() { return nullptr; }
  virtual Object *GetHeadTwo() { return nullptr; }
  virtual Object *GetOtherHead() { return nullptr; }
  virtual Object *GetActiveHead() { return nullptr; }
  virtual Object *GetDockedHead() { return nullptr; }

  virtual Vec<double> GetSpringOrientation() { return {}; }
  virtual Vec<double> GetBoundObjectOrientation() { return {}; }
  virtual bool Unbind() { return false; };
};
#endif
