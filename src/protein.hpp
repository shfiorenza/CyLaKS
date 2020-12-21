#ifndef _CYLAKS_PROTEIN_HPP_
#define _CYLAKS_PROTEIN_HPP_
#include "binding_head.hpp"
#include "linear_spring.hpp"
#include "system_namespace.hpp"
// FIXME shouldnt include this
#include "binding_site.hpp"
#include "protofilament.hpp"

class Protein : public Object {
protected:
  int dist_cutoff_{10};
  int n_neighbors_bind_ii_{0};
  Vec<BindingSite *> neighbors_bind_ii_;

public:
  size_t active_index_{0};
  size_t n_heads_active_{0};
  BindingHead head_one_, head_two_;
  LinearSpring spring_;

  bool tethered_{false};
  Protein *partner_{nullptr};

private:
  void InitializeNeighborList();

public:
  Protein() {}
  void Initialize(size_t sid, size_t id) {
    Object::Initialize(sid, id);
    head_one_.Initialize(sid, id, _r_xlink_head, this, &head_two_);
    head_two_.Initialize(sid, id, _r_xlink_head, this, &head_one_);
    spring_.Initialize(sid, id, this);
    neighbors_bind_ii_.resize(2 * dist_cutoff_ + 1);
  }

  bool HasSatellite();
  void UntetherSatellite();

  int GetNumHeadsActive() { return n_heads_active_; }
  BindingHead *GetActiveHead();
  BindingHead *GetHeadOne() { return &head_one_; }
  BindingHead *GetHeadTwo() { return &head_two_; }

  // void UpdateNeighborList();
  // Object *GetWeightedNeighbor();

  void UpdateNeighbors_Bind_II();
  double GetWeight_Bind_II(BindingSite *neighb);

  double GetTotalWeight_Bind_I_Teth();
  double GetTotalWeight_Bind_II();
  double GetTotalWeight_Bind_II_Teth();
  BindingSite *GetNeighbor_Bind_I_Teth();
  BindingSite *GetNeighbor_Bind_II();
  BindingSite *GetNeighbor_Bind_II_Teth();

  virtual double GetAnchorCoordinate(int i_dim) {

    if (n_heads_active_ != 2) {
      Sys::ErrorExit("Protein::GetAnchorCoord()");
    }
    return (head_one_.site_->pos_[i_dim] + head_two_.site_->pos_[i_dim]) / 2;
  }
  virtual int GetDirectionTowardRest(BindingHead *head) {
    if (n_heads_active_ == 1) {
      return 1;
    } else if (n_heads_active_ == 2) {
      // printf("%g > %g?\n", head->site_->pos_[0], GetAnchorCoordinate(0));
      if (head->site_->pos_[0] > GetAnchorCoordinate(0)) {
        return -1;
      } else if (head->site_->pos_[0] < GetAnchorCoordinate(0)) {
        return 1;
      } else {
        return 0;
      }
    } else {
      Sys::ErrorExit("Protein::GetDirToRest()\n");
    }
    return 0;
  }

  virtual void UpdateExtension() {

    if (n_heads_active_ != 2) {
      return;
    }
    double r_x{head_one_.site_->pos_[0] - head_two_.site_->pos_[0]};
    double r_y{head_one_.site_->pos_[1] - head_two_.site_->pos_[1]};
    double r{sqrt(Square(r_x) + Square(r_y))};
    double sine{r_y / r};
    double cosine{r_x / r};
    double dr{r - Params::Xlinks::r_0};
    double f_mag{fabs(dr) * Params::Xlinks::k_spring};
    Vec<double> r_vec{r_x, r_y};

    // Add forces
    double anchor_x{GetAnchorCoordinate(0)};
    if (head_one_.site_->pos_[0] < anchor_x) {
      head_one_.site_->filament_->force_[0] += f_mag * cosine;
      head_two_.site_->filament_->force_[0] += -f_mag * cosine;
    } else {
      head_one_.site_->filament_->force_[0] += -f_mag * cosine;
      head_two_.site_->filament_->force_[0] += f_mag * cosine;
    }
    double anchor_y{GetAnchorCoordinate(1)};
    if (head_one_.site_->pos_[1] < anchor_y) {
      head_one_.site_->filament_->force_[1] += f_mag * sine;
      head_two_.site_->filament_->force_[1] += -f_mag * sine;
    } else {
      head_one_.site_->filament_->force_[1] += -f_mag * sine;
      head_two_.site_->filament_->force_[1] += f_mag * sine;
    }
    // Add torques
    double rod_com_x1{head_one_.site_->filament_->pos_[0]};
    double rod_com_y1{head_one_.site_->filament_->pos_[1]};
    double d_x1{head_one_.site_->pos_[0] - rod_com_x1};
    double d_y1{head_one_.site_->pos_[1] - rod_com_y1};
    double lever_arm1{sqrt(Square(d_x1) + Square(d_y1))};
    if (head_one_.site_->pos_[0] > rod_com_x1) {
      if (anchor_y > rod_com_y1) {
        head_one_.site_->filament_->torque_ += f_mag * lever_arm1 * sine;
      } else {
        head_one_.site_->filament_->torque_ += -f_mag * lever_arm1 * sine;
      }
    } else {
      if (anchor_y > rod_com_y1) {
        head_one_.site_->filament_->torque_ += -f_mag * lever_arm1 * sine;
      } else {
        head_one_.site_->filament_->torque_ += f_mag * lever_arm1 * sine;
      }
    }
    double rod_com_x2{head_two_.site_->filament_->pos_[0]};
    double rod_com_y2{head_two_.site_->filament_->pos_[1]};
    double d_x2{head_two_.site_->pos_[0] - rod_com_x2};
    double d_y2{head_two_.site_->pos_[1] - rod_com_y2};
    double lever_arm2{sqrt(Square(d_x2) + Square(d_y2))};
    if (head_two_.site_->pos_[0] > rod_com_x2) {
      if (anchor_y > rod_com_y2) {
        head_two_.site_->filament_->torque_ += f_mag * lever_arm2 * sine;
      } else {
        head_two_.site_->filament_->torque_ += -f_mag * lever_arm2 * sine;
      }
    } else {
      if (anchor_y > rod_com_y2) {
        head_two_.site_->filament_->torque_ += -f_mag * lever_arm2 * sine;
      } else {
        head_two_.site_->filament_->torque_ += f_mag * lever_arm2 * sine;
      }
    }
  }

  virtual bool Bind(BindingSite *site, BindingHead *head);
  virtual bool Unbind(BindingHead *head);

  virtual double GetWeight_Diffuse(BindingHead *head, int dir);
  virtual bool Diffuse(BindingHead *head, int dir);

  virtual bool Tether();
  virtual bool Untether();
};

#endif