#ifndef _BASE_OBJECT
#define _BASE_OBJECT
#include <vector>
template <typename DATA_T> using Vec = std::vector<DATA_T>;

class BaseObject {
private:
  unsigned int id_;
  double pos_;
  bool visible_;

public:
  BaseObject(int id, bool visible, double pos) {
    id_ = id;
    visible_ = visible;
    pos_ = pos;
  }
  BaseObject(int id) {
    id_ = id;
    visible_ = true;
    pos_ = 0.0;
  }

  int GetID() { return id_; }

  double GetPos() { return pos_; }
  void SetPos(double new_pos) { pos_ = new_pos; }
};
#endif