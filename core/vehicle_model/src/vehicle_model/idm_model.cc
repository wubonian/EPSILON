#include "vehicle_model/idm_model.h"

#include "common/math/calculations.h"
#include "odeint-v2/boost/numeric/odeint.hpp"

namespace odeint = boost::numeric::odeint;

namespace simulator {

IntelligentDriverModel::IntelligentDriverModel() { UpdateInternalState(); }

IntelligentDriverModel::IntelligentDriverModel(const Param &parm)
    : param_(parm) {
  UpdateInternalState();
}

IntelligentDriverModel::~IntelligentDriverModel() {}

/* 基于当前的state, 与delta_t, 自车加速度为IDM模型计算结果, 前车加速度为0
   迭代更新内部的state: internal_state_ */
void IntelligentDriverModel::Step(double dt) {
  odeint::integrate(boost::ref(*this), internal_state_, 0.0, dt, dt);
  // Linear(internal_state_, dt, &internal_state_);

  // printf("[Internal]%lf, %lf, %lf, %lf\n", internal_state_[0],
  //        internal_state_[1], internal_state_[2], internal_state_[3]);

  state_.s = internal_state_[0];
  state_.v = internal_state_[1];
  state_.s_front = internal_state_[2];
  state_.v_front = internal_state_[3];
  UpdateInternalState();
}

/* 基于当前的state, 与delta_t, 假设自车加速度为IDM模型算出来的加速度, 前车加速度为0
   计算出delta_t时间后, 自车与前车的state */
void IntelligentDriverModel::Linear(const InternalState &x, const double dt,
                                    InternalState *x_out) {
  State cur_state;
  cur_state.s = x[0];
  cur_state.v = x[1];
  cur_state.s_front = x[2];
  cur_state.v_front = x[3];

  decimal_t acc;
  common::IntelligentDriverModel::GetIIdmDesiredAcceleration(param_, cur_state,
                                                             &acc);

  std::cout << "acc = " << acc << std::endl;
  acc = std::max(acc,
                 -std::min(param_.kHardBrakingDeceleration, cur_state.v / dt));

  std::cout << "acc = " << acc << std::endl;

  (*x_out)[0] = x[0] + cur_state.v * dt + 0.5 * acc * dt * dt;
  (*x_out)[1] = cur_state.v + acc * dt;
  (*x_out)[2] = x[2] + x[3] * dt;
  (*x_out)[3] = x[3];
}

/* 对象函数的重载, 基于当前的state, 计算对应各state的导数 */
void IntelligentDriverModel::operator()(const InternalState &x,
                                        InternalState &dxdt, const double dt) {
  State cur_state;
  cur_state.s = x[0];
  cur_state.v = x[1];
  cur_state.s_front = x[2];
  cur_state.v_front = x[3];

  decimal_t acc;
  common::IntelligentDriverModel::GetAccDesiredAcceleration(param_, cur_state,
                                                            &acc);
  acc = std::max(acc,
                 -std::min(param_.kHardBrakingDeceleration, cur_state.v / dt));
  dxdt[0] = cur_state.v;
  dxdt[1] = acc;
  dxdt[2] = cur_state.v_front;
  dxdt[3] = 0.0;  // assume other vehicle keep the current velocity
}

const IntelligentDriverModel::State &IntelligentDriverModel::state(void) const {
  return state_;
}

void IntelligentDriverModel::set_state(
    const IntelligentDriverModel::State &state) {
  state_ = state;
  UpdateInternalState();
}

void IntelligentDriverModel::UpdateInternalState(void) {
  internal_state_[0] = state_.s;
  internal_state_[1] = state_.v;
  internal_state_[2] = state_.s_front;
  internal_state_[3] = state_.v_front;
}

}  // namespace simulator