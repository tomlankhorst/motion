#ifndef _MOTION_PROFILE_H_
#define _MOTION_PROFILE_H_

#include <cstdint>
#include <cmath>
#include <Eigen/Dense>

namespace motion { namespace profile {

using namespace Eigen;

template<typename T=float>
class cubic
{
  using TS = Matrix<T,3,3>;
  using TC = Matrix<T,4,1>;
  using TT = Matrix<T,1,4>;

  enum rows { 
    t_row = 0,
    v_row,
    q_row
  }

  enum cols {
    cur_col = 0,
    nul_col,
    f_col
  }

  /**
   * State matrix
   *
   * t t0 tf
   * v v0 vf
   * q q0 qf
   */
  TS S; 

  /**
   * Coefficients vector
   *
   * a0 a1 a2 a3 .'
   */
  TC C;

  T time()
  {
    auto t = S(t_row, cur_col);
    auto t0 = S(t_row, nul_col);
    auto tf = S(t_row, f_col);
    return std::min(tf, std::max(t0, t));
  }

  const TT qs() 
  {
    auto t = time();
    TT s << 1, t, std::pow(t,2), std::pow(t,3);
    return s;
  }

  const TT vs() 
  {
    auto t = time();
    TT s << 0, 1, 2*t, 3*std::pow(t,2);
    return s;
  }

  void upd_coef()
  {
    auto t0 = S(t_row, nul_col);
    auto tf = S(t_row, f_col);
    auto v0 = S(v_row, nul_col);
    auto vf = S(v_row, f_col);
    auto q0 = S(q_row, nul_col);
    auto qf = S(q_row, f_col);

    Matrix<T,4,4> A <<  1, t0,  std::pow(t0,2),   std::pow(t0,3), 
                        0, 1,   2*t0,             3*std::pow(t0,2), 
                        1, tf,  std::pow(tf,2),   std::pow(tf,3), 
                        0, 1,   2*tf,             3*std::pow(tf,2);

    Matrix<T,4,1> Q << q0, v0, qf, vf;

    C = A.colPivHouseholderQr().solve(Q);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  cubic()
  {
    // Initialize the state to zeros
    S = TS::Zeros();
    C = TC::Zeros();
  }
  T q_at(T t)
  {
    S(t_row,cur_col) = t;
    return q();
  }
  T q()
  {
    auto q = qs()*C;
    return q(0);
  }
  T v_at(T t)
  {
    S(t_row,cur_col) = t;
    return v();
  }
  T v()
  {
    auto v = vs()*C;
    return v(0);
  }
  void set(T t, T qf, T vf = 0)
  {
    auto q0 = q();
    S(q_row, nul_col) = q0;
    auto v0 = v();
    S(v_row, nul_col) = v0;

    S(q_row, f_col) = qf;
    S(v_row, f_col) = vf;

    upd_coef();
  }
};

}}

#endif // _MOTION_PROFILE_H_
