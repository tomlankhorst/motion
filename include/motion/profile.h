#ifndef _MOTION_PROFILE_H_
#define _MOTION_PROFILE_H_

#include <cstdint>
#include <cmath>
#include <Eigen/Dense>

namespace motion { namespace profile {

using namespace Eigen;

template<typename T>
class q_profile {
public:
    virtual T q_at(T) = 0;
    virtual T q() = 0;
};

template<typename T>
class v_profile : public q_profile<T> {
    virtual T v_at(T) = 0;
    virtual T v() = 0;
};

template<typename T=float>
class linear : public q_profile<T>
{
  using TS = Matrix<T,2,3>;
  using TC = Matrix<T,2,1>;
  using TT = Matrix<T,1,2>;

  enum rows {
    t_row = 0,
    q_row
  };

  enum cols {
    cur_col = 0,
    nul_col,
    f_col
  };

  // State-matrix
  TS S;

  // Coef-vector
  TC C;

  T time()
  {
    auto t = S(t_row, cur_col);
    auto t0 = S(t_row, nul_col);
    auto tf = S(t_row, f_col);

    auto tu = tf > t0 ? tf : t0;
    auto tl = tf > t0 ? t0 : tf;

    return std::min(tu, std::max(tl, t));
  }

  const TT qs()
  {
    auto t = time();
    TT s;
    s << 1, t;
    return s;
  }

  void upd_coef()
  {
    auto t0 = S(t_row, nul_col);
    auto tf = S(t_row, f_col);
    auto q0 = S(q_row, nul_col);
    auto qf = S(q_row, f_col);

    Matrix<T,2,2> A;
    Matrix<T,2,1> Q;

    A <<  1, t0,
          1, tf;

    Q << q0, qf;

    C = A.colPivHouseholderQr().solve(Q);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  linear()
  {
    // Initialize the state to zeros
    S = TS::Zero();
    C = TC::Zero();
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
  void set(T tf, T qf)
  {
    auto t0 = S(t_row, cur_col);
    auto q0 = q();
    set(t0, tf, q0, qf);
  }
  void set(T t0, T tf, T q0, T qf)
  {
    S.block(0,1,2,2) << t0, tf, q0, qf;
    upd_coef();
  }
};

template<typename T=float>
class cubic : public v_profile<T>
{
  using TS = Matrix<T,3,3>;
  using TC = Matrix<T,4,1>;
  using TT = Matrix<T,1,4>;

  enum rows {
    t_row = 0,
    q_row,
    v_row
  };

  enum cols {
    cur_col = 0,
    nul_col,
    f_col
  };

  /**
   * State matrix
   *
   * t t0 tf
   * q q0 qf
   * v v0 vf
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

    auto tu = tf > t0 ? tf : t0;
    auto tl = tf > t0 ? t0 : tf;

    return std::min(tu, std::max(tl, t));
  }

  const TT qs() 
  {
    auto t = time();
    TT s;
    s << 1, t, std::pow(t,2), std::pow(t,3);
    return s;
  }

  const TT vs() 
  {
    auto t = time();
    TT s;
    s << 0, 1, 2*t, 3*std::pow(t,2);
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

    Matrix<T,4,4> A;
    Matrix<T,4,1> Q;

    A <<  1, t0,  std::pow(t0,2),   std::pow(t0,3),
          0, 1,   2*t0,             3*std::pow(t0,2),
          1, tf,  std::pow(tf,2),   std::pow(tf,3),
          0, 1,   2*tf,             3*std::pow(tf,2);

    Q << q0, v0, qf, vf;

    C = A.colPivHouseholderQr().solve(Q);
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  cubic()
  {
    // Initialize the state to zeros
    S = TS::Zero();
    C = TC::Zero();
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
  void set(T tf, T qf, T vf = 0)
  {
    auto q0 = q();
    auto v0 = v();
    auto t0 = S(t_row, cur_col);

    set(t0, tf, q0, qf, v0, vf);
  }
  void set(T t0, T tf, T q0, T qf, T v0=0, T vf=0)
  {
    S.block(0,1,3,2) << t0, tf, q0, qf, v0, vf;
    upd_coef();
  }
};

}}

#endif // _MOTION_PROFILE_H_
