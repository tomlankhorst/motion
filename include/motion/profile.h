#ifndef _MOTION_PROFILE_H_
#define _MOTION_PROFILE_H_

#include <cstdint>
#include <cmath>
#include <Eigen/Dense>

namespace motion { namespace profile {

using namespace Eigen;

/**
 * @enum S-row t(ime), v(elocity), a(cceleration)
 */
enum Sr { rt=0, rq, rv, ra };

/**
 * @enum S-col c(urrent), 0(ul), f(inal)
 */
enum Sc { cc=0, c0,cf };

template <typename T, size_t order>
class profile {
protected:
  // First order has one state (q), third has two (q,v)...
  using TS = Matrix<T,1+(order+1)/2,3>;

  // Coefficients is order+1
  using TC = Matrix<T,order+1,1>;

  // Phi vector has order+1 values
  using TT = Matrix<T,1,order+1>;

  /**
   * State matrix with `col` and `row` indices
   */
  TS S;

  /**
   * Coefficient vector
   */
  TC C;

  /**
   * Current time within range
   *
   * @return T time within range t0 tf
   */
  T tr()
  {
    auto t = S(rt, cc);
    auto t0 = S(rt,c0);
    auto tf = S(rt,cf);

    auto tu = tf > t0 ? tf : t0;
    auto tl = tf > t0 ? t0 : tf;

    return std::min(tu, std::max(tl, t));
  }

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /**
   * Initialize state and coefficients to zero
   */
  profile()
  {
    S = TS::Zero();
    C = TC::Zero();
  }

  /**
   * Takes a time and returns the appropriate value
   * @param t time
   * @return T value
   */
  virtual T q_at(T t) = 0;

  /**
   * Returns the value for the set time
   * @return T value
   */
  virtual T q() = 0;
};

/**
 * Linear motion profile
 *
 * @tparam T datatype
 */
template<typename T=float>
class linear : public profile<T, 1>
{
  using p = profile<T,1>;
  using typename p::TT;
  using p::S;
  using p::C;

  const TT qs()
  {
    TT s;
    s << 1, p::tr();
    return s;
  }

  void upd_coef()
  {
    Matrix<T,2,2> A;
    Matrix<T,2,1> Q;

    auto t0 = S(rt,c0);
    auto tf = S(rt,cf);

    A <<  1, t0,
          1, tf;

    Q << S(rq,c0), S(rq,cf);

    C = A.colPivHouseholderQr().solve(Q);
  }

public:
  /**
   * @inherit
   */
  T q_at(T t)
  {
    S(rt,cc) = t;
    return q();
  }

  /**
   * @inherit
   */
  T q()
  {
    auto q = qs()*C;
    return q(0);
  }

  /**
   * Set a new final time and q
   * @param tf
   * @param qf
   */
  void set(T tf, T qf)
  {
    auto t0 = S(rt, cc);
    auto q0 = q();
    set(t0, tf, q0, qf);
  }

  /**
   * Set both a final and start time and q
   * @param t0
   * @param tf
   * @param q0
   * @param qf
   */
  void set(T t0, T tf, T q0, T qf)
  {
    S.block(0,1,2,2) << t0, tf, q0, qf;
    upd_coef();
  }

};

/**
 * Cubic motion profile
 *
 * @tparam T datatype
 */
template<typename T=float>
class cubic : public profile<T,3>
{
  using p = profile<T,3>;
  using typename p::TT;
  using p::S;
  using p::C;

  const TT qs()
  {
    TT s;
    auto t = p::tr();
    s << 1, t, std::pow(t,2), std::pow(t,3);
    return s;
  }

  const TT vs() 
  {
    TT s;
    auto t = p::tr();
    s << 0, 1, 2*t, 3*std::pow(t,2);
    return s;
  }

  void upd_coef()
  {
    Matrix<T,4,4> A;
    Matrix<T,4,1> Q;

    auto t0 = S(rt,c0);
    auto tf = S(rt,cf);

    A <<  1, t0,  std::pow(t0,2),   std::pow(t0,3),
          0, 1,   2*t0,             3*std::pow(t0,2),
          1, tf,  std::pow(tf,2),   std::pow(tf,3),
          0, 1,   2*tf,             3*std::pow(tf,2);

    Q << S(rq,c0), S(rv,c0), S(rq,cf), S(rv,cf);

    C = A.colPivHouseholderQr().solve(Q);
  }

public:

  /**
   * @inherit
   */
  T q_at(T t)
  {
    S(rt,cc) = t;
    return q();
  }

  /**
   * @inherit
   */
  T q()
  {
    auto q = qs()*C;
    return q(0);
  }

  /**
   * Get v at t
   * @param t time
   * @return q
   */
  T v_at(T t)
  {
    S(rt,cc) = t;
    return v();
  }

  /**
   * Get v at current t
   * @return v
   */
  T v()
  {
    auto v = vs()*C;
    return v(0);
  }

  /**
   * Set new final t, q, v with current time as zero
   * @param tf
   * @param qf
   * @param vf
   */
  void set(T tf, T qf, T vf = 0)
  {
    set(S(rt,cc), tf, q(), qf, v(), vf);
  }

  /**
   * Set new null and final t, q, v
   * @param t0
   * @param tf
   * @param q0
   * @param qf
   * @param v0
   * @param vf
   */
  void set(T t0, T tf, T q0, T qf, T v0=0, T vf=0)
  {
    S.block(0,1,3,2) << t0, tf, q0, qf, v0, vf;
    upd_coef();
  }
};

}}

#endif // _MOTION_PROFILE_H_
