#pragma once

// the ode code is taken from odient library:
// https://github.com/headmyshoulder/odeint-v2

#include <array>
#include <deal.II/lac/vector.h>
#include <vector>

template <class state_type>
inline void resize(const state_type &in, state_type &out) {
  // standard implementation works for containers
  out.resize(in.size());
}

// specialization for std::array
template <class T, std::size_t N>
inline void resize(const std::array<T, N> &, std::array<T, N> &) {
  /* arrays don't need resizing */
}

template <typename T>
inline void resize(const dealii::Vector<T> &in, dealii::Vector<T> &out) {
  out.reinit(in.size());
}

struct container_algebra {
  template <class S1, class S2, class S3, class Op>
  void for_each3(S1 &s1, S2 &s2, S3 &s3, Op op) const {
    using std::begin;
    using std::end;

    auto first1 = begin(s1);
    auto last1 = end(s1);
    auto first2 = begin(s2);
    auto first3 = begin(s3);
    for (; first1 != last1;)
      op(*first1++, *first2++, *first3++);
  }

  template <class S1, class S2, class S3, class S4, class S5, class S6,
            class Op>
  void for_each6(S1 &s1, S2 &s2, S3 &s3, S4 &s4, S5 &s5, S6 &s6, Op op) const {
    using std::begin;
    using std::end;

    auto first1 = begin(s1);
    auto last1 = end(s1);
    auto first2 = begin(s2);
    auto first3 = begin(s3);
    auto first4 = begin(s4);
    auto first5 = begin(s5);
    auto first6 = begin(s6);
    for (; first1 != last1;)
      op(*first1++, *first2++, *first3++, *first4++, *first5++, *first6++);
  }
};

struct default_operations {

  template <class F1 = double, class F2 = F1> struct scale_sum2 {
    typedef void result_type;

    const F1 alpha1;
    const F2 alpha2;

    scale_sum2(F1 a1, F2 a2) : alpha1(a1), alpha2(a2) {}

    template <class T0, class T1, class T2>
    void operator()(T0 &t0, const T1 &t1, const T2 &t2) const {
      t0 = alpha1 * t1 + alpha2 * t2;
    }
  };

  template <class F1 = double, class F2 = F1, class F3 = F2, class F4 = F3,
            class F5 = F4>
  struct scale_sum5 {
    typedef void result_type;

    const F1 alpha1;
    const F2 alpha2;
    const F3 alpha3;
    const F4 alpha4;
    const F5 alpha5;

    scale_sum5(F1 a1, F2 a2, F3 a3, F4 a4, F5 a5)
        : alpha1(a1), alpha2(a2), alpha3(a3), alpha4(a4), alpha5(a5) {}

    template <class T0, class T1, class T2, class T3, class T4, class T5>
    void operator()(T0 &t0, const T1 &t1, const T2 &t2, const T3 &t3,
                    const T4 &t4, const T5 &t5) const {
      t0 = alpha1 * t1 + alpha2 * t2 + alpha3 * t3 + alpha4 * t4 + alpha5 * t5;
    }
  };
};

//euler
// TODO add operation and container:
template <class state_type, class value_type = double,
          class time_type = value_type>
class euler {
public:
  template <typename System>
  void do_step(System &system, state_type &x, time_type t, time_type dt) {
    resize(x, k);
    system(x, k, t);
    for (unsigned int i = 0; i < x.size(); ++i)
      x[i] += dt * k[i];
  }

private:
  state_type k;
};

// rungekutta method:
template <class state_type, class value_type = double,
          class time_type = value_type, class algebra = container_algebra,
          class operations = default_operations>
class runge_kutta4 {
public:
  template <typename System>
  void do_step(System &system, state_type &x, time_type t, time_type dt) {
    adjust_size(x);
    const value_type one = 1;
    const time_type dt2 = dt / 2, dt3 = dt / 3, dt6 = dt / 6;

    typedef typename operations::template scale_sum2<value_type, time_type>
        scale_sum2;

    typedef typename operations::template scale_sum5<
        value_type, time_type, time_type, time_type, time_type>
        scale_sum5;

    system(x, k1, t);
    m_algebra.for_each3(x_tmp, x, k1, scale_sum2(one, dt2));

    system(x_tmp, k2, t + dt2);
    m_algebra.for_each3(x_tmp, x, k2, scale_sum2(one, dt2));

    system(x_tmp, k3, t + dt2);
    m_algebra.for_each3(x_tmp, x, k3, scale_sum2(one, dt));

    system(x_tmp, k4, t + dt);
    m_algebra.for_each6(x, x, k1, k2, k3, k4,
                        scale_sum5(one, dt6, dt3, dt3, dt6));
  }

private:
  state_type x_tmp, k1, k2, k3, k4;
  algebra m_algebra;

  void adjust_size(const state_type &x) {
    resize(x, x_tmp);
    resize(x, k1);
    resize(x, k2);
    resize(x, k3);
    resize(x, k4);
  }
};
