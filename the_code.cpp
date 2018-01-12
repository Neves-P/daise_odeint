// Copyright Timothy H. Keitt 2015
// See license for odeintr package

// [[Rcpp::depends(odeintr)]]

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include "boost/numeric/odeint.hpp"
namespace odeint = boost::numeric::odeint;

;

namespace odeintr
{
  static const std::size_t N = 9;

  typedef std::vector<double> state_type;
  
  static state_type state(N);
  
  typedef odeint::runge_kutta_dopri5<state_type> stepper_type;
  
  static auto stepper = odeint::make_dense_output(1e-06, 1e-06, stepper_type());
  
  typedef std::vector<double> vec_type;
  static std::vector<vec_type> rec_x(N);
  static vec_type rec_t;
  
  double gamvec_il3_1 = 0.2, gamvec_il3_one = 0.2, laavec_il1_plusone_1 = 0.3, laavec_il3_one = 0.3, laavec_il3_plusone_1 = 0.3, lacvec_il1_1 = 0.4, lacvec_il1_plusone_1 = 0.4, lacvec_il3_1 = 0.4, lacvec_il3_one = 0.4, lacvec_il3_plusone_1 = 0.4, lacvec_il4_plusone_1 = 0.4, muvec_il2_plusone_1 = 0.1, muvec_il3_1 = 0.1, muvec_il3_one = 0.1, muvec_il3_plusone_1 = 0.1, nn_in1_1 = 19, nn_in2_1 = 1, nn_in3_plusone_1 = 11, gamvec_il3_2 = 0.2, laavec_il1_plusone_2 = 0.3, laavec_il3_plusone_2 = 0.3, lacvec_il1_2 = 0.4, lacvec_il1_plusone_2 = 0.4, lacvec_il3_2 = 0.4, lacvec_il3_plusone_2 = 0.4, lacvec_il4_plusone_2 = 0.4, muvec_il2_plusone_2 = 0.1, muvec_il3_2 = 0.1, muvec_il3_plusone_2 = 0.1, nn_in1_2 = 20, nn_in2_2 = 2, nn_in3_plusone_2 = 12, gamvec_il3_3 = 0.2, laavec_il1_plusone_3 = 0.3, laavec_il3_plusone_3 = 0.3, lacvec_il1_3 = 0.4, lacvec_il1_plusone_3 = 0.4, lacvec_il3_3 = 0.4, lacvec_il3_plusone_3 = 0.4, lacvec_il4_plusone_3 = 0.4, muvec_il2_plusone_3 = 0.1, muvec_il3_3 = 0.1, muvec_il3_plusone_3 = 0.1, nn_in1_3 = 21, nn_in2_3 = 3, nn_in3_plusone_3 = 13, gamvec_il3_4 = 0.2, laavec_il1_plusone_4 = 0.3, laavec_il3_plusone_4 = 0.3, lacvec_il1_4 = 0.4, lacvec_il1_plusone_4 = 0.4, lacvec_il3_4 = 0.4, lacvec_il3_plusone_4 = 0.4, lacvec_il4_plusone_4 = 0.4, muvec_il2_plusone_4 = 0.1, muvec_il3_4 = 0.1, muvec_il3_plusone_4 = 0.1, nn_in1_4 = 22, nn_in2_4 = 4, nn_in3_plusone_4 = 14;;
  
  #include "utils.h"
  
  static void
  sys(const state_type x, state_type &dxdt, const double t)
  {
    dxdt[0] = laavec_il1_plusone_1 * xx2_ix1_1 + lacvec_il4_plusone_1 * xx2_ix4_1 + muvec_il2_plusone_1 * xx2_ix3_1 + lacvec_il1_1 * nn_in1_1 * xx1_ix1_1 + -(muvec_il3_1 + lacvec_il3_1) * nn_in2_1 * xx1_ix3_1 + gamvec_il3_1 * xx1_ix3_1; dxdt[1] = laavec_il1_plusone_2 * xx2_ix1_2 + lacvec_il4_plusone_2 * xx2_ix4_2 + muvec_il2_plusone_2 * xx2_ix3_2 + lacvec_il1_2 * nn_in1_2 * xx1_ix1_2 + -(muvec_il3_2 + lacvec_il3_2) * nn_in2_2 * xx1_ix3_2 + gamvec_il3_2 * xx1_ix3_2; dxdt[2] = laavec_il1_plusone_3 * xx2_ix1_3 + lacvec_il4_plusone_3 * xx2_ix4_3 + muvec_il2_plusone_3 * xx2_ix3_3 + lacvec_il1_3 * nn_in1_3 * xx1_ix1_3 + -(muvec_il3_3 + lacvec_il3_3) * nn_in2_3 * xx1_ix3_3 + gamvec_il3_3 * xx1_ix3_3; dxdt[3] = laavec_il1_plusone_4 * xx2_ix1_4 + lacvec_il4_plusone_4 * xx2_ix4_4 + muvec_il2_plusone_4 * xx2_ix3_4 + lacvec_il1_4 * nn_in1_4 * xx1_ix1_4 + -(muvec_il3_4 + lacvec_il3_4) * nn_in2_4 * xx1_ix3_4 + gamvec_il3_4 * xx1_ix3_4; dxdt[4] = gamvec_il3_1 * xx1_ix3_1 + lacvec_il1_plusone_1 * nn_in1_1 * xx2_ix1_1 + muvec_il2_plusone_1 * nn_in2_1 * xx2_ix2_1 + -(muvec_il3_plusone_1 + lacvec_il3_plusone_1) * nn_in3_plusone_1 * xx2_ix3_1 + laavec_il3_plusone_1 * xx2_ix3_1; dxdt[5] = gamvec_il3_2 * xx1_ix3_2 + lacvec_il1_plusone_2 * nn_in1_2 * xx2_ix1_2 + muvec_il2_plusone_2 * nn_in2_2 * xx2_ix2_2 + -(muvec_il3_plusone_2 + lacvec_il3_plusone_2) * nn_in3_plusone_2 * xx2_ix3_2 + laavec_il3_plusone_2 * xx2_ix3_2; dxdt[6] = gamvec_il3_3 * xx1_ix3_3 + lacvec_il1_plusone_3 * nn_in1_3 * xx2_ix1_3 + muvec_il2_plusone_3 * nn_in2_3 * xx2_ix2_3 + -(muvec_il3_plusone_3 + lacvec_il3_plusone_3) * nn_in3_plusone_3 * xx2_ix3_3 + laavec_il3_plusone_3 * xx2_ix3_3; dxdt[7] = gamvec_il3_4 * xx1_ix3_4 + lacvec_il1_plusone_4 * nn_in1_4 * xx2_ix1_4 + muvec_il2_plusone_4 * nn_in2_4 * xx2_ix2_4 + -(muvec_il3_plusone_4 + lacvec_il3_plusone_4) * nn_in3_plusone_4 * xx2_ix3_4 + laavec_il3_plusone_4 * xx2_ix3_4; dxdt[8] = -(laavec_il3_one + lacvec_il3_one + gamvec_il3_one + muvec_il3_one) * xx3; ;
  }

  static void
  obs(const state_type x, const double t)
  {
    for (int i = 0; i != N; ++i)
      rec_x[i].push_back(x[i]);
    rec_t.push_back(t);
  }
  
}; // namespace odeintr

static void
reserve(odeintr::vec_type::size_type n)
{
  odeintr::rec_t.reserve(n);
  for (auto &i : odeintr::rec_x) i.reserve(n);
}

// [[Rcpp::export]]
Rcpp::List y_get_output()
{
  Rcpp::List out;
  out("Time") = Rcpp::wrap(odeintr::rec_t);
  for (int i = 0; i != odeintr::N; ++i)
  {
    auto cnam = std::string("X") + std::to_string(i + 1);
    out(cnam) = Rcpp::wrap(odeintr::rec_x[i]);
  }
  out.attr("class") = "data.frame";
  int rows_out = odeintr::rec_t.size();
  auto rn = Rcpp::IntegerVector::create(NA_INTEGER, -rows_out);
  out.attr("row.names") = rn;
  return out;
};

// [[Rcpp::export]]
void y_set_state(Rcpp::NumericVector new_state)
{
  if (new_state.size() != odeintr::N)
    Rcpp::stop("Invalid initial state");
  std::copy(new_state.begin(),
            new_state.end(),
            odeintr::state.begin());
}

// [[Rcpp::export]]
std::vector<double>
y_get_state()
{
  return odeintr::state;
}

// [[Rcpp::export]]
void y_reset_observer()
{
  for (auto &i : odeintr::rec_x) i.resize(0);
  odeintr::rec_t.resize(0);  
}

// [[Rcpp::export]]
Rcpp::List y_adap(Rcpp::NumericVector init,
                             double duration,
                             double step_size = 1.0,
                             double start = 0.0)
{
  y_set_state(init);
  y_reset_observer(); reserve(duration / step_size);
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size,
                             odeintr::obs);
  return y_get_output();
}

// [[Rcpp::export]]
Rcpp::List y_at(Rcpp::NumericVector init,
                           std::vector<double> times,
                           double step_size = 1.0,
                           double start = 0.0)
{
  y_set_state(init);
  y_reset_observer(); reserve(times.size());
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, times[0], step_size);
  odeint::integrate_times(odeintr::stepper, odeintr::sys, odeintr::state,
                          times.begin(), times.end(), step_size, odeintr::obs);
  return y_get_output();
}

// [[Rcpp::export]]
Rcpp::List
y_continue_at(std::vector<double> times, double step_size = 1.0)
{
  double start = odeintr::rec_t.back();
  y_reset_observer(); reserve(odeintr::rec_t.size() + times.size());
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, times[0], step_size);
  odeint::integrate_times(odeintr::stepper, odeintr::sys, odeintr::state,
                          times.begin(), times.end(), step_size, odeintr::obs);
  return y_get_output();
}

// [[Rcpp::export]]
Rcpp::List y(Rcpp::NumericVector init,
                       double duration,
                       double step_size = 1.0,
                       double start = 0.0)
{
  y_set_state(init);
  y_reset_observer(); reserve(duration / step_size);
  odeint::integrate_const(odeintr::stepper, odeintr::sys, odeintr::state,
                          start, start + duration, step_size,
                          odeintr::obs);
  return y_get_output();
}

// [[Rcpp::export]]
std::vector<double>
y_no_record(Rcpp::NumericVector init,
                       double duration,
                       double step_size = 1.0,
                       double start = 0.0)
{
  y_set_state(init);
  odeint::integrate_adaptive(odeintr::stepper, odeintr::sys, odeintr::state,
                             start, start + duration, step_size);
  return y_get_state();
}

// [[Rcpp::export]]
void y_set_params(double gamvec_il3_1 = 0.2, double gamvec_il3_one = 0.2, double laavec_il1_plusone_1 = 0.3, double laavec_il3_one = 0.3, double laavec_il3_plusone_1 = 0.3, double lacvec_il1_1 = 0.4, double lacvec_il1_plusone_1 = 0.4, double lacvec_il3_1 = 0.4, double lacvec_il3_one = 0.4, double lacvec_il3_plusone_1 = 0.4, double lacvec_il4_plusone_1 = 0.4, double muvec_il2_plusone_1 = 0.1, double muvec_il3_1 = 0.1, double muvec_il3_one = 0.1, double muvec_il3_plusone_1 = 0.1, double nn_in1_1 = 19, double nn_in2_1 = 1, double nn_in3_plusone_1 = 11, double gamvec_il3_2 = 0.2, double laavec_il1_plusone_2 = 0.3, double laavec_il3_plusone_2 = 0.3, double lacvec_il1_2 = 0.4, double lacvec_il1_plusone_2 = 0.4, double lacvec_il3_2 = 0.4, double lacvec_il3_plusone_2 = 0.4, double lacvec_il4_plusone_2 = 0.4, double muvec_il2_plusone_2 = 0.1, double muvec_il3_2 = 0.1, double muvec_il3_plusone_2 = 0.1, double nn_in1_2 = 20, double nn_in2_2 = 2, double nn_in3_plusone_2 = 12, double gamvec_il3_3 = 0.2, double laavec_il1_plusone_3 = 0.3, double laavec_il3_plusone_3 = 0.3, double lacvec_il1_3 = 0.4, double lacvec_il1_plusone_3 = 0.4, double lacvec_il3_3 = 0.4, double lacvec_il3_plusone_3 = 0.4, double lacvec_il4_plusone_3 = 0.4, double muvec_il2_plusone_3 = 0.1, double muvec_il3_3 = 0.1, double muvec_il3_plusone_3 = 0.1, double nn_in1_3 = 21, double nn_in2_3 = 3, double nn_in3_plusone_3 = 13, double gamvec_il3_4 = 0.2, double laavec_il1_plusone_4 = 0.3, double laavec_il3_plusone_4 = 0.3, double lacvec_il1_4 = 0.4, double lacvec_il1_plusone_4 = 0.4, double lacvec_il3_4 = 0.4, double lacvec_il3_plusone_4 = 0.4, double lacvec_il4_plusone_4 = 0.4, double muvec_il2_plusone_4 = 0.1, double muvec_il3_4 = 0.1, double muvec_il3_plusone_4 = 0.1, double nn_in1_4 = 22, double nn_in2_4 = 4, double nn_in3_plusone_4 = 14)
{ 
  odeintr::gamvec_il3_1 = gamvec_il3_1;
odeintr::gamvec_il3_one = gamvec_il3_one;
odeintr::laavec_il1_plusone_1 = laavec_il1_plusone_1;
odeintr::laavec_il3_one = laavec_il3_one;
odeintr::laavec_il3_plusone_1 = laavec_il3_plusone_1;
odeintr::lacvec_il1_1 = lacvec_il1_1;
odeintr::lacvec_il1_plusone_1 = lacvec_il1_plusone_1;
odeintr::lacvec_il3_1 = lacvec_il3_1;
odeintr::lacvec_il3_one = lacvec_il3_one;
odeintr::lacvec_il3_plusone_1 = lacvec_il3_plusone_1;
odeintr::lacvec_il4_plusone_1 = lacvec_il4_plusone_1;
odeintr::muvec_il2_plusone_1 = muvec_il2_plusone_1;
odeintr::muvec_il3_1 = muvec_il3_1;
odeintr::muvec_il3_one = muvec_il3_one;
odeintr::muvec_il3_plusone_1 = muvec_il3_plusone_1;
odeintr::nn_in1_1 = nn_in1_1;
odeintr::nn_in2_1 = nn_in2_1;
odeintr::nn_in3_plusone_1 = nn_in3_plusone_1;
odeintr::gamvec_il3_2 = gamvec_il3_2;
odeintr::laavec_il1_plusone_2 = laavec_il1_plusone_2;
odeintr::laavec_il3_plusone_2 = laavec_il3_plusone_2;
odeintr::lacvec_il1_2 = lacvec_il1_2;
odeintr::lacvec_il1_plusone_2 = lacvec_il1_plusone_2;
odeintr::lacvec_il3_2 = lacvec_il3_2;
odeintr::lacvec_il3_plusone_2 = lacvec_il3_plusone_2;
odeintr::lacvec_il4_plusone_2 = lacvec_il4_plusone_2;
odeintr::muvec_il2_plusone_2 = muvec_il2_plusone_2;
odeintr::muvec_il3_2 = muvec_il3_2;
odeintr::muvec_il3_plusone_2 = muvec_il3_plusone_2;
odeintr::nn_in1_2 = nn_in1_2;
odeintr::nn_in2_2 = nn_in2_2;
odeintr::nn_in3_plusone_2 = nn_in3_plusone_2;
odeintr::gamvec_il3_3 = gamvec_il3_3;
odeintr::laavec_il1_plusone_3 = laavec_il1_plusone_3;
odeintr::laavec_il3_plusone_3 = laavec_il3_plusone_3;
odeintr::lacvec_il1_3 = lacvec_il1_3;
odeintr::lacvec_il1_plusone_3 = lacvec_il1_plusone_3;
odeintr::lacvec_il3_3 = lacvec_il3_3;
odeintr::lacvec_il3_plusone_3 = lacvec_il3_plusone_3;
odeintr::lacvec_il4_plusone_3 = lacvec_il4_plusone_3;
odeintr::muvec_il2_plusone_3 = muvec_il2_plusone_3;
odeintr::muvec_il3_3 = muvec_il3_3;
odeintr::muvec_il3_plusone_3 = muvec_il3_plusone_3;
odeintr::nn_in1_3 = nn_in1_3;
odeintr::nn_in2_3 = nn_in2_3;
odeintr::nn_in3_plusone_3 = nn_in3_plusone_3;
odeintr::gamvec_il3_4 = gamvec_il3_4;
odeintr::laavec_il1_plusone_4 = laavec_il1_plusone_4;
odeintr::laavec_il3_plusone_4 = laavec_il3_plusone_4;
odeintr::lacvec_il1_4 = lacvec_il1_4;
odeintr::lacvec_il1_plusone_4 = lacvec_il1_plusone_4;
odeintr::lacvec_il3_4 = lacvec_il3_4;
odeintr::lacvec_il3_plusone_4 = lacvec_il3_plusone_4;
odeintr::lacvec_il4_plusone_4 = lacvec_il4_plusone_4;
odeintr::muvec_il2_plusone_4 = muvec_il2_plusone_4;
odeintr::muvec_il3_4 = muvec_il3_4;
odeintr::muvec_il3_plusone_4 = muvec_il3_plusone_4;
odeintr::nn_in1_4 = nn_in1_4;
odeintr::nn_in2_4 = nn_in2_4;
odeintr::nn_in3_plusone_4 = nn_in3_plusone_4;
}
// [[Rcpp::export]]
Rcpp::List y_get_params()
{
  Rcpp::List out;
  out["gamvec_il3_1"] = odeintr::gamvec_il3_1;
out["gamvec_il3_one"] = odeintr::gamvec_il3_one;
out["laavec_il1_plusone_1"] = odeintr::laavec_il1_plusone_1;
out["laavec_il3_one"] = odeintr::laavec_il3_one;
out["laavec_il3_plusone_1"] = odeintr::laavec_il3_plusone_1;
out["lacvec_il1_1"] = odeintr::lacvec_il1_1;
out["lacvec_il1_plusone_1"] = odeintr::lacvec_il1_plusone_1;
out["lacvec_il3_1"] = odeintr::lacvec_il3_1;
out["lacvec_il3_one"] = odeintr::lacvec_il3_one;
out["lacvec_il3_plusone_1"] = odeintr::lacvec_il3_plusone_1;
out["lacvec_il4_plusone_1"] = odeintr::lacvec_il4_plusone_1;
out["muvec_il2_plusone_1"] = odeintr::muvec_il2_plusone_1;
out["muvec_il3_1"] = odeintr::muvec_il3_1;
out["muvec_il3_one"] = odeintr::muvec_il3_one;
out["muvec_il3_plusone_1"] = odeintr::muvec_il3_plusone_1;
out["nn_in1_1"] = odeintr::nn_in1_1;
out["nn_in2_1"] = odeintr::nn_in2_1;
out["nn_in3_plusone_1"] = odeintr::nn_in3_plusone_1;
out["gamvec_il3_2"] = odeintr::gamvec_il3_2;
out["laavec_il1_plusone_2"] = odeintr::laavec_il1_plusone_2;
out["laavec_il3_plusone_2"] = odeintr::laavec_il3_plusone_2;
out["lacvec_il1_2"] = odeintr::lacvec_il1_2;
out["lacvec_il1_plusone_2"] = odeintr::lacvec_il1_plusone_2;
out["lacvec_il3_2"] = odeintr::lacvec_il3_2;
out["lacvec_il3_plusone_2"] = odeintr::lacvec_il3_plusone_2;
out["lacvec_il4_plusone_2"] = odeintr::lacvec_il4_plusone_2;
out["muvec_il2_plusone_2"] = odeintr::muvec_il2_plusone_2;
out["muvec_il3_2"] = odeintr::muvec_il3_2;
out["muvec_il3_plusone_2"] = odeintr::muvec_il3_plusone_2;
out["nn_in1_2"] = odeintr::nn_in1_2;
out["nn_in2_2"] = odeintr::nn_in2_2;
out["nn_in3_plusone_2"] = odeintr::nn_in3_plusone_2;
out["gamvec_il3_3"] = odeintr::gamvec_il3_3;
out["laavec_il1_plusone_3"] = odeintr::laavec_il1_plusone_3;
out["laavec_il3_plusone_3"] = odeintr::laavec_il3_plusone_3;
out["lacvec_il1_3"] = odeintr::lacvec_il1_3;
out["lacvec_il1_plusone_3"] = odeintr::lacvec_il1_plusone_3;
out["lacvec_il3_3"] = odeintr::lacvec_il3_3;
out["lacvec_il3_plusone_3"] = odeintr::lacvec_il3_plusone_3;
out["lacvec_il4_plusone_3"] = odeintr::lacvec_il4_plusone_3;
out["muvec_il2_plusone_3"] = odeintr::muvec_il2_plusone_3;
out["muvec_il3_3"] = odeintr::muvec_il3_3;
out["muvec_il3_plusone_3"] = odeintr::muvec_il3_plusone_3;
out["nn_in1_3"] = odeintr::nn_in1_3;
out["nn_in2_3"] = odeintr::nn_in2_3;
out["nn_in3_plusone_3"] = odeintr::nn_in3_plusone_3;
out["gamvec_il3_4"] = odeintr::gamvec_il3_4;
out["laavec_il1_plusone_4"] = odeintr::laavec_il1_plusone_4;
out["laavec_il3_plusone_4"] = odeintr::laavec_il3_plusone_4;
out["lacvec_il1_4"] = odeintr::lacvec_il1_4;
out["lacvec_il1_plusone_4"] = odeintr::lacvec_il1_plusone_4;
out["lacvec_il3_4"] = odeintr::lacvec_il3_4;
out["lacvec_il3_plusone_4"] = odeintr::lacvec_il3_plusone_4;
out["lacvec_il4_plusone_4"] = odeintr::lacvec_il4_plusone_4;
out["muvec_il2_plusone_4"] = odeintr::muvec_il2_plusone_4;
out["muvec_il3_4"] = odeintr::muvec_il3_4;
out["muvec_il3_plusone_4"] = odeintr::muvec_il3_plusone_4;
out["nn_in1_4"] = odeintr::nn_in1_4;
out["nn_in2_4"] = odeintr::nn_in2_4;
out["nn_in3_plusone_4"] = odeintr::nn_in3_plusone_4;
  return out;
}
;




