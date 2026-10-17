// Native evaluator for the frozen empirical two-axis PCHIP cutoff surface.
// No distribution fitting, model-file loading or AP-MI estimator changes.
#include "apmi_null_model.h"
#include "apmi_pchip_data.h"

#include <cmath>
#include <string>

namespace
{
const double MINIMUM_P = 1e-7;
const double MAXIMUM_P = 0.05;
const double SCOPE_BOUNDARY_P = 0.0001;
const int MODERATE_KNOTS = 22;

int signOf(double value)
{
   return value > 0.0 ? 1 : value < 0.0 ? -1 : 0;
}

double endpointDerivative(double h0, double h1, double slope0, double slope1)
{
   // SciPy PchipInterpolator._edge_case, including its shape-preserving cap.
   double derivative = ((2.0 * h0 + h1) * slope0 - h0 * slope1) / (h0 + h1);
   if (signOf(derivative) != signOf(slope0))
      derivative = 0.0;
   else if (signOf(slope0) != signOf(slope1) &&
            std::fabs(derivative) > 3.0 * std::fabs(slope0))
      derivative = 3.0 * slope0;
   return derivative;
}

double pchipAt(const double* x, const double* y, int count, double query)
{
   if (count < 2 || !std::isfinite(query) || query < x[0] || query > x[count - 1])
      throw std::string("AP-MI PCHIP cannot extrapolate outside its probability knots.");
   double h[ApmiPchipData::NP];
   double slope[ApmiPchipData::NP];
   double derivative[ApmiPchipData::NP];
   for (int index = 0; index < count - 1; ++index)
   {
      h[index] = x[index + 1] - x[index];
      if (!(h[index] > 0.0) || !std::isfinite(h[index]))
         throw std::string("AP-MI PCHIP has invalid probability coordinates.");
      slope[index] = (y[index + 1] - y[index]) / h[index];
   }
   if (count == 2)
      derivative[0] = derivative[1] = slope[0];
   else
   {
      for (int index = 1; index < count - 1; ++index)
      {
         const double previous = slope[index - 1];
         const double next = slope[index];
         if (previous == 0.0 || next == 0.0 || signOf(previous) != signOf(next))
            derivative[index] = 0.0;
         else
         {
            const double weight1 = 2.0 * h[index] + h[index - 1];
            const double weight2 = h[index] + 2.0 * h[index - 1];
            const double weightedReciprocal =
               (weight1 / previous + weight2 / next) / (weight1 + weight2);
            derivative[index] = 1.0 / weightedReciprocal;
         }
      }
      derivative[0] = endpointDerivative(h[0], h[1], slope[0], slope[1]);
      derivative[count - 1] = endpointDerivative(h[count - 2], h[count - 3],
                                                 slope[count - 2], slope[count - 3]);
   }
   int interval = 0;
   while (interval < count - 2 && query >= x[interval + 1])
      ++interval;
   const double dx = query - x[interval];
   // SciPy CubicHermiteSpline's coefficient construction and Horner query.
   const double t = (derivative[interval] + derivative[interval + 1] -
                      2.0 * slope[interval]) / h[interval];
   const double cubic = t / h[interval];
   const double quadratic = (slope[interval] - derivative[interval]) / h[interval] - t;
   return ((cubic * dx + quadratic) * dx + derivative[interval]) * dx + y[interval];
}

int exactTrainingIndex(int m)
{
   for (int index = 0; index < ApmiPchipData::NM; ++index)
      if (ApmiPchipData::M[index] == m)
         return index;
   return -1;
}

int mInterval(double x)
{
   if (!std::isfinite(x) || x < ApmiPchipData::X[0] ||
       x > ApmiPchipData::X[ApmiPchipData::NM - 1])
      throw std::string("AP-MI PCHIP cannot extrapolate outside its sample-size knots.");
   int interval = 0;
   while (interval < ApmiPchipData::NM - 2 && x >= ApmiPchipData::X[interval + 1])
      ++interval;
   return interval;
}

double mLogCutoff(int probabilityIndex, int interval, double x)
{
   const double dx = x - ApmiPchipData::X[interval];
   return ((ApmiPchipData::COEFF[probabilityIndex][0][interval] * dx +
            ApmiPchipData::COEFF[probabilityIndex][1][interval]) * dx +
            ApmiPchipData::COEFF[probabilityIndex][2][interval]) * dx +
            ApmiPchipData::COEFF[probabilityIndex][3][interval];
}
}

namespace Apmi
{
double pchipCutoff(int m, double p, int nparLimit)
{
   if (m < 50 || m > 10000)
      throw std::string("AP-MI PCHIP requires an integer effective sample size m in [50,10000].");
   if (nparLimit != 40)
      throw std::string("AP-MI PCHIP is calibrated only for npar=40.");
   if (!std::isfinite(p) || p < MINIMUM_P || p > MAXIMUM_P)
      throw std::string("AP-MI PCHIP requires finite p in [1e-7,0.05].");

   // Descending P; p=0.0001 is a literal shared knot of separate scopes.
   const int first = p >= SCOPE_BOUNDARY_P ? 0 : MODERATE_KNOTS - 1;
   const int last = p >= SCOPE_BOUNDARY_P ? MODERATE_KNOTS - 1 : ApmiPchipData::NP - 1;
   const int count = last - first + 1;
   double predicted[ApmiPchipData::NP];
   double logPredicted[ApmiPchipData::NP];
   const int exact = exactTrainingIndex(m);
   const double logM = std::log(static_cast<double>(m));
   const double x = (logM - ApmiPchipData::LOG_M_MIN) / ApmiPchipData::LOG_M_WIDTH;
   const int interval = exact < 0 ? mInterval(x) : -1;
   for (int index = first; index <= last; ++index)
   {
      if (exact >= 0)
      {
         predicted[index] = ApmiPchipData::OBSERVED[index][exact];
         logPredicted[index] = ApmiPchipData::LOG_OBSERVED[index][exact];
      }
      else
      {
         predicted[index] = std::exp(mLogCutoff(index, interval, x));
         // Preserve the declared exp/log round trip between the two stages.
         logPredicted[index] = std::log(predicted[index]);
      }
      if (!std::isfinite(predicted[index]) || predicted[index] <= 0.0 ||
          predicted[index] > logM || !std::isfinite(logPredicted[index]))
         throw std::string("AP-MI PCHIP predicted an invalid or out-of-support cutoff.");
      if (index > first && predicted[index] < predicted[index - 1])
         throw std::string("AP-MI PCHIP probability knots cross; no repair is permitted.");
   }
   for (int index = first; index <= last; ++index)
      if (p == ApmiPchipData::P[index])
         return predicted[index];

   double probabilityCoordinates[ApmiPchipData::NP];
   double cutoffCoordinates[ApmiPchipData::NP];
   for (int index = 0; index < count; ++index)
   {
      const int source = last - index;
      probabilityCoordinates[index] = ApmiPchipData::LOG_P[source];
      cutoffCoordinates[index] = logPredicted[source];
   }
   const double cutoff = std::exp(pchipAt(probabilityCoordinates, cutoffCoordinates,
                                         count, std::log(p)));
   if (!std::isfinite(cutoff) || cutoff <= 0.0 || cutoff > logM)
      throw std::string("AP-MI PCHIP interpolated an invalid or out-of-support cutoff.");
   return cutoff;
}

const char* pchipPrimarySha256() { return ApmiPchipData::PRIMARY_SHA256; }
const char* pchipDataSha256() { return ApmiPchipData::DATA_SHA256; }
const char* pchipSummarySha256() { return ApmiPchipData::SUMMARY_SHA256; }
}
