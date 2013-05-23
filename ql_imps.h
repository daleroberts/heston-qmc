// 
// Dale Roberts <dale.o.roberts@gmail.com>
//
#ifndef volmodels_ql_imps_h
#define volmodels_ql_imps_h

#include <ql/quantlib.hpp>

using namespace QuantLib;

boost::shared_ptr<YieldTermStructure>
flatRate(const Date& today, const boost::shared_ptr<Quote>& forward, const DayCounter &dc);

boost::shared_ptr<YieldTermStructure>
flatRate(const boost::shared_ptr<Quote>& forward, const DayCounter &dc);

boost::shared_ptr<YieldTermStructure>
flatRate(Rate forward, const DayCounter &dc);

boost::shared_ptr<QuantLib::BlackVolTermStructure>
flatVol(const QuantLib::Date& today, const boost::shared_ptr<Quote>& vol, const QuantLib::DayCounter& dc);

boost::shared_ptr<QuantLib::BlackVolTermStructure>
flatVol(const QuantLib::Date& today, Real vol, const DayCounter &dc);

Real HestonCallAnalytic(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K);

Real HestonCallMC(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K, int discretization, int trials, int steps, int seed);

Real BlackScholesCall(Real r, Real sigma, Real S0, Real T, Real K);

Real HestonFdBarrierDownOutCall(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K, Real H);

#endif
