// // 
// Dale Roberts <dale.o.roberts@gmail.com>
//
// Call QuantLib implementations

#include <stdio.h>
#include <ql/quantlib.hpp>
#include "ql_imps.h"

using namespace QuantLib;

boost::shared_ptr<YieldTermStructure>
flatRate(const Date& today, const boost::shared_ptr<Quote>& forward, const DayCounter &dc) {
    return boost::shared_ptr<YieldTermStructure>(new FlatForward(today, Handle<Quote>(forward),dc));
}

boost::shared_ptr<YieldTermStructure>
flatRate(const boost::shared_ptr<Quote>& forward, const DayCounter &dc) {
    return boost::shared_ptr<YieldTermStructure>(new FlatForward(0, NullCalendar(), Handle<Quote>(forward), dc));
}

boost::shared_ptr<YieldTermStructure>
flatRate(Rate forward, const DayCounter &dc) {
    return flatRate(boost::shared_ptr<Quote>(new SimpleQuote(forward)), dc);
}

boost::shared_ptr<QuantLib::BlackVolTermStructure>
flatVol(const QuantLib::Date& today, const boost::shared_ptr<Quote>& vol, const QuantLib::DayCounter& dc) {
    return boost::shared_ptr<QuantLib::BlackVolTermStructure>(new QuantLib::BlackConstantVol(today, NullCalendar(), Handle<Quote>(vol), dc));
}

boost::shared_ptr<QuantLib::BlackVolTermStructure>
flatVol(const QuantLib::Date& today, Real vol, const DayCounter &dc) {
    return flatVol(today, boost::shared_ptr<Quote>(new SimpleQuote(vol)), dc);
}

Real HestonCallAnalytic(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K) {
    
    // Square root process
    //  dv(t) = \kappa (\theta - v(t))dt + \sigma \sqrt{v(t)}dW(t)
    // v(0) = v0.
    
    Calendar calendar = NullCalendar();
    DayCounter dayCounter = SimpleDayCounter();
    
    Date todaysDate = calendar.adjust(Date::todaysDate());
    Settings::instance().evaluationDate() = todaysDate;

    Date maturityDate = calendar.advance(todaysDate, T*Years);
    printf("\nT:%f\n", dayCounter.yearFraction(todaysDate, maturityDate));
                
    Option::Type type(Option::Call);
    Rate riskFreeRate  = r;
    
    boost::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type, K));
    boost::shared_ptr<Exercise> europeanExercise(new EuropeanExercise(maturityDate));
    VanillaOption option(payoff, europeanExercise);
    
    Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
    
    Handle<YieldTermStructure> riskFreeTS(flatRate(riskFreeRate, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0, dayCounter));
    
    boost::shared_ptr<HestonProcess> process(new HestonProcess(riskFreeTS, dividendTS, underlyingH, v0, kappa, theta, sigma, rho));
    
    boost::shared_ptr<AnalyticHestonEngine> engine(new AnalyticHestonEngine(boost::shared_ptr<HestonModel>(new HestonModel(process)),64));
    
    option.setPricingEngine(engine);
        
    return (Real) option.NPV();
}

Real HestonCallMC(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K, int discretization, int trials, int steps, int seed) {

    HestonProcess::Discretization d;
    switch (discretization) 
    {
        case 1:
        d = HestonProcess::PartialTruncation;
        break;

        case 2:
        d = HestonProcess::FullTruncation;
        break;

        case 3:
        d = HestonProcess::Reflection;
        break;

        case 4:
        d = HestonProcess::NonCentralChiSquareVariance;
        break;

        case 5:
        d = HestonProcess::QuadraticExponential;
        break;

        case 6:
        d = HestonProcess::QuadraticExponentialMartingale;
        break;
    }

    Calendar calendar = NullCalendar();
    DayCounter dayCounter = SimpleDayCounter();
    
    Date todaysDate = calendar.adjust(Date::todaysDate());
    Settings::instance().evaluationDate() = todaysDate;

    Date maturityDate = calendar.advance(todaysDate, T*Years);
        
    Option::Type type(Option::Call);
    Rate riskFreeRate  = r;
    
    boost::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type, K));
    boost::shared_ptr<Exercise> europeanExercise(new EuropeanExercise(maturityDate));
    VanillaOption option(payoff, europeanExercise);
    
    Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
    
    Handle<YieldTermStructure> riskFreeTS(flatRate(riskFreeRate, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0, dayCounter));
    
    boost::shared_ptr<HestonProcess> process(new HestonProcess(riskFreeTS, dividendTS, underlyingH, v0, kappa, theta, sigma, rho, d));

//    boost::shared_ptr<PricingEngine> engine =
//    MakeMCEuropeanHestonEngine<PseudoRandom>(process)
//    .withSteps(10)
//    .withAntitheticVariate()
//    .withSamples(5000)
//    .withSeed(1234);

    boost::shared_ptr<PricingEngine> engine =
        MakeMCEuropeanHestonEngine<PseudoRandom>(process)
        .withSteps(steps)
        .withSamples(trials)
        .withSeed(seed);
    
    option.setPricingEngine(engine);
        
    return (Real) option.NPV();
}

Real BlackScholesCall(Real r, Real sigma, Real S0, Real T, Real K) {
    Calendar calendar = NullCalendar();
    DayCounter dayCounter = SimpleDayCounter();
    
    Date todaysDate = calendar.adjust(Date::todaysDate());
    Date maturityDate = calendar.advance(todaysDate, T*Years);

    Settings::instance().evaluationDate() = todaysDate;
    Option::Type type(Option::Call);
    
    boost::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type, K));
    boost::shared_ptr<Exercise> europeanExercise(new EuropeanExercise(maturityDate));
    VanillaOption option(payoff, europeanExercise);
    
    Handle<Quote> S0H(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
    
    Handle<BlackVolTermStructure> volatilityTS(flatVol(todaysDate, sigma, dayCounter));
    Handle<YieldTermStructure> riskFreeTS(flatRate(r, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0, dayCounter));
    
    boost::shared_ptr<BlackScholesMertonProcess> process(new BlackScholesMertonProcess(S0H, dividendTS, riskFreeTS, volatilityTS));
    boost::shared_ptr<PricingEngine> engine(new AnalyticEuropeanEngine(process));
    
    option.setPricingEngine(engine);

    return (Real) option.NPV();
}

Real HestonFdBarrierDownOutCall(Real r, Real rho, Real kappa, Real theta, Real sigma, Real v0, Real S0, Real T, Real K, Real H) 
{
    Calendar calendar = NullCalendar();
    DayCounter dayCounter = SimpleDayCounter();
    
    Date todaysDate = calendar.adjust(Date::todaysDate());
    Settings::instance().evaluationDate() = todaysDate;

    Date maturityDate = calendar.advance(todaysDate, T*Years);
                
    Option::Type type(Option::Call);
    Rate riskFreeRate  = r;
    
    boost::shared_ptr<StrikedTypePayoff> payoff(new PlainVanillaPayoff(type, K));
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(maturityDate));
    
    Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
    
    Handle<YieldTermStructure> riskFreeTS(flatRate(riskFreeRate, dayCounter));
    Handle<YieldTermStructure> dividendTS(flatRate(0, dayCounter));
    
    boost::shared_ptr<HestonProcess> process(new HestonProcess(riskFreeTS, dividendTS, underlyingH, v0, kappa, theta, sigma, rho));
    
    boost::shared_ptr<PricingEngine> engine(new FdHestonBarrierEngine(
                    boost::shared_ptr<HestonModel>(new HestonModel(process)),
                    200,400,100));
    
    BarrierOption option(Barrier::DownOut, H, 0.0, payoff, exercise);
    option.setPricingEngine(engine);
        
    return (Real) option.NPV();
}
