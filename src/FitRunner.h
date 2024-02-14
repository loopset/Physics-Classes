#ifndef FitRunner_h
#define FitRunner_h

#include "TFitResult.h"

#include "Fit/BinData.h"
#include "Fit/Fitter.h"

#include "FitModel.h"

#include <string>
#include <unordered_map>
#include <vector>
namespace Fitters
{
class Runner
{
public:
    typedef std::vector<double> DoubleVec;
    typedef std::vector<bool> BoolVec;
    typedef std::vector<std::pair<double, double>> PairVec;
    typedef std::unordered_map<std::string, DoubleVec> Init;
    typedef std::unordered_map<std::string, PairVec> Bounds;
    typedef std::unordered_map<std::string, BoolVec> Fixed;
    typedef std::unordered_map<std::string, DoubleVec> Step;

private:
    ROOT::Fit::Fitter fFitter;
    Model* fModel {};

public:
    Runner() = default;
    Runner(Model* model) : fModel(model) {}
    // Setters
    void SetModelFunc(Model* model);
    void SetModelWrap(const ROOT::Math::WrappedMultiTF1& wrap);
    void SetInitial(const Init& pars);
    void SetBounds(const Bounds& bounds);
    void SetFixed(const Fixed& fixed);
    void SetStep(const Step& step);

    // Getters
    ROOT::Fit::Fitter& GetFitter() { return fFitter; }
    TFitResult GetFitResult() const { return TFitResult {fFitter.Result()}; }

    // Other methods
    bool Fit(const ROOT::Fit::BinData& data);

private:
    void ParametersAtLimit();
    bool CompareDoubles(double a, double b, double tol = 0.0001) const;
};
} // namespace Fitters

#endif
