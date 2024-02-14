#ifndef FitRunner_cxx
#define FitRunner_cxx

#include "FitRunner.h"

#include "Fit/BinData.h"
#include "Math/WrappedMultiTF1.h"

#include "FitModel.h"

#include <iostream>

void Fitters::Runner::SetModelFunc(Fitters::Model* model)
{
    fFitter.SetFunction(*model, false);
}

void Fitters::Runner::SetModelWrap(const ROOT::Math::WrappedMultiTF1& wrap)
{
    fFitter.SetFunction(wrap, false);
}

void Fitters::Runner::SetInitial(const Init& pars)
{
    for(const auto& [func, vals] : pars)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            auto idx {fModel->GetIdxFromLabel(func, p)};
            fFitter.Config().ParSettings(idx).SetValue(vals[p]);
        }
    }
}

void Fitters::Runner::SetBounds(const Bounds& bounds)
{
    for(const auto& [func, vals] : bounds)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            auto idx {fModel->GetIdxFromLabel(func, p)};
            auto [min, max] {vals[p]};
            if(min == -11 || max == -11)
                continue;
            fFitter.Config().ParSettings(idx).SetLimits(min, max);
        }
    }
}

void Fitters::Runner::SetFixed(const Fixed& fixed)
{
    for(const auto& [func, vals] : fixed)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            auto idx {fModel->GetIdxFromLabel(func, p)};
            if(vals[p])
            {
                fFitter.Config().ParSettings(idx).Fix();
            }
        }
    }
}

void Fitters::Runner::SetStep(const Step& step)
{
    for(const auto& [func, vals] : step)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            auto idx {fModel->GetIdxFromLabel(func, p)};
            if(vals[p] == -11)
                continue;
            fFitter.Config().ParSettings(idx).SetStepSize(vals[p]);
        }
    }
}

bool Fitters::Runner::Fit(const ROOT::Fit::BinData& data)
{
    // Print settings
    fFitter.Config().MinimizerOptions().Print();
    // Perform fit
    auto ret {fFitter.Fit(data)};
    // Print
    fFitter.Result().Print(std::cout);
    // Check parameters at limit
    ParametersAtLimit();
    return ret;
}

bool Fitters::Runner::CompareDoubles(double a, double b, double tol) const
{
    const auto greatedMagnitude {std::max(std::abs(a), std::abs(b))};
    bool comp {std::abs(a - b) < tol * greatedMagnitude};
    return comp;
}


void Fitters::Runner::ParametersAtLimit()
{
    auto res {fFitter.Result()};
    for(int i = 0; i < res.Parameters().size(); i++)
    {
        double min {};
        double max {};
        bool isBound {res.IsParameterBound(i)};
        res.ParameterBounds(i, min, max);
        if(isBound)
        {
            std::string name {res.GetParameterName(i)};
            bool closeToMin {CompareDoubles(res.Parameter(i), min)};
            bool closeToMax {CompareDoubles(res.Parameter(i), max)};
            if(closeToMin)
            {
                std::cout << "\033[1m\033[31m"
                          << "Parameter " << name << " reached LOWER limit of " << res.Parameter(i) << "\033[0m"
                          << '\n';
            }
            if(closeToMax)
            {
                std::cout << "\033[1m\033[31m"
                          << "Parameter " << name << " reached UPPER limit of " << res.Parameter(i) << "\033[0m"
                          << '\n';
            }
        }
    }
}
#endif // !FitRunner_cxx
