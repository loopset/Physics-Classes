#include "FitRunner.h"

#include "TFile.h"
#include "TFitResult.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

void Fitters::Runner::SetFCN()
{
    // Toy double pars[NDim] that just serve as initialization
    std::vector<double> pars(fObj.NDim(), 0);
    // Objective func is managed by us
    // But model func is cloned when set
    fFitter.SetFCN(fObj, *(fObj.GetModel()), &(pars.front()), fObj.GetData()->GetSize(), true);
}

void Fitters::Runner::SetInitial(const Init& pars)
{
    for(const auto& [func, vals] : pars)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            auto idx {fObj.GetModel()->GetIdxFromLabel(func, p)};
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
            auto idx {fObj.GetModel()->GetIdxFromLabel(func, p)};
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
            auto idx {fObj.GetModel()->GetIdxFromLabel(func, p)};
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
            auto idx {fObj.GetModel()->GetIdxFromLabel(func, p)};
            if(vals[p] == -11)
                continue;
            fFitter.Config().ParSettings(idx).SetStepSize(vals[p]);
        }
    }
}

bool Fitters::Runner::Fit()
{
    // Print settings
    fFitter.Config().MinimizerOptions().Print();
    fObj.GetModel()->Print();
    fObj.Print();
    // Perform fit
    auto ret {fFitter.FitFCN()};
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

void Fitters::Runner::Write(const std::string& file) const
{
    // open file
    auto f {std::make_unique<TFile>(file.c_str(), "recreate")};
    // save parameter names (not saved by TFitResult :/)
    std::vector<std::string> names;
    for(int i = 0; i < fObj.NDim(); i++)
        names.push_back(fObj.GetModel()->ParameterName(i));
    f->WriteObject(&names, "ParNames");
    // and now fit result
    TFitResult res {fFitter.Result()};
    f->WriteObject(&res, "FitResult");
}
