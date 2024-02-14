#ifndef FitFCN_cxx
#define FitFCN_cxx
#include "FitFCN.h"

#include "TMath.h"

double Fitters::ObjFCN::DoEval(const double* p) const
{
    // Do a chi2 fit
    double res {};
    for(int i = 0, size = fData->Size(); i < size; i++)
    {
        const auto x {fData->GetCoordComponent(i, 0)};
        const auto yexp {fData->Value(i)};
        const auto yfit {(*fModel)(x, p)};
        auto diff {yexp - yfit};
        auto sigma {(yexp < 1) ? 1 : TMath::Sqrt(yexp)};
        res += TMath::Power(diff / sigma, 2);
    }
    return res;
}

#endif // !FitFCN_cxx
