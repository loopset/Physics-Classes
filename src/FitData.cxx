#include "FitData.h"

#include <cstdlib>
#include <stdexcept>

void Fitters::Data::Fill(const TH1D& h)
{
    for(int bin = 1; bin <= h.GetNbinsX(); bin++)
    {
        auto x {h.GetXaxis()->GetBinCenter(bin)};
        auto y {h.GetBinContent(bin)};
        if(fXLow <= x && x <= fXUp)
        {
            fX.push_back(x);
            fY.push_back(y);
        }
    }
    // Set size
    fSize = fX.size();
    // Set bin width
    if(fSize < 2)
        throw std::runtime_error("Data::Fill(): size of data is < 2, a required minimum!");
    fBinWidth = std::abs(fX[1] - fX[0]);
}

int Fitters::Data::GetBin(double x) const
{
    int bin {};
    if(x <= fXLow)
        return 0;
    else if(x >= fXUp)
        return fX.size() - 1;
    else
        bin = int(fX.size() * (x - fXLow) / (fXUp - fXLow));
    if(bin < 0 || bin >= fX.size())
        throw std::runtime_error("Data::GetBin(): wrong bin calculation");
    return bin;
}

double Fitters::Data::Integral(double xmin, double xmax) const
{
    auto low {GetBin(xmin)};
    auto up {GetBin(xmax)};
    double ret {};
    for(auto bin = low; bin <= up; bin++)
        ret += fY.at(bin);
    return ret;
}
