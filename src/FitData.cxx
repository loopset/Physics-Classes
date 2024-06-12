#include "FitData.h"

#include "PhysColors.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>

void Fitters::Data::Fill(const TH1D& h)
{
    // Before filling, check range
    auto hmin {h.GetXaxis()->GetXmin()};
    auto hmax {h.GetXaxis()->GetXmax()};
    if(fXLow < hmin || hmax < fXUp)
    {
        std::cout << BOLDRED << "Fitters::Data::Fill(): given range is wider than histogram one!" << '\n';
        std::cout << "  defaulting to hist range of: [" << hmin << ", " << hmax << "]" << RESET << '\n';
        fXLow = hmin;
        fXUp = hmax;
    }
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

void Fitters::Data::Print() const
{
    std::cout << BOLDGREEN << "---- Fitters::Data ----" << RESET << '\n';
    for(int i = 0; i < GetSize(); i++)
        std::cout << "X : " << fX[i] << " Y : " << fY[i] << '\n';
}
