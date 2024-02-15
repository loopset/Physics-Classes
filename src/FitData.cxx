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
