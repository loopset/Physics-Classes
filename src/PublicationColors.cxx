#ifndef PublicationColors_cxx
#define PublicationColors_cxx

#include "PublicationColors.h"
#include "TColor.h"
#include <stdexcept>
#include <string>
#include <vector>

PlotUtils::PublicationColors::PublicationColors(const std::string& colorset)
{
    DefineColors(colorset);
}

void PlotUtils::PublicationColors::DefineColors(const std::string& colorset)
{
    //CVD-friendly, defined in https://personal.sron.nl/~pault/#sec:qualitative
    
    std::vector<double> r, g, b;
    if(colorset == "bright")
    {
        r = {68, 102, 34, 204, 238, 170, 187};
        g = {119, 204, 136, 187, 102, 51, 187};
        b = {170, 238, 51, 68, 119, 119, 187};
    }
    else
        throw std::runtime_error("No colorset named " + colorset);
    //Define TColors
    for(int c = 0; c < r.size(); c++)
    {
        auto i {TColor::GetFreeColorIndex()};
        fColors.push_back(new TColor(i, r[c] / 256, g[c] / 256, b[c] / 256));
    }
}

int PlotUtils::PublicationColors::operator[](int i)
{
    return fColors.at(i)->GetNumber();
}

#endif
