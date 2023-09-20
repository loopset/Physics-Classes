#ifndef PublicationColors_h
#define PublicationColors_h

#include "TStyle.h"
#include "TColor.h"

#include <string>
#include <vector>

namespace PlotUtils
{
    class PublicationColors
    {
    private:
        std::vector<TColor*> fColors;
        
    public:
        PublicationColors(const std::string& colorset = "bright");
        int operator[](int i);
    private:
        void DefineColors(const std::string& colorset);
    };
}

#endif
