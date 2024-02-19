#ifndef FitData_h
#define FitData_h

#include "TH1.h"

#include <vector>
namespace Fitters
{
class Data
{
private:
    std::vector<double> fX {};
    std::vector<double> fY {};
    double fXLow {};
    double fXUp {};
    double fBinWidth {};
    unsigned int fSize {};

public:
    Data(const TH1D& h, double xlow, double xup) : fXLow(xlow), fXUp(xup) { Fill(h); };

    // Getters
    double GetXLow() const { return fXLow; }
    double GetXUp() const { return fXUp; }
    double GetBinWidth() const { return fBinWidth; }
    double GetX(unsigned int i) { return fX[i]; }
    double GetY(unsigned int i) { return fY[i]; }
    unsigned int GetSize() const { return fSize; }
    double GetXLowEdge(unsigned int i) { return fX[i] - 0.5 * fBinWidth; }
    double GetXUpEdge(unsigned int i) { return fX[i] + 0.5 * fBinWidth; }
    int GetBin(double x) const;
    double Integral(double xmin, double xmax) const;

private:
    void Fill(const TH1D& h);
};
}; // namespace Fitters

#endif // !FitData_h
