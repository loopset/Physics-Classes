#ifndef FitCompositeTF1_h
#define FitCompositeTF1_h

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>
namespace Fitters
{
class CompositeTF1
{
public:
    using TF1Vec = std::vector<TF1*>;
    using ParVec = std::vector<std::vector<double>>;

private:
    std::vector<std::string> fKeys {};
    TF1Vec fFuncs {};
    TF1* fComposite {};
    // Store results of fit
    TGraphErrors* fGlobal {};
    std::vector<TH1D*> fHs {};
    // Some settings
    std::pair<double, double> fRange {};
    int fNPars {};

public:
    CompositeTF1() = default;
    CompositeTF1(std::initializer_list<TF1*> ptrs);

    // Setters
    TGraphErrors* GetGlobal() const { return fGlobal; }
    std::vector<TH1D*> GetHs() const { return fHs; }

    // Getters
    TF1* GetComposite() { return fComposite; }
    TF1* GetFunc(int idx) { return fFuncs[idx]; }

    // Others
    ParVec UnpackPars(double* p);
    double Eval(double x);
    double Eval(double* x, double* p);
    std::vector<double> Integral();
    void Draw(const std::string& opts = "same");
    void InitDraw(TH1* hmodel = nullptr); // create functions to be drawn assuming fComposite has been fitted


private:
    std::vector<double> PackPars();
    void Init();
};
} // namespace Fitters

#endif // FitCompositeTF1_h
