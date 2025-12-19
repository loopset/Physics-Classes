#ifndef FitCompositeTF1_h
#define FitCompositeTF1_h

#include "TF1.h"

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
    std::pair<double, double> fRange {};
    int fNPars {};

public:
    CompositeTF1() = default;
    CompositeTF1(std::initializer_list<TF1*> ptrs);

    // Setters
    // Getters
    TF1* GetComposite() { return fComposite; }
    TF1* GetFunc(int idx) { return fFuncs[idx]; }

    // Others
    ParVec UnpackPars(double* p);
    double Eval(double x);
    double Eval(double* x, double* p);
    std::vector<double> Integral();
    void Draw(const std::string& opts = "same");


private:
    std::vector<double> PackPars();
    void Init();
};
} // namespace Fitters

#endif // FitCompositeTF1_h
