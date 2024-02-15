#ifndef FitObjective_h
#define FitObjective_h

#include "Math/IFunction.h"

#include "FitData.h"
#include "FitModel.h"

#include <memory>

namespace Fitters
{
class Objective : public ROOT::Math::IBaseFunctionMultiDimTempl<double>
{
private:
    // Pointer to model
    std::shared_ptr<Model> fModel;
    // Pointer to data
    std::shared_ptr<Data> fData;
    // Use integral on bin or not
    bool fUseDivisions {};
    // N of division on bin
    int fNdiv {20};
    // Use built-in ROOT integrator
    bool fUseIntegral {};

public:
    Objective() = default;
    Objective(const std::shared_ptr<Data>& data, const std::shared_ptr<Model>& model) : fData(data), fModel(model) {}
    Objective(const Data& data, const Model& model)
        : fData(std::make_shared<Data>(data)),
          fModel(std::shared_ptr<Model>(dynamic_cast<Model*>(model.Clone())))
    {
    }
    ~Objective() override = default;

    // Override methods
    unsigned int NDim() const override { return fModel->NPar(); }
    Objective* Clone() const override { return new Objective {fData, fModel}; }

    // Getters
    std::shared_ptr<Data> GetData() const { return fData; }
    std::shared_ptr<Model> GetModel() const { return fModel; }
    bool GetUseIntegral() const { return fUseIntegral; }
    bool GetUseDivisions() const { return fUseDivisions; }
    int GetNdiv() const { return fNdiv; }

    // Setters
    void SetUseIntegral(bool use) { fUseIntegral = use; }
    void SetUseDivisions(bool use) { fUseDivisions = use; }
    void SetNdiv(int div) { fNdiv = div; }

    // Others
    void Print() const;

private:
    double DoEval(const double* p) const override;
    double DoEvalWithIntegral(int i, const double* p) const;
    double DoEvalWithDivisions(double x, const double* p) const;
    double DoEvalSigma(double nexp, double nfit) const;
};
} // namespace Fitters

#endif // !FitObjective_h
