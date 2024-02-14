#ifndef FitFCN_h
#define FitFCN_h

#include "Fit/BinData.h"
#include "Math/IFunction.h"

#include "FitModel.h"

#include <memory>

namespace Fitters
{
class ObjFCN : public ROOT::Math::IBaseFunctionMultiDimTempl<double>
{
private:
    // Pointer to model
    std::shared_ptr<Fitters::Model> fModel;
    // Pointer to data
    std::shared_ptr<ROOT::Fit::BinData> fData;

public:
    ObjFCN() = default;
    ObjFCN(const std::shared_ptr<ROOT::Fit::BinData>& data, const std::shared_ptr<Model>& model)
        : fData(data),
          fModel(model)
    {
    }
    ObjFCN(const ROOT::Fit::BinData& data, const Model& model)
        : fData(std::make_shared<ROOT::Fit::BinData>(data)),
          fModel(std::shared_ptr<Model>(dynamic_cast<Model*>(model.Clone())))
    {
    }
    ~ObjFCN() override = default;

    // Override methods
    unsigned int NDim() const override { return fModel->NPar(); }
    ObjFCN* Clone() const override { return new ObjFCN {fData, fModel}; }

private:
    double DoEval(const double* p) const override;
};
} // namespace Fitters

#endif // !FitFCN_h
