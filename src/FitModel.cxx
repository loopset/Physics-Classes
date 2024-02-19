#include "FitModel.h"

#include "TMath.h"
#include "TRegexp.h"
#include "TSpline.h"

#include "PhysColors.h"

#include <ios>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

Fitters::Model::Model(int ngaus, int nvoigt, const std::vector<TH1D>& ps, bool withCte)
    : fNGauss(ngaus),
      fNVoigt(nvoigt),
      fNPS(ps.size()),
      fPS(ps),
      fCte(withCte)
{
    fPars = std::vector<double>(NPar());
    fParNames = std::vector<std::string>(NPar());
    fChart = std::vector<std::pair<std::string, unsigned int>>(NPar());
    InitFuncParNames();
    InitParNames();
    InitSplines();
}

void Fitters::Model::InitSplines()
{
    for(const auto& h : fPS)
        fSpePS.push_back(std::make_shared<TSpline3>(&h, "b2,e2", 0, 0));
}

double Fitters::Model::EvalPS(unsigned int i, double x) const
{
    if(fUseSpline)
        return fSpePS[i]->Eval(x);
    else
    {
        auto bin {fPS[i].GetXaxis()->FindBin(x)};
        return fPS[i].GetBinContent(bin);
    }
}

void Fitters::Model::InitFuncParNames()
{
    fFuncParNames = {
        {"g", {"_Amp", "_Mean", "_Sigma"}},
        {"v", {"_Amp", "_Mean", "_Sigma", "_Lg"}},
        {"ps", {"_Amp"}},
        {"cte", {"_Amp"}},
    };
}

void Fitters::Model::InitParNames()
{
    // Gaus
    for(unsigned int i = 0; i < fNGauss; i++)
    {
        for(unsigned int p = 0; p < fNParGauss; p++)
        {
            unsigned int idx {i * fNParGauss + p};
            // Set parameter name
            fParNames[idx] = "g" + std::to_string(i) + fFuncParNames["g"][p];
            // Push to chart
            fChart[idx] = {"g", i};
        }
    }
    unsigned int offset {static_cast<unsigned int>(fNGauss * fNParGauss)};
    // Voigt
    for(unsigned int i = 0; i < fNVoigt; i++)
    {
        for(unsigned int p = 0; p < fNParVoigt; p++)
        {
            unsigned int idx {i * fNParVoigt + p + offset};
            // Set name
            fParNames[idx] = "v" + std::to_string(i) + fFuncParNames["v"][p];
            // Push to chart
            fChart[idx] = {"v", i};
        }
    }
    offset += fNVoigt * fNParVoigt;
    // Phase spaces
    for(unsigned int i = 0; i < fNPS; i++)
    {
        for(unsigned int p = 0; p < fNParPS; p++)
        {
            unsigned int idx {i * fNParPS + p + offset};
            // Set name
            fParNames[idx] = "ps" + std::to_string(i) + fFuncParNames["ps"][p];
            // Chart
            fChart[idx] = {"ps", i};
        }
    }
    offset += fNPS * fNParPS;
    // Constant
    if(fCte)
    {
        // Set name
        fParNames[offset] = "cte0_Amp";
        // Push to chart
        fChart[offset] = {"cte", 0};
    }
}

std::pair<std::string, int> Fitters::Model::GetTypeIdx(const std::string& name) const
{

    std::string func {TString(name)(TRegexp("[a-z]+"))};
    int idx {std::stoi(TString(name)(TRegexp("[0-9]+")))};
    return {func, idx};
}

unsigned int Fitters::Model::GetIdxFromLabel(const std::string& typeIdx, unsigned int par) const
{
    // Attach label of parameter
    auto [name, idx] {GetTypeIdx(typeIdx)};
    auto label {typeIdx + fFuncParNames.at(name).at(par)};
    // And find index in vector
    for(unsigned int i = 0, size = fParNames.size(); i < size; i++)
    {
        if(fParNames[i] == label)
            return i;
    }
    throw std::runtime_error("Model::GetIdxFromLabel(): could not locate idx of label");
}

Fitters::Model::ParVec Fitters::Model::UnpackParameters(const double* pars) const
{
    ParPack gaus(fNGauss);
    ParPack voigt(fNVoigt);
    ParPack ps(fNPS);
    ParPack cte(fCte);
    for(int p = 0; p < NPar(); p++)
    {
        auto [type, idx] {fChart[p]};
        if(type == "g")
            gaus[idx].push_back(pars[p]);
        else if(type == "v")
            voigt[idx].push_back(pars[p]);
        else if(type == "ps")
            ps[idx].push_back(pars[p]);
        else if(type == "cte")
            cte[idx].push_back(pars[p]);
        else
            throw std::runtime_error("Model::UnpackParameters(): received wrong type of func");
    }
    return {std::move(gaus), std::move(voigt), std::move(ps), std::move(cte)};
}

unsigned int Fitters::Model::NPar() const
{
    return fNGauss * fNParGauss + fNVoigt * fNParVoigt + fNPS * fNParPS + static_cast<int>(fCte);
}

double Fitters::Model::EvalWithPacks(double x, ParPack& gaus, ParPack& voigt, ParPack& phase, ParPack& cte) const
{
    double ret {};
    // Run for every function
    // 1-->Gaussians
    for(int g = 0; g < fNGauss; g++)
        ret += gaus[g][0] * TMath::Gaus(x, gaus[g][1], gaus[g][2]);
    // 2-->Voigts
    for(int v = 0; v < fNVoigt; v++)
        ret += voigt[v][0] * TMath::Voigt(x - voigt[v][1], voigt[v][2], voigt[v][3]);
    // 3-->Phase spaces
    for(int ps = 0; ps < fNPS; ps++)
    {
        auto val {EvalPS(ps, x)};
        ret += phase[ps].front() * val;
    }
    // Sum cte contribution
    if(fCte)
        ret += cte[0].front();
    return ret;
}

double Fitters::Model::DoEvalPar(const double* xx, const double* p) const
{
    double x {xx[0]};
    // Unpack parameters = c-like array to our data structure
    // Check of null based on ROOT's forum answer
    // https://root-forum.cern.ch/t/error-on-doevalpar-with-custom-fit-model-and-root-fitter/58158/3?u=loopset
    const double* pars {(p) ? p : fPars.data()};
    auto pack = UnpackParameters(pars);
    auto gaus = pack[0];
    auto voigt = pack[1];
    auto phase = pack[2];
    auto cte = pack[3];
    // Call the public function
    return EvalWithPacks(x, gaus, voigt, phase, cte);
}

void Fitters::Model::Print() const
{
    std::cout << BOLDYELLOW << "···· Model func settings ····" << '\n';
    std::cout << "-> NGauss : " << fNGauss << '\n';
    std::cout << "-> NVoigt : " << fNVoigt << '\n';
    std::cout << "-> NPS    : " << fNPS << '\n';
    std::cout << "-> UseSpline ? " << std::boolalpha << fUseSpline << '\n';
    std::cout << "······························" << RESET << '\n';
}
