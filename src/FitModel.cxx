#ifndef FitModel_cxx
#define FitModel_cxx

#include "FitModel.h"

#include "TMath.h"
#include "TRegexp.h"

#include <string>
#include <vector>

Fitters::Model::Model(int ngaus, int nvoigt, const std::vector<TH1*>& ps, bool withCte)
    : fNGauss(ngaus),
      fNVoigt(nvoigt),
      fNPS(ps.size()),
      fPS(ps),
      fCte(withCte)
{
    fPars = std::vector<double>(NPar());
    fParNames = std::vector<std::string>(NPar());
    InitParNames();
}

double Fitters::Model::EvalPS(unsigned int i, double x) const
{
    auto* h {fPS[i]};
    auto bin {h->FindBin(x)};
    return h->GetBinContent(bin);
}

void Fitters::Model::InitParNames()
{
    // Gaus
    std::vector<std::string> gl {"_Amp", "_Mean", "_Sigma"};
    for(unsigned int i = 0; i < fNGauss; i++)
    {
        for(unsigned int p = 0; p < fNParGauss; p++)
        {
            unsigned int idx {i * fNParGauss + p};
            fParNames[idx] = "g" + std::to_string(i) + gl[p];
        }
    }
    unsigned int offset {static_cast<unsigned int>(fNGauss * fNParGauss)};
    // Voigt
    std::vector<std::string> vl {"_Amp", "_Mean", "_Sigma", "_Lg"};
    for(unsigned int i = 0; i < fNVoigt; i++)
    {
        for(unsigned int p = 0; p < fNParVoigt; p++)
        {
            unsigned int idx {i * fNParVoigt + p + offset};
            fParNames[idx] = "v" + std::to_string(i) + vl[p];
        }
    }
    offset += fNVoigt * fNParVoigt;
    // Phase spaces
    std::vector<std::string> pl {"_Amp"};
    for(unsigned int i = 0; i < fNPS; i++)
    {
        for(unsigned int p = 0; p < fNParPS; p++)
        {
            unsigned int idx {i * fNParPS + p + offset};
            fParNames[idx] = "ps" + std::to_string(i) + pl[p];
        }
    }
    offset += fNPS * fNParPS;
    // Constant
    if(fCte)
        fParNames[offset] = "cte0_Amp";
}

std::pair<std::string, int> Fitters::Model::GetTypeIdx(const std::string& name) const
{

    std::string func {TString(name)(TRegexp("[a-z]+"))};
    int idx {std::stoi(TString(name)(TRegexp("[0-9]+")))};
    return {func, idx};
}

Fitters::Model::ParPack Fitters::Model::UnpackParameters(const double* pars) const
{
    ParGroup gaus {};
    ParGroup voigt {};
    ParGroup ps {};
    ParGroup cte {};
    for(int p = 0; p < NPar(); p++)
    {
        std::cout << "i : " << p << '\n';
        auto [type, idx] {GetTypeIdx(ParameterName(p))};
        std::cout << "Type : " << type << " idx : " << idx << '\n';
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
    return {gaus, voigt, ps, cte};
}

unsigned int Fitters::Model::NPar() const
{
    return fNGauss * fNParGauss + fNVoigt * fNParVoigt + fNPS * fNParPS + static_cast<int>(fCte);
}


double Fitters::Model::DoEvalPar(const double* xx, const double* p) const
{
    double x {xx[0]};
    // Unpack parameters = c-like array to our data structure
    // auto pack = UnpackParameters(p);
    if(!p)
        std::cout << "No p" << '\n';
    // auto gaus = pack[0];
    // auto voigt = pack[1];
    // auto phase = pack[2];
    // auto cte = pack[3];
    double sum {0};
    // Run for every iteration
    // 1-->Gaussians
    // std::vector<double> evalgauss(fNGauss);
    // for(int g = 0; g < fNGauss; g++)
    // {
    //     // evalgauss[g] += gaus[g][0] * TMath::Gaus(x, gaus[g][1], gaus[g][2]);
    // }
    // // 2-->Voigts
    // std::vector<double> evalvoigt(fNVoigt);
    // for(int v = 0; v < fNVoigt; v++)
    // {
    //     // evalvoigt[v] += voigt[v][0] * TMath::Voigt(x - voigt[v][1], voigt[v][2], voigt[v][3]);
    // }
    // // 3-->Phase spaces
    // std::vector<double> evalps(fNPS);
    // for(int ps = 0; ps < fNPS; ps++)
    // {
    //     // auto val {EvalPS(ps, x)};
    //     // evalps[ps] = phase[ps].front() * val;
    // }
    // // Sum all contributions but cte
    // double sum {};
    // for(auto& vec : {&evalgauss, &evalvoigt, &evalps})
    //     for(const auto& e : *vec)
    //         sum += e;
    // // Sum cte contribution
    // if(fCte)
    //     sum += cte[0].front();
    return sum;
}

#endif // !FitModel_cxx
