#ifndef ActFitters_h
#define ActFitters_h

#include "Fitter.h"

#include "Fit/FitResult.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TH1.h"
#include "TString.h"
#include "TRegexp.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

Fitters::SpectrumData::SpectrumData(double min, double max, TH1* h)
    : fBinWidth(h->GetBinWidth(1))//assume bin with constant
{
    SetRangeFromHisto(min, max, h);
    FillVectors(h);
}

void Fitters::SpectrumData::SetRangeFromHisto(double min, double max, TH1* h)
{
    //This is to ensure a matching of the bins along all the code with the original TH1
    auto binlow {h->GetXaxis()->FindBin(min)};
    auto low {h->GetXaxis()->GetBinLowEdge(binlow)};
    auto binup {h->GetXaxis()->FindBin(max)};
    auto up {h->GetXaxis()->GetBinUpEdge(binup)};
    fRange = {low, up};
}

void Fitters::SpectrumData::FillVectors(TH1* h)
{
    for(int b = 1; b <= h->GetNbinsX(); b++)
    {
        auto x {h->GetBinCenter(b)};
        auto y {h->GetBinContent(b)};
        if(fRange.first <= x && x <= fRange.second)
        {
            fX.push_back(x);
            fY.push_back(y);
        }
    }
}

void Fitters::SpectrumData::FillPS(TH1* h)
{
    fPS.push_back({});
    for(int b = 1; b <= h->GetNbinsX(); b++)
    {
        auto x {h->GetBinCenter(b)};
        auto y {h->GetBinContent(b)};
        if(fRange.first <= x && x <= fRange.second)
        {
            fPS.back().push_back(y);
        }
    }
    //asset it has same size as main data
    if(fPS.back().size() != fX.size())
        throw std::runtime_error("Phase Space array size does not match main spectra one!");
}

int Fitters::SpectrumData::FindBin(double x) const
{
    int bin {};
    if(x <= fRange.first)
        return 0;
    else if (x >= fRange.second)
        return fX.size() - 1;
    else
    {
        bin = int(fX.size() * (x - fRange.first) / (fRange.second - fRange.first));
    }
    if(bin >= fX.size())
        std::cout<<"Greater than size"<<'\n';
    return bin;
}

double Fitters::SpectrumData::Integral(double xmin, double xmax) const
{
    double ret {};
    auto low {FindBin(xmin)};
    auto up {FindBin(xmax)};
    for(int bin = low; bin <= up; bin++)
    {
        //std::cout<<"Bin = "<<bin<<" y = "<<fY[bin]<<'\n';
        ret += fY[bin];
    }
    return ret;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Fitters::SpectrumFunction::SpectrumFunction(int ngaus, int nvoigt, const SpectrumData* data)
    : fNGauss(ngaus), fNVoigt(nvoigt), fNPS(data->GetNPS()), fData(data)
{
    CreateChart();
}

int Fitters::SpectrumFunction::GetNParFunc(const std::string &key)
{
    if(key == "g")
        return fNParGauss;
    else if(key == "v")
        return fNParVoigt;
    else if(key == "ps")
        return fNParPS;
    else if(key == "cte")
        return 1;
    else
        throw std::runtime_error("Received wrong key in GetNParFunc()");
}

int Fitters::SpectrumFunction::GetNFuncs() const
{
    return fNGauss + fNVoigt + fNPS + (int)fCte;
}

void Fitters::SpectrumFunction::CreateChart()
{
    fChart.clear();
    //Gaussians
    for(int g = 0; g < fNGauss; g++)
        for(int p = 0; p < fNParGauss; p++)
            fChart.push_back({"g", g});
    //Voigts
    for(int v = 0; v < fNVoigt; v++)
        for(int p = 0; p < fNParVoigt; p++)
            fChart.push_back({"v", v});
    //Phase spaces
    for(int ps = 0; ps < fNPS; ps++)
        for(int p = 0; p < fNParPS; p++)
            fChart.push_back({"ps", ps});
    if(fCte)
        fChart.push_back({"cte", 0});
}

Fitters::SpectrumFunction::ParPack Fitters::SpectrumFunction::UnpackCParameters(const double* pars) const
{
    ParGroup gaus {}; ParGroup voigt {}; ParGroup ps {};
    ParGroup cte {};
    for(int p = 0; p < GetNPars(); p++)
    {
        auto [type, idx] = fChart[p];
        if(type == "g")
            gaus[idx].push_back(pars[p]);
        else if(type == "v")
            voigt[idx].push_back(pars[p]);
        else if(type == "ps")
            ps[idx].push_back(pars[p]);
        else if(type == "cte")
            cte[idx].push_back(pars[p]);
        else
            throw std::runtime_error("Received wrong type of func in UnpackCParameters!");
    }
    return {gaus, voigt, ps, cte};
}


double Fitters::SpectrumFunction::EvalInBin(double bin, const double* pars, bool customx) const
{
    int niter {};
    double step {};
    double xStart {};
    if(!customx)//bin is in fact a bin
    {
        niter = fNDivBin;
        step = fData->GetBinWidth() / fNDivBin;
        xStart = fData->GetX()[bin] - 0.5 * fData->GetBinWidth();
    }
    else//bin is punctual
    {
        niter = 1;
        step = 0;
        xStart = bin;
    }
    //Unpack parameters = c-like array to our data structure
    auto pack = UnpackCParameters(pars);
    auto gaus = pack[0]; auto voigt = pack[1];
    auto phase = pack[2]; auto cte = pack[3];
    //Run for every iteration
    //1-->Gaussians
    std::vector<double> evalgauss(fNGauss);
    for(int g = 0; g < fNGauss; g++)
    {
        for(int i = 0; i < niter; i++)
        {
            auto xi {xStart + (i + 0.5) * step};
            evalgauss[g] += gaus[g][0] * TMath::Gaus(xi, gaus[g][1], gaus[g][2]);
        }
        evalgauss[g] /= niter;
    }
    //2-->Voigts
    std::vector<double> evalvoigt(fNVoigt);
    for(int v = 0; v < fNVoigt; v++)
    {
        for(int i = 0; i < niter; i++)
        {
            auto xi {xStart + (i + 0.5) * step};
            evalvoigt[v] += voigt[v][0] * TMath::Voigt(xi - voigt[v][1], voigt[v][2], voigt[v][3]);
        }
        evalvoigt[v] /= niter;
        //std::cout<<"voigt = "<<evalvoigt[v]<<'\n';
    }
    //3-->Phase spaces
    std::vector<double> evalps(fNPS);
    for(int ps = 0; ps < fNPS; ps++)
    {
        double val {};
        if(!customx)
            val = fData->GetPS(ps)[bin];
        else
            val = fData->GetPS(ps)[fData->FindBin(bin)];
        evalps[ps] = phase[ps].front() * val;
    }
    //Sum all contributions but cte
    double sum {};
    for(const auto& vec : {evalgauss, evalvoigt, evalps})
        for(const auto& e : vec)
            sum += e;
    //Sum cte contribution
    if(fCte)
        sum += cte[0].front();
    return sum;
}

double Fitters::SpectrumFunction::operator()(const double *pars) const
{
    // for(int i = 0; i < GetNPars(); i++)
    //     std::cout<<"Par "<<i<<" = "<<pars[i]<<'\n';
    double chi2 {};
    for(int bin = 0; bin < fData->GetNBinsX(); bin++)
    {
        double yexp {fData->GetY().at(bin)};
        double yfit {EvalInBin(bin, pars)};
        //std::cout<<"Exp = "<<yexp<<" fit = "<<yfit<<'\n';
        double diff {yexp - yfit};
        double sigma {EvalSigma(yexp, yfit)};
        chi2 += TMath::Power(diff / sigma, 2);
        if(!std::isfinite(chi2))
        {
            std::cout<<"Nan in bin = "<<bin<<" with content = "<<yexp<<'\n';
            std::cout<<"yfit = "<<yfit<<'\n';
            std::cout<<"Sigma = "<<sigma<<'\n';
            throw std::runtime_error("Nan");
        }
            
    }
    return chi2;
}

double Fitters::SpectrumFunction::EvalSigma(double nexp, double nfit) const
{
    double sigma {};
    if(nexp == 0)
        sigma = 1.84;
    if(nexp == 1 && nfit <= 1)
        sigma = 0.827;
    if(nexp == 1 && nfit > 1)
        sigma = 2.3;
    if(nexp == 2 && nfit <= 2)
        sigma = 1.292;
    if(nexp == 2 && nfit > 2)
        sigma = 2.64;
    if(nexp > 2 && (nfit > nexp))
        sigma = TMath::Sqrt(nexp) + 1;//upper limit. see eq.(13) in Some remarks on the error analysis in the case of poor statistics, by K.h. Schmidt
    if(nexp > 2 && (nfit <= nexp))
        sigma = TMath::Sqrt(nexp);//lower limit
    return sigma;
}

void Fitters::SpectrumFunction::PrintChart() const
{
    for(int idx = 0; idx < fChart.size(); idx++)
    {
        std::cout<<"Idx = "<<idx<<" par = "<<fChart[idx].first<<fChart[idx].second<<'\n';
    }
}

int Fitters::SpectrumFunction::GetNPars() const
{
    return fNGauss * fNParGauss + fNVoigt * fNParVoigt + fNPS * fNParPS + ((fCte) ? 1 : 0);
}

Fitters::SpectrumFitter::SpectrumFitter(Fitters::SpectrumFunction* func, const std::string& minimizer,
                                        const std::string& algorithm, int strategy)
    : fFunc(func)
{
    //Set fitter
    fFitter.Config().SetMinimizer(minimizer.c_str(), algorithm.c_str());
    //Options
    auto opts {fFitter.Config().MinimizerOptions()};
    //opts.SetMaxFunctionCalls(2);
    opts.SetPrintLevel(0);
    opts.SetStrategy(strategy);
    std::cout<<"== SpectrumFitter settings =="<<'\n';
    opts.Print();
    fFitter.Config().SetMinimizerOptions(opts);
}

bool Fitters::SpectrumFitter::Fit(bool minos)
{
    //Init parameters
    double pars[fFunc->GetNPars()];
    InitCParameters(pars);
    //Set fcn
    fFitter.SetFCN(fFunc->GetNPars(), *fFunc, pars, fFunc->GetDataSize(), true);//true if chi2; false otherwise
    //Config fit
    ConfigFit();
    fFitter.Config().SetUpdateAfterFit();
    //Run
    bool ok = fFitter.FitFCN();
    if(minos)
        fFitter.CalculateMinosErrors();
    fFitResult = fFitter.Result();
    fFitResult.Print(std::cout);
    PrintParametersAtLimit();
    return ok;
}

void Fitters::SpectrumFitter::InitCParameters(double* pars)
{
    auto chart {fFunc->GetChart()};
    int counter {};
    for(auto& [key, vec] : fInitPars)
    {
        //std::cout<<"Key = "<<key<<" vec size = "<<vec.size()<<'\n';
        auto [begin, end] = LocatePars(key);
        AssertDimensions(key, vec);
        int idx {};
        for(int p = begin; p <= end; p++)
        {
            //std::cout<<"Idx = "<<p<<" val = "<<vec[idx]<<'\n';
            pars[p] = vec[idx];
            idx++;
            counter++;
        }
    }
    //std::cout<<"Counter = "<<counter<<" size = "<<fFunc->GetNPars()<<'\n';
    if(counter != fFunc->GetNPars())
        throw std::runtime_error("Passed InitPars size does not match SpectrumFunction specs");
}

void Fitters::SpectrumFitter::ConfigFit()
{
    ConfigLabels();
    ConfigBounds();
    ConfigFixed();
    ConfigSteps();
}

void Fitters::SpectrumFitter::ConfigLabels()
{
    //Labels depends on fInitPars! But we ensured previously that these were correctly set (same dimensions as GetNPar)
    std::vector<std::string> labels {"_Amp", "_Mean", "_Sigma", "_Lg"};
    auto chart {fFunc->GetChart()};
    int counter {};
    for(auto& [key, _] : fInitPars)
    {
        auto [begin, end] = LocatePars(key);
        int idx {};
        for(int i = begin; i <= end; i++)
        {
            fFitter.Config().ParSettings(i)
                .SetName(key + labels[idx]);
            idx++;
            counter++;
        }
    }
}

void Fitters::SpectrumFitter::ConfigBounds()
{
    auto chart {fFunc->GetChart()};
    int counter {};
    for(auto& [key, vec] : fInitBounds)
    {
        auto [begin, end] = LocatePars(key);
        AssertDimensions(key, vec);
        int idx {};
        for(int i = begin; i <= end; i++)
        {
            auto [min, max] = vec[idx];
            if(min == -11 || max == -11)//disable bounds for that parameter
                continue;
            fFitter.Config().ParSettings(i)
                .SetLimits(min, max);
            idx++;
            counter++;
        }
    }
}

void Fitters::SpectrumFitter::ConfigFixed()
{
    auto chart {fFunc->GetChart()};
    int counter {};
    for(auto& [key, vec] : fFixedPars)
    {


        auto [begin, end] = LocatePars(key);
        AssertDimensions(key, vec);
        int idx {};
        for(int i = begin; i <= end; i++)
        {
            auto isFix {vec[idx]};
            if(isFix)
            {
                fFitter.Config().ParSettings(i)
                    .Fix();
            }
            idx++;
            counter++;
        }
    }
}

void Fitters::SpectrumFitter::ConfigSteps()
{
    auto chart {fFunc->GetChart()};
    int counter {};
    for(auto& [key, vec] : fStepPars)
    {
        auto [begin, end] = LocatePars(key);
        AssertDimensions(key, vec);
        int idx {};
        for(int i = begin; i <= end; i++)
        {
            auto step {vec[idx]};
            if(step == -11)//allow us to disable stepping
                continue;
            fFitter.Config().ParSettings(i)
                .SetStepSize(step);
            idx++;
            counter++;
        }
    }
}

std::pair<int, int> Fitters::SpectrumFitter::LocatePars(const std::string& key)
{
    auto chart {fFunc->GetChart()};
    std::string func {TString(key)(TRegexp("[a-z]+"))};
    int idx {std::stoi(TString(key)(TRegexp("[0-9]+")))};
    //std::cout<<"Locating key = "<<key<<" as func = "<<func<<" idx = "<<idx<<'\n';
    auto lambda = [&](const std::pair<std::string, int>& vals){return vals.first == func && vals.second == idx;};
    auto itfirst {std::find_if(chart.begin(), chart.end(), lambda)};
    auto first {std::distance(chart.begin(), itfirst)};
    auto itlast {std::find_if(chart.rbegin(), chart.rend(), lambda)};
    auto last {std::distance(chart.begin(), itlast.base()) - 1};//for last we have to decrease by one
    //std::cout<<" First = "<<first<<" last = "<<last<<'\n';
    if(itfirst == chart.end())
        throw std::runtime_error("Function " + key + " not found in Chart -> Update your SpectrumFunctions specs");
    return {first, last};
}

template<typename T>
void Fitters::SpectrumFitter::AssertDimensions(const std::string& key, const std::vector<T>& vec)
{
    std::string type {TString(key)(TRegexp("[a-z]+"))};
    auto length {vec.size()};
    //std::cout<<"(end - begin) = "<<length <<'\n'; 
    if(length != fFunc->GetNParFunc(type))
        throw std::runtime_error("Size of " + key + " pars does not match specifications");        
}

void Fitters::SpectrumFitter::PrintParametersAtLimit()
{
    auto res = fFitResult.Parameters();
    for(int i = 0; i < res.size(); i++)
    {
        double min {}; double max {};
        bool isBound {fFitResult.IsParameterBound(i)};
        fFitResult.ParameterBounds(i, min, max);
        if(isBound)
        {
            std::string name {fFitResult.GetParameterName(i)};
            //std::cout<<"Checking limits for "<<name<<'\n';
            //std::cout<<"Min = "<<min<<" Max = "<<max<<'\n';
            bool closeToMin {CompareDoubles(res[i], min)};
            bool closeToMax {CompareDoubles(res[i], max)};
            if(closeToMin)
            {
                std::cout<<"\033[1m\033[31m"<<"Parameter "<<name<<" reached LOWER limit of "<<res[i]<<"\033[0m"<<'\n';
            }
            if(closeToMax)
            {
                std::cout<<"\033[1m\033[31m"<<"Parameter "<<name<<" reached UPPER limit of "<<res[i]<<"\033[0m"<<'\n';
            }
        }        
    }
}

bool Fitters::SpectrumFitter::CompareDoubles(double a, double b, double tol)
{
    const auto greatedMagnitude {std::max(std::abs(a), std::abs(b))};
    bool comp {std::abs(a - b) < tol * greatedMagnitude};
    return comp;
}

void Fitters::SpectrumFitter::WriteToFile(const std::string &file)
{
    //Create map to save
    InitPars ret;
    //Types dont match: pack stores keys as ints (since the different funcs maps are vector elements)
    //But init pars needs everything in a sole map
    auto pack {fFunc->UnpackCParameters(fFitResult.GetParams())};
    auto gaus = pack[0]; auto voigt = pack[1];
    auto phase = pack[2]; auto cte = pack[3];
    //Gauss
    for(const auto& [i, pars] : gaus)
    {
        std::string key {"g" + std::to_string(i)};
        ret[key] = pars;
    }
    //Voigt
    for(const auto& [i, pars] : voigt)
    {
        std::string key {"v" + std::to_string(i)};
        ret[key] = pars;
    }
    //Phase space
    for(const auto& [i, pars] : phase)
    {
        std::string key {"ps" + std::to_string(i)};
        ret[key] = pars;
    }
    //Cte
    for(const auto& [i, pars] : cte)
    {
        std::string key {"cte" + std::to_string(i)};
        ret[key] = pars;
    }
    //Write
    auto* f {new TFile(file.c_str(), "recreate")};
    //Use TTree to avoid requiring dictionary for std collection
    auto* t {new TTree("InitPars", "FitResults in a fashion that is understandable by our classes")};
    std::string key {};
    t->Branch("key", &key);
    std::vector<double> values {};
    t->Branch("values", &values);
    for(auto& [k, v] : ret)
    {
        key = k;
        values = v;
        t->Fill();
    }
    t->Write();
    f->Close(); delete f;
}

void Fitters::SpectrumFitter::PrintShiftedMeans()
{
    //Values
    auto pars {std::vector<double>(fFitResult.GetParams(),
                                   fFitResult.GetParams() + fFunc->GetNPars())};
    //Only for gaussian
    std::cout<<"//////////////// Shifted gaussians in E_x //////////////////////"<<'\n';
    double shift {};
    for(int p = 0; p < pars.size(); p++)
    {
        TString label {fFitResult.GetParameterName(p)};
        if(label == "g0_Mean")
        {
            shift = pars[p];
            std::cout<<"-> Shift in E_x = "<<shift<<" MeV"<<'\n';
            continue;
        }
        if(label.Contains("_Mean"))
        {
            std::cout<<"-> Shifted "<<label<<" = "<<pars[p] - shift<<" MeV"<<'\n';
        }
    }
    std::cout<<"///////////////////////////////////////////////////////////////"<<'\n';
}

Fitters::SpectrumPlotter::SpectrumPlotter(Fitters::SpectrumData* data, Fitters::SpectrumFunction* func, ROOT::Fit::FitResult res)
    : fData(data), fFunc(func), fRes(res)
{}

TGraph* Fitters::SpectrumPlotter::GetGlobalFitGraph()
{
    auto* gr {new TGraph()};
    for(double x = fData->GetRange().first, xmax = fData->GetRange().second; x < xmax; x += fData->GetBinWidth() / 10)
    {
        auto y {fFunc->EvalAfterFitByX(x, fRes.GetParams())};
        //std::cout<<"x = "<<x<<" y = "<<y<<'\n';
        gr->SetPoint(gr->GetN(), x, y);
    }
    return gr;
}

TGraphErrors* Fitters::SpectrumPlotter::GetFitResiduals()
{
    auto* gr {new TGraphErrors()};
    gr->SetTitle(";E_{ex} [MeV];Residuals");
    for(int bin = 0; bin < fData->GetNBinsX(); bin++)
    {
        //original data
        auto [xexp, yexp] = fData->GetDataPair(bin);
        double sigma {TMath::Sqrt(yexp)};
        //fit
        auto [_, yfit] = fFunc->EvalAfterFit(bin, fRes.GetParams());
        double val {(yexp - yfit) / sigma};
        if(!std::isfinite(val))
        {
            sigma = 0;
            val = yexp - yfit;
        }
        gr->SetPoint(gr->GetN(), xexp, val);
        gr->SetPointError(gr->GetN() - 1, 0, sigma);
    }
    return gr;
}

std::unordered_map<std::string, TF1*> Fitters::SpectrumPlotter::GetIndividualFuncs()
{
    std::unordered_map<std::string, TF1*> ret;
    auto pack = fFunc->UnpackCParameters(fRes.GetParams());
    auto gaus = pack[0]; auto voigt = pack[1];
    auto phase = pack[2]; auto cte = pack[3];
    //Gauss
    for(const auto& [i, pars] : gaus)
    {
        std::string key {"g" + std::to_string(i)};
        ret[key] = new TF1(key.c_str(), "gaus",
                           fData->GetRange().first, fData->GetRange().second);
        ret[key]->SetParameters(&pars[0]);
    }
    //Voigt
    for(const auto& [i, pars] : voigt)
    {
        std::string key {"v" + std::to_string(i)};
        ret[key] = new TF1(key.c_str(), "[0] * TMath::Voigt(x - [1], [2], [3])",
                           fData->GetRange().first, fData->GetRange().second);
        ret[key]->SetParameters(&pars[0]);
    }
    //Cte
    for(const auto& [i, pars] : cte)
    {
        std::string key {"cte" + std::to_string(i)};
        ret[key] = new TF1(key.c_str(), "[0]",
                           fData->GetRange().first, fData->GetRange().second);
        ret[key]->SetParameters(&pars[0]);
    }
    return ret;
}

void Fitters::SpectrumPlotter::FillHistoFromFunc(TH1F* h, TF1* f)
{
    for(int bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        auto x {h->GetXaxis()->GetBinCenter(bin)};
        auto y {f->Eval(x)};
        h->SetBinContent(bin, y);
    }
}

std::unordered_map<std::string, TH1F*> Fitters::SpectrumPlotter::GetIndividualHists()
{
    std::unordered_map<std::string, TH1F*> ret;
    auto pack = fFunc->UnpackCParameters(fRes.GetParams());
    auto gaus = pack[0]; auto voigt = pack[1];
    auto phase = pack[2]; auto cte = pack[3];
    //Gauss
    for(const auto& [i, pars] : gaus)
    {
        std::string key {"g" + std::to_string(i)};
        //Function
        auto f = new TF1(key.c_str(), "gaus",
                           fData->GetRange().first, fData->GetRange().second);
        f->SetParameters(&pars[0]);
        //Histogram
        ret[key] = new TH1F(("h" + key).c_str(), key.c_str(), fData->GetNBinsX(), fData->GetRange().first, fData->GetRange().second);
        FillHistoFromFunc(ret[key], f);
        delete f;
    }
    //Voigt
    for(const auto& [i, pars] : voigt)
    {
        std::string key {"v" + std::to_string(i)};
        //Function
        auto f = new TF1(key.c_str(), "[0] * TMath::Voigt(x - [1], [2], [3])",
                           fData->GetRange().first, fData->GetRange().second);
        f->SetParameters(&pars[0]);
        //Histogram
        ret[key] = new TH1F(("h" + key).c_str(), key.c_str(), fData->GetNBinsX(), fData->GetRange().first, fData->GetRange().second);
        FillHistoFromFunc(ret[key], f);
        delete f;
    }
    //Cte
    for(const auto& [i, pars] : cte)
    {
        std::string key {"cte" + std::to_string(i)};
        //Function
        auto f = new TF1(key.c_str(), "[0]",
                           fData->GetRange().first, fData->GetRange().second);
        f->SetParameters(&pars[0]);
        //Histogram
        ret[key] = new TH1F(("h" + key).c_str(), key.c_str(), fData->GetNBinsX(), fData->GetRange().first, fData->GetRange().second);
        FillHistoFromFunc(ret[key], f);
        delete f;
    }
    // //Check
    // for(const auto& x : fData->GetX())
    // {
    //     std::cout<<"X in data = "<<x<<'\n';
    // }
    // for(int bin = 1; bin <= ret["g0"]->GetNbinsX(); bin++)
    // {
    //     std::cout<<"X in histos = "<<ret["g0"]->GetXaxis()->GetBinCenter(bin)<<'\n';
    // }
    return ret;
}

TH1F* Fitters::SpectrumPlotter::GetIndividualPS(int idx, const TString& hname)
{
    auto pack = fFunc->UnpackCParameters(fRes.GetParams());
    auto gaus = pack[0]; auto voigt = pack[1];
    auto phase = pack[2]; auto cte = pack[3];
    //Init histogram
    TString newname {};
    if(hname.Length() > 0)
        newname = hname;
    else
        newname = TString::Format("hPSFit%d", idx);
    auto* hret = new TH1F(newname,
                          TString::Format("PS number %d fitted", idx),
                          fData->GetNBinsX(), fData->GetRange().first, fData->GetRange().second);
    //Fill it
    for(int bin = 0; bin < fData->GetNBinsX(); bin++)
    {
        int hbin {hret->GetXaxis()->FindFixBin(fData->GetX(bin))};
        double counts {phase.at(idx).front() * fData->GetYPS(idx, bin)};
        double sigma {TMath::Sqrt(counts)};
        hret->SetBinContent(hbin, counts);
        hret->SetBinError(hbin, sigma);
    }
    return hret;
}
#endif
