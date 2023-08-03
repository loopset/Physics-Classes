#ifndef AngularDistribution_cxx
#define AngularDistribution_cxx

#include "AngularDistribution.h"

#include "Fitter.h"
#include "Rtypes.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TString.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

void AngularDistribution::AngularFitter::SetData(double xmin, double xmax, std::vector<TH1F*> hEx)
{
    for(auto& h : hEx)
    {
        std::cout<<"AngularDistribution: Adding data from histo titled "<<h->GetTitle()<<'\n';
        fData.push_back(Fitters::SpectrumData(xmin, xmax, h));
        fHistos.push_back(h);
    }
}

void AngularDistribution::AngularFitter::ReadFromFile(const std::string &file, const std::string& treename)
{
    auto* f {new TFile(file.c_str())};
    auto* t {f->Get<TTree>(treename.c_str())};
    if(!t)
        throw std::runtime_error("Error! Could not load TTree with InitPars config for AngularFitter");
    std::string* key {}; std::vector<double>* values {};
    t->SetBranchAddress("key", &key);
    t->SetBranchAddress("values", &values);
    std::cout<<"========================================"<<'\n';
    std::cout<<"Reading InitPars from file "<<file<<'\n';
    //Read number of functions
    int ngaus {}; int nvoigt {};
    int nps {}; int nctes {};
    for(int i = 0; i < t->GetEntries(); i++)
    {
        t->GetEntry(i);
        fInitPars[*key] = *values;
        //debug
        std::cout<<"Key = "<<*key<<'\n';
        for(auto& val : *values)
            std::cout<<"  ->Val = "<<val<<'\n';
        //Read number of functions
        TString k {*key};
        if(k.Contains("g"))
            ngaus++;
        else if(k.Contains("v"))
            nvoigt++;
        else if(k.Contains("ps"))
            nps++;
        else if(k.Contains("cte"))
            nctes++;
        else
            throw std::runtime_error("Received wrong key when reading file -> Check saved keys!");
    }
    std::cout<<"----------------------------------------"<<'\n';
    std::cout<<"Setting function model with parameters..."<<'\n';
    std::cout<<"NGauss = "<<ngaus<<'\n';
    std::cout<<"NVoigt = "<<nvoigt<<'\n';
    std::cout<<"NPS    = "<<nps<<'\n';
    std::cout<<"NCtes  = "<<nctes<<'\n';
    std::cout<<"========================================"<<'\n';
    f->Close(); delete f;
    SetFuncModel(ngaus, nvoigt, nctes);
}

void AngularDistribution::AngularFitter::AddPhaseSpace(TH1 *hPS)
{
    std::cout<<"AngularDistribution: Adding PS titled "<<hPS->GetTitle()<<'\n';
    for(auto& data : fData)
        data.AddPhaseSpace(hPS);
}

void AngularDistribution::AngularFitter::SetFuncModel(int ngaus, int nvoigt, bool cte)
{
    for(auto& data : fData)
    {
        fFuncs.push_back(Fitters::SpectrumFunction(ngaus, nvoigt, &data));
        if(cte)
            fFuncs.back().EnableCteBackground();
    }
}

void AngularDistribution::AngularFitter::Run()
{
    for(int f = 0; f < fFuncs.size(); f++)
    {
        std::cout<<"******* AngularDistribution "<<f<<" ********"<<'\n';
        Fitters::SpectrumFitter fitter {&fFuncs[f]};
        fitter.SetInitPars(fInitPars);
        fitter.SetInitBounds(GetInitBounds());
        fitter.SetFixedPars(GetFixedPars());
        fitter.Fit();
        fRes.push_back(fitter.GetFitResult());
        std::cout<<"*********************************************"<<'\n';
    }
}

Fitters::SpectrumFitter::InitBounds AngularDistribution::AngularFitter::GetInitBounds()
{
    //Bounds are only set for amplitudes
    Fitters::SpectrumFitter::InitBounds bounds {};
    for(const auto& [key, vals] : fInitPars)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            std::pair<double, double> pair;
            if(p == 0)
                pair = {0, 1000};
            else
                pair = {-11, -11};

            bounds[key].push_back(pair);
        }
    }
    return bounds;
}

Fitters::SpectrumFitter::FixedPars AngularDistribution::AngularFitter::GetFixedPars()
{
    //Fix all but amplitudes
    Fitters::SpectrumFitter::FixedPars fixed {};
    for(const auto& [key, vals] : fInitPars)
    {
        for(int p = 0; p < vals.size(); p++)
        {
            bool boo {true};
            if(p == 0)
                boo = false;
            fixed[key].push_back(boo);
        }
    }
    return fixed;
}

TCanvas* AngularDistribution::AngularFitter::Print(double xmin, double xmax)
{
    auto* cres {new TCanvas("cres", "Angular distributions")};
    cres->DivideSquare(fData.size());
    for(int i = 0; i < fData.size(); i++)
    {
        //Init painter
        Fitters::SpectrumPlotter plotter {&fData[i], &fFuncs[i], fRes[i]};
        auto* gfit {plotter.GetGlobalFitGraph()};
        auto fits {plotter.GetIndividualFuncs()};
        TH1F* psfit {};
        if(fData[i].GetNPS())
        {
            psfit = plotter.GetIndividualPS(0, TString::Format("hPSFitI%dCM%d", 0, i));//just 0 for now
            // style for PS
            psfit->SetLineColor(kOrange);
        }
        //Set styles
        // Global fit
        gfit->SetLineColor(kRed); gfit->SetLineWidth(2);
        
        //Draw!
        cres->cd(i + 1);
        //Set range
        if(xmin != 0 && xmax != 0)
            fHistos[i]->GetXaxis()->SetRangeUser(xmin, xmax);
        fHistos[i]->Draw();
        gfit->Draw("same");
        int idx {};
        for(auto& [_, g] : fits)
        {
            int color {idx + 6};
            if(color == 10)
                color = 46;//skip white color
            g->SetLineWidth(2);
            g->SetLineColor(color);
            g->Draw("same");
            idx++;
        }
        if(psfit)
            psfit->Draw("hist same");
    }
    return cres;
}

void AngularDistribution::AngularFitter::ComputeFuncIntegral(int idx)
{
    //Initialize row of vector
    fIntegrals.push_back({});
    //Unpack parameters
    auto pack {fFuncs[idx].UnpackCParameters(fRes[idx].GetParams())};
    auto gaus = pack[0]; auto voigt = pack[1];
    auto xmin {fData[idx].GetRange().first}; auto xmax {fData[idx].GetRange().second};
    //We are only interested in gaus and voigt fit parameters
     //Gauss
    for(const auto& [i, pars] : gaus)
    {
        std::string key {"g" + std::to_string(i)};
        auto f = new TF1(key.c_str(), "gaus", xmin, xmax);
        f->SetParameters(&pars[0]);
        //Integrate
        fIntegrals.back()[key] = f->Integral(xmin, xmax) / fData[idx].GetBinWidth();
        delete f;
    }
    //Voigt
    for(const auto& [i, pars] : voigt)
    {
        std::string key {"v" + std::to_string(i)};
        auto f = new TF1(key.c_str(), "[0] * TMath::Voigt(x - [1], [2], [3])",
                         xmin, xmax);
        f->SetParameters(&pars[0]);
        //Integrate
        fIntegrals.back()[key] = f->Integral(xmin, xmax) / fData[idx].GetBinWidth();
        delete f;
    }
}

void AngularDistribution::AngularFitter::ComputeCountSum(int idx, int nsigma)
{
    //Initialize row of vector
    fIntegrals.push_back({});
    //Unpack parameters
    auto pack {fFuncs[idx].UnpackCParameters(fRes[idx].GetParams())};
    auto gaus = pack[0]; auto voigt = pack[1];
    auto xmin {fData[idx].GetRange().first}; auto xmax {fData[idx].GetRange().second};
    //Direct integral in data
    //Gauss
    for(const auto& [i, pars] : gaus)
    {
        std::string key {"g" + std::to_string(i)};
        //Set integral range
        double mean {pars[1]}; double sigma {pars[2]};
        double low {mean - nsigma * sigma}; double up {mean + nsigma * sigma};
        auto integral {fData[idx].Integral(low, up)};
        //Scale corresponding to sigma value
        double scale {};
        if(nsigma == 1)
            scale = 0.68;
        else if(nsigma == 2)
            scale = 0.95;
        else
            throw std::runtime_error("No nsigma correction factor implemented");
        fIntegrals.back()[key] = integral / scale;
        //Print
        std::cout<<"===== Gauss sum number "<<idx<<" ====="<<'\n';
        std::cout<<"Low = "<<low<<" Up = "<<up<<'\n';
        std::cout<<"Integral (not scaled) = "<<integral<<'\n';
        std::cout<<"=============================="<<'\n';
    }
    //Voigt
    for(const auto& [i, pars] : voigt)
    {
        std::string key {"v" + std::to_string(i)};
        //Set integral range
        double mean {pars[1]}; double sigma {pars[2]};
        double low {mean - nsigma * sigma}; double up {mean + nsigma * sigma};
        auto integral {fData[idx].Integral(low, up)};
        //Scale corresponding to sigma value
        double scale {};
        if(nsigma == 1)
            scale = 0.68;
        else if(nsigma == 2)
            scale = 0.95;
        else
            throw std::runtime_error("No nsigma correction factor implemented");
        fIntegrals.back()[key] = integral / scale;
        //Print
        std::cout<<"===== Voigt sum number "<<idx<<" ====="<<'\n';
        std::cout<<"Low = "<<low<<" Up = "<<up<<'\n';
        std::cout<<"Integral (not scaled) = "<<integral<<'\n';
        std::cout<<"=============================="<<'\n';
    }
}

void AngularDistribution::AngularFitter::ComputeYield(const std::string& mode, int nsigma)
{
    fIntegrals.clear();
    for(int f = 0; f < fFuncs.size(); f++)
    {
        if(mode == "integral")
            ComputeFuncIntegral(f);
        else if(mode == "sum")
            ComputeCountSum(f, nsigma);
        else
            throw std::runtime_error("Error: modes available are integral and sum");
    }
}

std::vector<double> AngularDistribution::AngularFitter::GetIntegralsForPeak(const std::string &key)
{
    std::vector<double> ret;
    for(const auto& map : fIntegrals)
    {
        for(const auto& [name, integral] : map)
        {
            if(name == key)
                ret.push_back(integral);
        }
    }
    return ret;
}

AngularDistribution::Efficiency::Efficiency(const std::string& file, const std::string& name)
{
    auto* f {new TFile(file.c_str())};
    fEff = f->Get<TEfficiency>(name.c_str());
    if(!fEff)
        throw std::runtime_error("Error: no TEfficiency named" + name + " in file " + file);
    //Init also spline and graph
    fGraph = fEff->CreateGraph();
    fSpe = new TSpline3("effspline", (TGraph*)fGraph, "b2, e2", 0, 0);
    fSpe->SetTitle(fEff->GetTitle());
    //Init uncertainty spline
    InitUncertaintySpline();
}

void AngularDistribution::Efficiency::InitUncertaintySpline()
{
    std::vector<double> x, uct;
    for(int p = 0; p < fGraph->GetN(); p++)
    {
        //X is the same
        x.push_back(fGraph->GetPointX(p));
        //Uncertainty is computed using TGraphAsymmErrors GetErrorY
        //See formula in its documentation
        uct.push_back(fGraph->GetErrorY(p));
    }
    fUSpe = new TSpline3("ueffspeline", &(x[0]), &(uct[0]), x.size(), "b2,e2", 0, 0);
    fUSpe->SetTitle("Uncertainty in efficiency spline");
}

double AngularDistribution::Efficiency::GetPointEff(double thetaCM)
{
    return fSpe->Eval(thetaCM);
}

double AngularDistribution::Efficiency::GetPointUncertainty(double thetaCM)
{
    return fUSpe->Eval(thetaCM);
}

double AngularDistribution::Efficiency::GetAveragedEff(double thetaCMmin, double thetaCMmax)
{
    std::vector<double> vals;
    for(int p = 0; p < fGraph->GetN(); p++)
    {
        if(thetaCMmin <= fGraph->GetPointX(p) && fGraph->GetPointX(p) < thetaCMmax)
            vals.push_back(fGraph->GetPointY(p));
    }
    return TMath::Mean(vals.begin(), vals.end());
}

TCanvas* AngularDistribution::Efficiency::GetCanvas(const std::string& opt) const
{
    auto* cEff {new TCanvas("cEff", "Efficiency canvas")};
    fEff->Draw(("apl" + opt).c_str());
    return cEff;
}

AngularDistribution::ThetaCMIntervals::ThetaCMIntervals(double min, double max, double step, TH1F* hmodel)
{
    int i {};
    for(auto theta = min; theta < max; theta += step)
    {
        //Init intervals
        fVals.push_back({theta, theta + step});
        //Init histograms
        fHistos.push_back(new TH1F(TString::Format("hCM%d", i), TString::Format("#theta_{CM} #in [%.0f, %.0f) #circ;E_{x} [MeV]", theta, theta + step),
                                   hmodel->GetNbinsX(), hmodel->GetXaxis()->GetXmin(), hmodel->GetXaxis()->GetXmax()));
        //Init angle solid elements
        fOmega.push_back(ComputeAngleSolidElement(theta, theta + step));
        i++;
    }
}

void AngularDistribution::ThetaCMIntervals::Fill(double thetaCM, double Ex)
{
    for(int i = 0; i < fVals.size(); i++)
    {
        auto min {fVals[i].first}; auto max {fVals[i].second};
        if(min <= thetaCM && thetaCM < max)
        {
            fHistos[i]->Fill(Ex);
            break;
        }
    }
}

void AngularDistribution::ThetaCMIntervals::FillHisto(int idx, double val)
{
    fHistos.at(idx)->Fill(val);
}

double AngularDistribution::ThetaCMIntervals::ComputeAngleSolidElement(double min, double max)
{
    return TMath::TwoPi() * (TMath::Cos(min * TMath::DegToRad()) - TMath::Cos(max * TMath::DegToRad()));
}

double AngularDistribution::ThetaCMIntervals::GetIntervalCentre(int idx)
{
    return (fVals.at(idx).second + fVals.at(idx).first) / 2;
}

TCanvas* AngularDistribution::ThetaCMIntervals::Draw()
{
    auto* civs {new TCanvas("civs", "Intervals canvas")};
    civs->DivideSquare(fHistos.size());
    for(int h = 0; h < fHistos.size(); h++)
    {
        civs->cd(h + 1);
        fHistos[h]->Draw();
    }
    return civs;
}

AngularDistribution::DiffCrossSection::DiffCrossSection(const std::string& peak, AngularFitter& ang,
                                                        ThetaCMIntervals& ivs,
                                                        Efficiency& eff,
                                                        PhysicsUtils::ExperimentInfo& exp)
{
    //Call to inner function
    PerformCalculation(peak, ang, ivs, eff, exp);
}

void AngularDistribution::DiffCrossSection::PerformCalculation(const std::string& peak, AngularFitter& ang,
                                                               ThetaCMIntervals& ivs,
                                                               Efficiency& eff,
                                                               PhysicsUtils::ExperimentInfo& exp)
{
    //Vector with integrals for all peaks
    auto N {ang.GetIntegralsForPeak(peak)};
    //Init TGraphErrors
    fxs = new TGraphErrors();
    fxs->SetName(TString::Format("gang%s", peak.c_str()));
    fxs->SetTitle(TString::Format("xs for %s state;#theta_{CM} [#circ];#frac{d#sigma}{d#Omega} [mb/sr]", peak.c_str()));
    for(int i = 0; i < N.size(); i++)
    {
        //1->Solid angle element
        auto Omega {ivs.GetOmega(i)};
        //2->Efficiency calculation
        auto epsilon {eff.GetAveragedEff(ivs.GetInterval(i).first, ivs.GetInterval(i).second)};
        //3->Compute differential cross section!
        auto xs {N[i] / (exp.GetNt() * exp.GetNb() * Omega * epsilon)};
        //4->Convert to mb / sr units
        xs *= 1e27;
        //5->Get uncertainty
        auto centre {ivs.GetIntervalCentre(i)};
        auto uncertainty {UncertaintyXS(N[i], exp.GetNt(), exp.GetNb(), Omega, epsilon, centre, eff, exp)};

        //Print info
        std::cout<<"====== XS computation number "<<i<<" ======"<<'\n';
        std::cout<<"N       = "<<N[i]<<'\n';
        std::cout<<"Nt      = "<<exp.GetNt()<<'\n';
        std::cout<<"Nb      = "<<exp.GetNb()<<'\n';
        std::cout<<"Omega   = "<<Omega<<'\n';
        std::cout<<"epsilon = "<<epsilon<<'\n';
        std::cout<<"xs      = "<<xs<<" +/- "<<uncertainty<<" mb / sr"<<'\n';

        //Fill graph
        fxs->SetPoint(fxs->GetN(), centre, xs);
        fxs->SetPointError(fxs->GetN() - 1, 0, uncertainty);
    }
    //Set basic style for graph
    fxs->SetLineWidth(2);
}

double AngularDistribution::DiffCrossSection::UncertaintyXS(double N, double Nt, double Nb, double Omega, double epsilon,
                                                            double thetaCMcentre, Efficiency& eff, PhysicsUtils::ExperimentInfo& exp)
{
    //1-> Uncertainty in N
    double coeffN {1. / (Nt * Nb * Omega * epsilon)};
    double uN {TMath::Sqrt(N)};
    //2-> Uncertainty in Nb
    double coeffNb {- N / (Nt * Omega * epsilon) / TMath::Power(Nb, 2)};
    double uNb {exp.GetUNb()};
    //3-> Uncertainty in epsilon
    double coeffEpsilon {- N / (Nt * Nb * Omega) / TMath::Power(epsilon, 2)};
    double uEpsilon {eff.GetPointUncertainty(thetaCMcentre)};

    //Add everything
    double sum {TMath::Sqrt(coeffN * coeffN * uN * uN +
                            coeffNb * coeffNb * uNb * uNb +
                            coeffEpsilon * coeffEpsilon * uEpsilon * uEpsilon)};
    //Convert to mb units
    sum *= 1e27;
    return sum;
}

TMultiGraph* AngularDistribution::CompareMethods(AngularFitter& ang, const std::string& peak,
                                                 ThetaCMIntervals& ivs, Efficiency& eff, PhysicsUtils::ExperimentInfo& exp)
{
    //1->Integral
    ang.ComputeYield("integral");
    DiffCrossSection xsintegral {peak, ang, ivs, eff, exp};
    auto* gint {xsintegral.GetExperimental()};
    gint->SetTitle("Integral of func");
    gint->SetMarkerStyle(24);
    gint->SetLineWidth(2);
    //2->1 sigma
    ang.ComputeYield("sum", 1);
    DiffCrossSection xs1sigma {peak, ang, ivs, eff, exp};
    auto* g1 {xs1sigma.GetExperimental()};
    g1->SetTitle("1-#sigma");
    g1->SetMarkerStyle(25);
    g1->SetLineWidth(2);
    //3->2 sigma
    ang.ComputeYield("sum", 2);
    DiffCrossSection xs2sigma {peak, ang, ivs, eff, exp};
    auto* g2 {xs2sigma.GetExperimental()};
    g2->SetTitle("2-#sigma");
    g2->SetMarkerStyle(26);
    g2->SetLineWidth(2);
    //Init
    auto* ret {new TMultiGraph()};
    ret->SetTitle("Comparaison of different methods;#theta_{CM} [#circ];#frac{d#sigma}{d#Omega} [mb/sr]");
    ret->Add(gint); ret->Add(g1); ret->Add(g2);
    return ret;
}
#endif
