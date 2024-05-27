#include "AngDifferentialXS.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"

#include "PhysColors.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

void Angular::DifferentialXS::Do(const std::string& peak)
{
    // Get counts by interval for peak
    auto N {fFitter->GetIgCountsFor(peak)};
    // Init TGraphErrors
    fXS[peak] = new TGraphErrors;
    fXS[peak]->SetNameTitle(TString::Format("xs%s", peak.c_str()),
                            TString::Format("%s;#theta_{CM} [#circ];d#sigma / d#Omega [mb / sr]", peak.c_str()));
    // Run for each interval
    for(int iv = 0; iv < N.size(); iv++)
    {
        // 0-> ThetaCM interval center
        auto thetaCM {fIvs->GetCenter(iv)};
        // 1-> Solid angle element
        auto Omega {fIvs->GetOmega(iv)};
        // 2-> Efficiency
        auto eps {fEff->GetMeanEff(peak, fIvs->GetLow(iv), fIvs->GetUp(iv))};
        // std::cout << "Mean eps : " << eps << " puntual : " << fEff->GetPointEff(peak, thetaCM) << '\n';
        // auto eps {fEff->GetPointEff(peak, thetaCM)};
        // Skip calculation if N[iv] is below threshold
        if(N[iv] < fNThresh)
        {
            std::cout << BOLDCYAN << "DifferentialXS::Do(): N = " << N[iv] << " < " << fNThresh << " for peak " << peak
                      << '\n';
            std::cout << "-> ThetaCM : " << thetaCM << '\n';
            std::cout << "-> Omega   : " << Omega << '\n';
            std::cout << "-> Eff     : " << eps << '\n';
            std::cout << "   Skipping point!" << '\n';
            std::cout << RESET;
            continue;
        }
        // 3-> Value!
        auto xs {N[iv] / (fExp->GetNt() * fExp->GetNb() * Omega * eps)};
        if(!std::isfinite(xs))
        {
            std::cout << BOLDRED << "DifferentialXS::Do(): NaN found for " << peak << '\n';
            std::cout << "-> xs      : " << xs << '\n';
            std::cout << "-> ThetaCM : " << thetaCM << '\n';
            std::cout << "-> Omega   : " << Omega << '\n';
            std::cout << "-> Eff     : " << eps << '\n';
            std::cout << "   Skipping point!" << '\n';
            std::cout << RESET;
            continue;
        }
        // 4-> Convert to mb / sr units
        xs *= 1e27;
        // Fill graph
        fXS[peak]->SetPoint(fXS[peak]->GetN(), thetaCM, xs);
        // With uncertainty
        auto unc {Uncertainty(peak, N[iv], fExp->GetNt(), fExp->GetNb(), Omega, eps, thetaCM)};
        fXS[peak]->SetPointError(fXS[peak]->GetN() - 1, 0, unc);
    }
}

double Angular::DifferentialXS::Uncertainty(const std::string& peak, double N, double Nt, double Nb, double Omega,
                                            double eps, double thetaCM)
{
    // 1-> N
    double coeffN {1. / (Nt * Nb * Omega * eps)};
    double uN {TMath::Sqrt(N)};
    // 2-> Nb
    double coeffNb {-N / (Nt * Omega * eps) / TMath::Power(Nb, 2)};
    double uNb {fExp->GetUNb()};
    // 3-> Eps
    double coeffEps {-N / (Nt * Nb * Omega) / TMath::Power(eps, 2)};
    double uEps {fEff->GetPointUEff(peak, thetaCM)};

    // Add everything
    double sum {
        TMath::Sqrt(coeffN * coeffN * uN * uN + coeffNb * coeffNb * uNb * uNb + coeffEps * coeffEps * uEps * uEps)};
    // Convert to mb/sr units
    sum *= 1e27;
    return sum;
}

void Angular::DifferentialXS::DoFor(const std::vector<std::string>& peaks)
{
    for(const auto& peak : peaks)
        Do(peak);
}

TCanvas* Angular::DifferentialXS::Draw(const TString& title) const
{
    static int cXSIdx {};
    auto* c {new TCanvas {TString::Format("cXS%d", cXSIdx), (title.Length()) ? title : "Angular::DifferentialXS"}};
    cXSIdx++;
    c->DivideSquare(fXS.size());
    int idx {1};
    for(const auto& [name, g] : fXS)
    {
        c->cd(idx);
        g->SetLineWidth(2);
        g->SetMarkerStyle(25);
        g->Draw("ap pmc plc");
        idx++;
    }
    return c;
}

void Angular::DifferentialXS::Write(const std::string& dir, const std::string& name) const
{
    // Run for each peak
    for(const auto& [peak, g] : fXS)
    {
        // 1-> Init output file
        std::ofstream streamer {(dir + peak + "_xs.dat")};
        if(!streamer)
            throw std::runtime_error("DifferentialXS::Write: cannot open directory " + dir);
        // 2-> Push back from graph points
        for(auto i = 0; i < g->GetN(); i++)
            streamer << g->GetPointX(i) << "  " << g->GetPointY(i) << "  " << g->GetErrorY(i) << '\n';
        // 3-> Close
        streamer.close();
    }
    // Save also in root format
    auto* fout {new TFile {(dir + name + ".root").c_str(), "recreate"}};
    for(const auto& [peak, g] : fXS)
        g->Write(("g" + peak).c_str());
    fout->Close();
}
