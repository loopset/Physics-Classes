#include "AngDifferentialXS.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TString.h"

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
        auto eps {fEff->GetPointEff(peak, thetaCM)};
        // 3-> Value!
        auto xs {N[iv] / (fExp->GetNt() * fExp->GetNb() * Omega * eps)};
        // 4-> Convert to mb / sr units
        xs *= 1e27;
        // Fill graph
        fXS[peak]->SetPoint(fXS[peak]->GetN(), thetaCM, xs);
        // With uncertainty
        auto unc {Uncertainty(N[iv], fExp->GetNt(), fExp->GetNb(), Omega, eps, thetaCM)};
        fXS[peak]->SetPointError(fXS[peak]->GetN() - 1, 0, unc);
    }
}

double Angular::DifferentialXS::Uncertainty(double N, double Nt, double Nb, double Omega, double eps, double thetaCM)
{
    // 1-> N
    double coeffN {1. / (Nt * Nb * Omega * eps)};
    double uN {TMath::Sqrt(N)};
    // 2-> Nb
    double coeffNb {-N / (Nt * Omega * eps) / TMath::Power(Nb, 2)};
    double uNb {fExp->GetUNb()};
    // 3-> Eps
    double coeffEps {-N / (Nt * Nb * Omega) / TMath::Power(eps, 2)};
    double uEps {0};

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

TCanvas* Angular::DifferentialXS::Draw() const
{
    auto* c {new TCanvas {"cDiffXS", "Differential XS"}};
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
