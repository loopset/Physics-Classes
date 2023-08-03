#ifndef PhysicsUtils_cxx
#define PhysicsUtils_cxx

#include "PhysicsUtils.h"

#include "TMath.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TSpline.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

double PhysicsUtils::ExperimentInfo::GetUNb() const
{
    if(fDiv == 0)
        std::cout<<"Warning: CATS division factor = 0! Set it!"<<'\n';
    double uct {TMath::Sqrt(fNtriggers) * fDiv *fDiv};
    return uct;
}

std::string PhysicsUtils::Uncertainty::GetStr() const
{
    return {};
}

PhysicsUtils::SigmaInterpolator::SigmaInterpolator(const std::string& file, const std::string& name)
{
    auto* f {new TFile(file.c_str())};
    fGraph = f->Get<TGraphErrors>(name.c_str());
    if(!fGraph)
        throw std::runtime_error("Could not read file or TGraphErrors in file " + file);
    f->Close(); delete f;
    //Init spline
    fSpe = new TSpline3("effspline", (TGraph*)fGraph, "b2, e2", 0, 0);
}

void PhysicsUtils::SigmaInterpolator::ComputeAndSetScalingFactorInGS(double expsimgags)
{
    double simusigma {fGraph->GetPointY(0)};
    double scaling {expsimgags / simusigma};
    SetScalingFactor(scaling);
    std::cout<<"=============================================="<<'\n';
    std::cout<<"Simulated sigma for g.s    = "<<simusigma<<'\n';
    std::cout<<"Experimental sigma for g.s = "<<expsimgags<<'\n';
    std::cout<<"Scaling factor             = "<<scaling<<'\n';
    std::cout<<"=============================================="<<'\n';
}

void PhysicsUtils::SigmaInterpolator::SetScalingFactor(double scaling)
{
    fFactor = scaling;
    fScaled = (TGraphErrors*)fGraph->Clone();
    fScaled->SetNameTitle("gscaled", "Scaled sigma graph;E_{x} [MeV];#sigma_{Scaled} [MeV]");
    fScaled->Scale(fFactor, "y");
    //overwrite spline
    if(fSpe)
        delete fSpe;
    fSpe = new TSpline3("effspline", (TGraph*)fScaled, "b2, e2", 0, 0);
}

#endif
