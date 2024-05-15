#include "PhysColors.h"

#include "TBox.h"
#include "TCanvas.h"
#include "TText.h"

PhysUtils::Colors* PhysUtils::Colors::fInstance = nullptr;

PhysUtils::Colors* PhysUtils::Colors::GetInstance()
{
    if(!fInstance)
        fInstance = new Colors;
    return fInstance;
}

void PhysUtils::Colors::Init()
{
    // CVD-friendly, defined in https://personal.sron.nl/~pault/#sec:qualitative

    std::vector<double> r, g, b;
    r = {68, 102, 34, 204, 238, 170, 187};
    g = {119, 204, 136, 187, 102, 51, 187};
    b = {170, 238, 51, 68, 119, 119, 187};
    for(int c = 0; c < r.size(); c++)
    {
        auto i {TColor::GetFreeColorIndex()};
        fColors.push_back(new TColor(i, r[c] / 256, g[c] / 256, b[c] / 256));
    }
};

int PhysUtils::Colors::operator[](int i) const
{
    if(i >= fColors.size())
        return 1; // default to black
    return fColors.at(i)->GetNumber();
}

void PhysUtils::Colors::Draw() const
{
    auto* c {new TCanvas {"cPhysColors", "PhysUtils::Colors canvas"}};
    double x1 {0};
    double y1 {0};
    double x2 {20};
    double y2 {1};

    gPad->SetFillColor(0);
    gPad->Clear();
    gPad->Range(x1, y1, x2, y2);

    TText text(0, 0, "");
    text.SetTextFont(61);
    text.SetTextSize(0.07);
    text.SetTextAlign(22);

    TBox box;

    // Draw color table boxes.
    double hs = (y2 - y1) / 1;
    double ws = (x2 - x1) / fColors.size();
    for(int i = 0; i < fColors.size(); i++)
    {
        double xlow = x1 + ws * (i + 0.1);
        double xup = x1 + ws * (i + 0.9);
        double ylow = y1 + hs * 0.1;
        double yup = y1 + hs * 0.9;
        box.SetFillStyle(1001);
        int color {fColors.at(i)->GetNumber()};
        box.SetFillColor(fColors.at(i)->GetNumber());
        box.DrawBox(xlow, ylow, xup, yup);
        box.SetFillStyle(0);
        box.SetLineColor(1);
        box.DrawBox(xlow, ylow, xup, yup);
        if(color == 1)
            text.SetTextColor(0);
        else
            text.SetTextColor(1);
        text.DrawText(0.5 * (xlow + xup), 0.5 * (ylow + yup), TString::Format("%d", i).Data());
    }
}
