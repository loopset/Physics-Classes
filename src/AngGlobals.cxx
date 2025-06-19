#include "AngGlobals.h"

bool Angular::gIsLab = false;

void Angular::ToggleIsLab()
{
    gIsLab = !gIsLab;
}

bool Angular::GetIsLab()
{
    return gIsLab;
}

bool Angular::gUseHessErrors = false;

void Angular::ToggleHessErrors()
{
    gUseHessErrors = !gUseHessErrors;
}

bool Angular::GetUseHessErrors()
{
    return gUseHessErrors;
}
