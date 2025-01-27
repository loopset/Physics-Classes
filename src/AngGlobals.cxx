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
