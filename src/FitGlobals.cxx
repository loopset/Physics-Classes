#include "FitGlobals.h"

bool Fitters::gIsLogL = false;

void Fitters::ToggleIsLogL()
{
    gIsLogL = !gIsLogL;
}

bool Fitters::GetIsLogL()
{
    return gIsLogL;
}
