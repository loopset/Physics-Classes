#ifndef AngGlobals_h
#define AngGlobals_h

namespace Angular
{

extern bool gIsLab; //!< Global extern variable to represent whether xs calculations are in CM or Lab

extern bool gUseHessErrors; //!< Activate Hessian errors in AngFitter

void ToggleIsLab();

bool GetIsLab();

void ToggleHessErrors();

bool GetUseHessErrors();

} // namespace Angular
#endif // !AngGlobals_h
