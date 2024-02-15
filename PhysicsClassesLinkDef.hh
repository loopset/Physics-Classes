#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

// FIT
#pragma link C++ namespace Fitters;

#pragma link C++ class Fitters::Data;
#pragma link C++ class Fitters::Model;
#pragma link C++ class Fitters::Objective;
#pragma link C++ class Fitters::Runner;
#pragma link C++ class Fitters::Plotter;

// Model plotter
#pragma link C++ namespace PlotUtils;

#pragma link C++ class PlotUtils::ModelToPlot;
#pragma link C++ class PlotUtils::ModelPointers;
#pragma link C++ class PlotUtils::ModelPlotter;


#endif
