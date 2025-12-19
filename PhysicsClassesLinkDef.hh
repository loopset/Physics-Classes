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
#pragma link C++ class Fitters::Interface+;
#pragma link C++ class Fitters::CompositeTF1;

// ANGULAR
#pragma link C++ namespace Angular;

#pragma link C++ global Angular::gIsLab;
#pragma link C++ class Angular::Intervals;
#pragma link C++ class Angular::Fitter;
#pragma link C++ class Angular::DifferentialXS;
#pragma link C++ class Angular::Comparator;

// INTERPOLATORS
#pragma link C++ namespace Interpolators;

#pragma link C++ class Interpolators::Efficiency;
#pragma link C++ class Interpolators::Sigmas;

// Utils
#pragma link C++ namespace PhysUtils;

#pragma link C++ class PhysUtils::Experiment;
#pragma link C++ class PhysUtils::Colors;
#pragma link C++ class PhysUtils::SpectroscopicFactor + ;
#pragma link C++ class std::vector < PhysUtils::SpectroscopicFactor> + ;
#pragma link C++ class PhysUtils::SFCollection + ;
#pragma link C++ class PhysUtils::State + ;
#pragma link C++ class std::vector < PhysUtils::State> + ;

// Model plotter
#pragma link C++ namespace PlotUtils;

#pragma link C++ class PlotUtils::ModelToPlot;
#pragma link C++ class PlotUtils::ModelPointers;
#pragma link C++ class PlotUtils::ModelPlotter;

// CALIBRATION
#pragma link C++ namespace Calibration;

#pragma link C++ class Calibration::Source;
#pragma link C++ class Calibration::Runner;

// OMPs
#pragma link C++ namespace PhysOMP;
#pragma link C++ class PhysOMP::OMP;
#pragma link C++ class PhysOMP::Daehnick;
#pragma link C++ class PhysOMP::Haixia;
#pragma link C++ class PhysOMP::DA1p;
#pragma link C++ class PhysOMP::Pang;
#pragma link C++ class PhysOMP::HT1p;
#pragma link C++ class PhysOMP::KoningDelaroche;
#pragma link C++ class PhysOMP::CH89;
#pragma link C++ class PhysOMP::BecchettiGreenless;

// Shell mode information
#pragma link C++ class PhysUtils::QuantumNumbers + ;
#pragma link C++ class PhysUtils::SMData;
#pragma link C++ class PhysUtils::SMParser;

#endif
