#include "copasi/copasi.h"
#include "CEuler.h"
#include "copasi/core/CDataVector.h"
#include "copasi/function/CFunction.h"
#include "copasi/randomGenerator/CRandom.h"
#include "CTrajectoryMethod.h"
#include "CTrajectoryProblem.h"
#include "copasi/math/CMathContainer.h"
#include "copasi/model/CState.h"
#include "copasi/model/CCompartment.h"
#include "copasi/model/CModel.h"
#include "copasi/utilities/CCopasiMethod.h"
#include "copasi/model/CModel.h"

//Konstruktoren von directmethod übernommen und rausgenommen, was nicht gebraucht wird - eigentlich nur mTargettime
CEulerMethod::CEulerMethod(const CDataContainer * pParent,
                           const CTaskEnum::Method & methodType,
                           const CTaskEnum::Task & taskType): 
        CTrajectoryMethod(pParent, methodType, taskType), //Basisklassenkosntruktot - war bei stochdirect so
        mTargetTime(0.0)
{
  initializeParameter();
}

CEulerMethod::CEulerMethod(const CEulerMethod & src, const CDataContainer * pParent)
: CTrajectoryMethod(src, pParent)
{
  initializeParameter();
}

CEulerMethod::~CEulerMethod()
{
  //einfach leer 
}

void CEulerMethod::initializeParameter()
{
  assertParameter("Step size", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.01);
}

void CEulerMethod::start() 
{
  CTrajectoryMethod::start();
  const CTrajectoryProblem * pTP = static_cast<const CTrajectoryProblem *>(mpProblem); //static cast gibt zugriff auf andere Basisklasse 
  mTargetTime = *mpContainerStateTime + pTP->getDuration(); //mTargettime wird die Dauer des Problems zugewiesen 
  mStepsize = getValue< double >("Step size"); //getValue oder getParameeter?
}


CTrajectoryMethod::Status CEulerMethod::step(const double & deltaT,
                                             const bool & /* final */) //nicht verhandelbar, weil so in der basisklasse
{
  C_FLOAT64 EndTime = *mpContainerStateTime + deltaT;
  //Calculation of step amount 

  //PROBLEM: wenn letzter Schritt größer ist 
  while (*mpContainerStateTime<EndTime)
  {
    doSingleStep(*mpContainerStateTime, *mpContainerStateTime + mStepsize);
  }
  return NORMAL; //muss glaub weil bei step in TrajectoryMethods normal FAILURE zurückgegeben wird standartmäßig
}


// Euler Schritt um durch math container durchzuitterieren y(t+1, y+1) = y(t,y) + h * f`(t,y) 
// ich muss auf die Startwerte zugreifen 
// und auf die rechte Seite der ODEs 
// do one step würde dann die startwerte in die rechte seite einsetzen, um die Steigung zu bekommen bedeutet yn + h * f(yn)
C_FLOAT64 CEulerMethod::doSingleStep(C_FLOAT64 startTime, const C_FLOAT64 & endTime)
{
  C_FLOAT64 h = mStepsize;
  // y müsste eigentlich mpContainer sein 
  

  

  // 2. Aktualisiere rechte Seite (d.h. berechne Ableitungen)
  //mpMathContainer->getValues();

  // 3. Zugriff auf Zustand und Ableitungen
  //std::vector<C_FLOAT64> & y = mpMathContainer->getState();
  //std::vector<C_FLOAT64> & ydot = mpMathContainer->getDerivatives();

  /* 4. Euler-Schritt: y = y + h * ydot
  for (size_t i = 0; i < y.size(); ++i)
  {
    y[i] += h * ydot[i];


  }

  */ // 5. Zeit vorschieben
  *mpContainerStateTime += h;

  return *mpContainerStateTime;
}