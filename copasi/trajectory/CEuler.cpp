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
mTargetTime(0.0), 
mData(), 
mpY(NULL), 
mpYdot(NULL),
mpYd(NULL)
{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
  initializeParameter();
}

CEulerMethod::CEulerMethod(const CEulerMethod & src,
                           const CDataContainer * pParent): 
CTrajectoryMethod(src, pParent),
mTargetTime(0.0), 
mData(), 
mpY(NULL), 
mpYdot(NULL)
{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
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
  mData.dim = (C_INT)(mContainerState.size() - mpContainer->getCountFixedEventTargets());
  pdeletev(mpY);
  mpY = new C_FLOAT64[mData.dim];
  pdeletev(mpYd);
  mpYd = new C_FLOAT64[mData.dim];
  mTargetTime = *mpContainerStateTime + pTP->getDuration(); //mTargettime wird die Dauer des Problems zugewiesen 
  mStepsize = getValue< double >("Step size"); //getValue oder getParameeter?
  mpYdot = mpContainer->getRate(*mpReducedModel).array() + mpContainer->getCountFixedEventTargets();
  memcpy(mpY, mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));
}


/**
   * from Radau5 - also used after start
   */
// from Radau5 but just evalF would also be fine I think ????
void CEulerMethod::EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot, C_FLOAT64 *, C_INT *)
{
  static_cast<Data *>((void *) n)->pMethod->evalF(t, y, ydot);
}

void CEulerMethod::evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot)
{
  //CVector< C_FLOAT64 > yTemp(mData.dim);
  //memcpy(yTemp.array(), mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));

  if (y != mpContainerStateTime)
    memcpy(mpContainerStateTime, y, mData.dim * sizeof(C_FLOAT64));

  mpContainer->updateSimulatedValues(*mpReducedModel);
  memcpy(ydot, mpYdot, mData.dim * sizeof(C_FLOAT64));

#ifdef DEBUG_NUMERICS
  std::cout << "State:     " << mpContainer->getState(false) << std::endl;
  std::cout << "Rate:      " << mpContainer->getRate(false) << std::endl;
#endif // DEBUG_NUMERICS

  //memcpy(mpContainerStateTime, yTemp.array(), mData.dim * sizeof(C_FLOAT64));

  return;
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
  evalF(&startTime, mpY, mpYd); 

  for (int i = 0; i < mData.dim; ++i)
{
   mpY[i] = mpY[i] + h * mpYd[i];
}
  // Zeit vorschieben
  *mpContainerStateTime += h;

  memcpy(mpContainerStateTime, mpY, mData.dim * sizeof(C_FLOAT64));
  mpContainer->updateSimulatedValues(false);

return *mpContainerStateTime;

}