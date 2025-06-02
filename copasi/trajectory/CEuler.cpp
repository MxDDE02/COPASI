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
#include <vector>
#include <iomanip>  // für std::setprecision



CEulerMethod::CEulerMethod(const CDataContainer * pParent,
                           const CTaskEnum::Method & methodType,
                           const CTaskEnum::Task & taskType): 
CTrajectoryMethod(pParent, methodType, taskType),
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
mIntervalSize(0.0),
mData(), 
mpY(NULL), 
mpYdot(NULL), 
mpYd(NULL)
{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
  initializeParameter();
}

CEulerMethod::~CEulerMethod()
{
  pdeletev(mpYd);
  pdeletev(mpY);
}

void CEulerMethod::initializeParameter()
{
  assertParameter("Step size", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.01);
  assertParameter("Interpolation", CCopasiParameter::Type::BOOL, true);;
}

void CEulerMethod::start()
{
  // 1. Call base class method to initialize container state and time
  CTrajectoryMethod::start();

  // 2. Retrieve the trajectory problem for duration and step size
  const CTrajectoryProblem * pTP = static_cast<const CTrajectoryProblem *>(mpProblem);

  // 3. Determine the dimension of the system (excluding fixed event targets)
  mData.dim = (C_INT)(mContainerState.size() - mpContainer->getCountFixedEventTargets());

  // 4. Get pointer to rate vector (excluding fixed event targets)
  mpYdot = mpContainer->getRate(*mpReducedModel).array() + mpContainer->getCountFixedEventTargets();

  // 5. Set simulation end time and initial step interval
  mTargetTime = *mpContainerStateTime + pTP->getDuration();
  mIntervalSize = *mpContainerStateTime + pTP->getStepSize();

  // 6. Retrieve the integration step size from parameters
  mStepsize = getValue< double >("Step size");

  // 7. Allocate memory for state (mpY) and derivative (mpYd) vectors
  mpY = new C_FLOAT64[mData.dim];
  mpYd = new C_FLOAT64[mData.dim];

  // 8. Copy current state into local state vector
  memcpy(mpY, mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));

  // 9. Ensure step size does not exceed total interval size 
  /*
  if (mStepsize > mIntervalSize)
  {
    mStepsize = mIntervalSize;
  }
  */
}

/* Uncomment this - if the integrator is written in Fortan 
void CEulerMethod::EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot, C_FLOAT64 *, C_INT *)
{
  static_cast<Data *>((void *) n)->pMethod->evalF(t, y, ydot);
}
*/

//Uncoment this - if the evaluation of the math container should be in a seperate function (unnecessary in Euler)
/*
void CEulerMethod::evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot)
{
  //The comments just make a backup version of the current mpContainerStateTime 
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
*/

CTrajectoryMethod::Status CEulerMethod::step(const double & deltaT, const bool & /* final */)
{
  const C_FLOAT64 outputTime = *mpContainerStateTime + deltaT;

  //Saving of the initial Math-Container in the mHistory - for interpolation
  TimeStatePair initial;
  initial.time = *mpContainerStateTime;
  initial.state.assign(mpY, mpY + mData.dim);
  initial.rate.assign(mpYdot, mpYdot + mData.dim);
  mHistory.push_back(initial);


  while (*mpContainerStateTime < outputTime)
  {
    doSingleStep(*mpContainerStateTime);
  }
  
  //Interpolation -> get the Trajectory Problem defined state at the requested time 
  std::vector<C_FLOAT64> interpolatedState = interpolateAt(outputTime);
  memcpy(mpContainerStateTime, interpolatedState.data(), mData.dim * sizeof(C_FLOAT64));


  *mpContainerStateTime = outputTime;

  return NORMAL;
}

C_FLOAT64 CEulerMethod::doSingleStep(C_FLOAT64 startTime)
{
  // Calclate the rates => evalF 
  if (mpY != mpContainerStateTime)
  memcpy(mpContainerStateTime, mpY, mData.dim * sizeof(C_FLOAT64));
  mpContainer->updateSimulatedValues(false);
  memcpy(mpYd, mpYdot, mData.dim * sizeof(C_FLOAT64)); 

  // Euler-Step: y = y + h * f(y)
    for (int i = 0; i < mData.dim; ++i)
    {
      mpY[i] = mpY[i] + mStepsize* mpYd[i];
    }

    // Update the time 
    *mpContainerStateTime += mStepsize;

    // writes the result in the math container 
    memcpy(mpContainerStateTime, mpY, mData.dim * sizeof(C_FLOAT64));
    //mpContainer->updateSimulatedValues(false);

    // Save each Euler-Step result in the mHistory 
    TimeStatePair ts;
    ts.time = *mpContainerStateTime;
    ts.state.assign(mpY, mpY + mData.dim);
    ts.rate.assign(mpYd, mpYd + mData.dim);
    mHistory.push_back(ts);
  }



std::vector<C_FLOAT64> CEulerMethod::interpolateAt(C_FLOAT64 t) const
{
  //does nothing if no integration has been done 
  if (mHistory.empty()) return {}; 
  // precaution: if the time is smaller then the starttime (=0) the starting state will be returned 
  if (t <= mHistory.front().time) return mHistory.front().state; 
  // prectaution: if the time is bigger then the endtime , the end state will be returned 
  if (t >= mHistory.back().time) return mHistory.back().state;
 
  // Find the interpolation - intervall borders 
  for (size_t i = 1; i < mHistory.size(); ++i)
  {
    // precaution: if a entry for the requested time already exists
    if (mHistory[i].time == t) return mHistory[i].state;
    //search for the first time point, which is bigger then t => right border 
    else if (mHistory[i].time > t)
    {
      // time points which should interpolate 
      const TimeStatePair& p0 = mHistory[i - 1]; //left border 
      const TimeStatePair& p1 = mHistory[i]; //right border 
      //if Euler interpolation is requested 
      //C_FLOAT64 h0 = t - p0.time; 

      //vector for saving the interpolated state at time t 
      std::vector<C_FLOAT64> interpolated(mData.dim);
      for (int j = 0; j < mData.dim; ++j)
      {
        //linear interpolation: euler step with different step size 
        //interpolated[j] = p0.state[j] + h0 * p0.rate[j];

        //linear interpolation y = y1 + ((x-x1)/(x2-x1)) * (y2-y1)
        interpolated[j] = p0.state[j] + ((t - p0.time) / (p1.time - p0.time)) * (p1.state[j] - p0.state[j]);

      } 

      return interpolated;
    }
  }
}
