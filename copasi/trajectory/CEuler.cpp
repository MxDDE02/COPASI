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

CEulerMethod::CEulerMethod(const CDataContainer * pParent,
                           const CTaskEnum::Method & methodType,
                           const CTaskEnum::Task & taskType): 
CTrajectoryMethod(pParent, methodType, taskType),
mData(), 
mpY(NULL), 
mpYdot(NULL),
mpYd(NULL), 
interpolated(NULL) 
{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
  initializeParameter();
}

CEulerMethod::CEulerMethod(const CEulerMethod & src,
                           const CDataContainer * pParent): 
CTrajectoryMethod(src, pParent),
mData(), 
mpY(NULL), 
mpYdot(NULL), 
mpYd(NULL), 
interpolated(NULL) 
{
  assert((void *) &mData == (void *) &mData.dim);
  mData.pMethod = this;
  initializeParameter();
}

CEulerMethod::~CEulerMethod()
{
  pdeletev(mpYd);
  pdeletev(mpY);
  pdeletev(interpolated); 
}

void CEulerMethod::initializeParameter()
{
  assertParameter("initial step size", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.01);
  assertParameter("absolute tolerance", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.000001);
  assertParameter("relative tolerance", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.0001);
  assertParameter("Maximal internal steps", CCopasiParameter::Type::INT, 1000000);

}

bool CEulerMethod::elevateChildren()
{
  initializeParameter();
  return true;
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
  int NumRoots = mpContainer->getRoots().size();

  // 6. Retrieve input parameters from the Trajectory Problem 
  mStepsize = getValue< double >("initial step size");
  euler_atolerance = getValue< double >("absolute tolerance");
  euler_rtolerance = getValue< double >("relative tolerance");
  steplimit = getValue< int >("Maximal internal steps");

  // 7. Allocate memory 
  mpY = new C_FLOAT64[mData.dim];
  mpYd = new C_FLOAT64[mData.dim];
  interpolated = new C_FLOAT64[mData.dim];

  // set outputime (to zero)
  outputTime = *mpContainerStateTime;

  // 8. Copy current state into local state vector
  memcpy(mpY, mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));

  
}

bool CEulerMethod::isValidProblem(const CCopasiProblem * pProblem)
{
  // if the problem is not a trajectory problem 
  if (!CTrajectoryMethod::isValidProblem(pProblem)) return false;

  // creates pTP to read - not change - the problem 
  const CTrajectoryProblem * pTP = dynamic_cast<const CTrajectoryProblem *>(pProblem);

  //checks if the duration is realistic 
  if (pTP->getDuration() < 0.0)
    {
      CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 9);
      return false;
    }
  // check if the initial step size is positive 
  if (getValue< C_FLOAT64 >("initial step size") < 0)
    {
      CCopasiMessage(CCopasiMessage::ERROR, "Stepsize must be positive");
      return false;
    }

  return true;
}

/* Uncomment this - if the integrator is written in Fortran 
void CEulerMethod::EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot, C_FLOAT64 *, C_INT *)
{
  static_cast<Data *>((void *) n)->pMethod->evalF(t, y, ydot);
}
*/

void CEulerMethod::evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot)
{
  //The comments just make a backup version of the current mpContainerStateTime 
  CVector< C_FLOAT64 > yTemp(mData.dim);
  memcpy(yTemp.array(), mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));

  if (y != mpContainerStateTime)
    memcpy(mpContainerStateTime, y, mData.dim * sizeof(C_FLOAT64));
  //this does the actual evaluation and puts it into ydot 
  mpContainer->updateSimulatedValues(*mpReducedModel);
  memcpy(ydot, mpYdot, mData.dim * sizeof(C_FLOAT64));
//this is old debugging 
#ifdef DEBUG_NUMERICS
  std::cout << "State:     " << mpContainer->getState(false) << std::endl;
  std::cout << "Rate:      " << mpContainer->getRate(false) << std::endl;
#endif // DEBUG_NUMERICS
//retreave the old value 
  memcpy(mpContainerStateTime, yTemp.array(), mData.dim * sizeof(C_FLOAT64));

  return;
}

void CEulerMethod::evalR(const C_FLOAT64 * t, const C_FLOAT64 *  /* y */,
                          const C_INT *  nr, C_FLOAT64 * r)
{
  *mpContainerStateTime = *t;
  mpContainer->updateRootValues(*mpReducedModel);

  CVectorCore< C_FLOAT64 > RootValues(*nr, r);
  RootValues = mpContainer->getRoots();
};

CTrajectoryMethod::Status CEulerMethod::step(const double & deltaT, const bool & /* final */)
{
  // set the inital parameters 
  outputTime = *mpContainerStateTime + deltaT;
  int internalsteps = 0; 

  //Saving of the initial Math-Container in the mHistory - for interpolation
  TimeStatePair initial;
  initial.time = *mpContainerStateTime;
  initial.state.assign(mpY, mpY + mData.dim);
  mHistoryinter.push_back(initial);

  //evaluate the roots at the start (x1, y1)
  const CVector< C_FLOAT64 >& rootsstart = mpContainer->getRoots();
  TimeStatePair initialroot;
  initialroot.time = *mpContainerStateTime;
  initialroot.state.assign(rootsstart.array(), rootsstart.array() + rootsstart.size());
  mHistoryrootsstart.push_back(initialroot);


  //the actual integration + stop if the internalsteps exceed the maximal step limit 
  while (*mpContainerStateTime < outputTime)
  {
    doOneStep(*mpContainerStateTime);
    internalsteps++; 
    if(internalsteps >= steplimit)
    {
     CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 12);
     return FAILURE; 
    }
  } 
  
  //Interpolation -> get the Trajectory Problem defined state at the requested time 
  std::vector<C_FLOAT64> interpolatedState = interpolateAttime(outputTime);
  //Update everything to the interpolated output 
  memcpy(mpContainerStateTime, interpolatedState.data(), mData.dim * sizeof(C_FLOAT64));
  memcpy(mpY, mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));
  *mpContainerStateTime = outputTime;
  //updating the state and the roots -> for root finding 
  mpContainer->updateSimulatedValues(false);
  mpContainer->updateRootValues(false);


  //evaluate the roots at the interpolated time (x2, y2)
  const CVector< C_FLOAT64 >& rootsinterpolated = mpContainer->getRoots(); //this is why we needed to update the roots previously 
  TimeStatePair interpolatedroot;
  interpolatedroot.time = *mpContainerStateTime;
  interpolatedroot.state.assign(rootsinterpolated.array(), rootsinterpolated.array() + rootsinterpolated.size());
  mHistoryrootsinterpolate.push_back(interpolatedroot);

  //evaluate if the sign changed -> if yeah we have get the root  
  std::vector<C_FLOAT64> rootTimes;

  for (size_t i = 0; i < mHistoryrootsinterpolate.size(); ++i)
  {
    const TimeStatePair& x0 = mHistoryrootsstart[i];
    const TimeStatePair& x1 = mHistoryrootsinterpolate[i];

    for (size_t j = 0; j < x0.state.size(); ++j)
    {
      if (x0.state[j] * x1.state[j] < 0)
      {
       // calculate time for root 
        C_FLOAT64 t_root = x0.time - x0.state[j] * ((x1.time - x0.time) / (x1.state[j] - x0.state[j]));

        // save the rootvalue (time)
        rootTimes.push_back(t_root);
      }
    }
  }
  if (!rootTimes.empty())
  {
    // find the first event (=smalles time point) with their index 
    C_FLOAT64 minVal = *std::min_element(rootTimes.begin(), rootTimes.end());
    size_t minIndex = std::distance(rootTimes.begin(), std::min_element(rootTimes.begin(), rootTimes.end()));


    // std::vector<C_FLOAT64> interpolatedrootState = interpolateAttime(minVal);
    // // //Update everything to the interpolated output 
    // memcpy(mpContainerStateTime, interpolatedrootState.data(), mData.dim * sizeof(C_FLOAT64));
    // memcpy(mpY, mpContainerStateTime, mData.dim * sizeof(C_FLOAT64));
    // *mpContainerStateTime = minVal;
    // //updating the state and the roots -> for root finding 
    //  mpContainer->updateSimulatedValues(false);
    //  mpContainer->updateRootValues(false);
    //  CMathUpdateSequence;
    //  mpContainer->applyUpdateSequence(minIndex); 
  }



  
  //clearing the mHistory for next steps 
  mHistoryrootsstart.clear(); 
  mHistoryrootsinterpolate.clear(); 
  mHistoryinter.clear();
  return NORMAL;
}


bool CEulerMethod::doOneStep(C_FLOAT64 startTime)
{
  // Error calculation, until error is small 
  bool accepted = false;
  while (!accepted)
  {
    std::vector<C_FLOAT64> y_original(mpY, mpY + mData.dim);
    C_FLOAT64 t_old = *mpContainerStateTime;

    std::vector<C_FLOAT64> fullstep(mData.dim);
    std::vector<C_FLOAT64> halfstep(mData.dim);
    std::vector<C_FLOAT64> deltaerror(mData.dim);
    std::vector<C_FLOAT64> scale(mData.dim);

    // == 1. calculate rates ==
    evalF(mpContainerStateTime, mpY, mpYd);

    // == 2. full step ==
    for (int i = 0; i < mData.dim; ++i)
      fullstep[i] = y_original[i] + mStepsize * mpYd[i];
    // == 2. full step ==
    for (int i = 0; i < mData.dim; ++i)
      fullstep[i] = y_original[i] + mStepsize * mpYd[i];

    // == 3. first half step ==
    for (int i = 0; i < mData.dim; ++i)
      mpY[i] = y_original[i] + (mStepsize / 2.0) * mpYd[i];
    // == 3. first half step ==
    for (int i = 0; i < mData.dim; ++i)
      mpY[i] = y_original[i] + (mStepsize / 2.0) * mpYd[i];

    *mpContainerStateTime = t_old + mStepsize / 2.0;

    evalF(mpContainerStateTime, mpY, mpYd);
    // memcpy(mpContainerStateTime, mpY, mData.dim * sizeof(C_FLOAT64));
    // mpContainer->updateSimulatedValues(false);
    // memcpy(yd_temp.data(), mpYdot, mData.dim * sizeof(C_FLOAT64));

    // == 4. second half step ==
    for (int i = 0; i < mData.dim; ++i)
      halfstep[i] = mpY[i] + (mStepsize / 2.0) * mpYd[i];

    // == 5. error estimation ==
    C_FLOAT64 localerror = 0.0;
    for (int i = 1; i < mData.dim; ++i)
    {
      deltaerror[i] = std::abs(halfstep[i] - fullstep[i]);
      scale[i] = euler_atolerance + std::max(std::abs(y_original[i]), std::abs(fullstep[i])) * euler_rtolerance;
      localerror += std::abs(deltaerror[i] / scale[i]);
    }
    localerror = localerror/(mData.dim-1); 

    // == 6. decition if step can be accepted ==
    if (localerror <= 1.0)
    {
      for (int i = 0; i < mData.dim; ++i)
        mpY[i] = fullstep[i];

      *mpContainerStateTime = t_old + mStepsize;
      memcpy(mpContainerStateTime, mpY, mData.dim * sizeof(C_FLOAT64));
      mpContainer->updateSimulatedValues(false);

      TimeStatePair ts;
      ts.time = *mpContainerStateTime;
      ts.state.assign(mpY, mpY + mData.dim);
      if(std::abs(ts.time - outputTime)<mStepsize)
      {
        mHistoryinter.push_back(ts);
      }

      accepted = true; 
    }
    else
    {
      // step is not accepted -> stepsize will get reduced 
      mStepsize = mStepsize * std::sqrt(euler_rtolerance / localerror); 
      //bring back original state
      memcpy(mpY, y_original.data(), mData.dim * sizeof(C_FLOAT64));
      *mpContainerStateTime = t_old;
      // if the stepsize is small so that Euler just takes ages - the intervalsteplimit takes care of that 
    }
  }
}

std::vector<C_FLOAT64> CEulerMethod::interpolateAttime(C_FLOAT64 t) const
{
  if (mHistoryinter.empty()) 
  {
    CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 32);
    return {}; 
  }

  if (t < mHistoryinter.front().time || t > mHistoryinter.back().time)
  {
    CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 32);
    return {}; 
  }

  for (size_t i = 1; i < mHistoryinter.size(); ++i)
  {
    if (mHistoryinter[i].time >= t)
    {
      const TimeStatePair& p0 = mHistoryinter[i - 1];
      const TimeStatePair& p1 = mHistoryinter[i];

      for (int j = 0; j < mData.dim; ++j)
      {
        interpolated[j] = p0.state[j] + ((t - p0.time) / (p1.time - p0.time)) * (p1.state[j] - p0.state[j]);
      }

      // retun: the vector of the interpolated state
      return std::vector<C_FLOAT64>(interpolated, interpolated + mData.dim);
    }
  }

  // should never happen 
  CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 32);
  return {};
}