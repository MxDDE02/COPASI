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
mdimension(), 
errorold(0),
mY(), 
mpYdot(NULL),
mYd(), 
mpY(NULL),
mpYd(NULL),
interpolated(),
pinterpolated(NULL), 
mNumRoot(0),
mRootsA(),
mRootsB(),
mRootsNonZero(),
mpRootValueOld(NULL),
mpRootValueNew(NULL),
mLastRootTime(-std::numeric_limits< C_FLOAT64 >::infinity()), 
mpRootValueCalculator(NULL)

{
  initializeParameter();
}

CEulerMethod::CEulerMethod(const CEulerMethod & src,
                           const CDataContainer * pParent): 
CTrajectoryMethod(src, pParent),
mdimension(), 
errorold(0),
mY(), 
mpYdot(NULL), 
mYd(), 
mpY(NULL),
mpYd(NULL),
interpolated(), 
pinterpolated(NULL), 
mNumRoot(src.mNumRoot),
mRootsA(src.mRootsA),
mRootsB(src.mRootsB),
mRootsNonZero(src.mRootsNonZero),
mpRootValueOld(NULL),
mpRootValueNew(NULL), 
mpRootValueCalculator(NULL)
{
  initializeParameter();
}

CEulerMethod::~CEulerMethod()
{
  //pdeletev(mpYd);
  //pdeletev(mpY);
  //pdeletev(pinterpolated); 
  if (mRootsFound.array() != NULL)
    {
      delete [] mRootsFound.array();
    }
 }

void CEulerMethod::initializeParameter()
{
  assertParameter("initial step size", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.000001);
  assertParameter("absolute tolerance", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.000000001);
  assertParameter("relative tolerance", CCopasiParameter::Type::DOUBLE, (C_FLOAT64) 0.000001);
  assertParameter("Maximal internal steps", CCopasiParameter::Type::INT, 1000000);
  assertParameter("PI controller for adaptive stepsize", CCopasiParameter::Type::BOOL, true);
  mpRootValueCalculator = new CBrent::EvalTemplate< CEulerMethod >(this, &CEulerMethod::rootValue);

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

  // 3. Determine the mdimension of the system (excluding fixed event targets)
  mdimension = mContainerState.size() - mpContainer->getCountFixedEventTargets(); 

  // 4. Get pointer to rate vector (excluding fixed event targets)
  mpYdot = mpContainer->getRate(*mpReducedModel).array() + mpContainer->getCountFixedEventTargets();

  // 6. Retrieve input parameters from the Trajectory Problem 
  mStepsize = getValue< double >("initial step size");
  euler_atolerance = getValue< double >("absolute tolerance");
  euler_rtolerance = getValue< double >("relative tolerance");
  steplimit = getValue< int >("Maximal internal steps");
  PI = getValue< bool >("PI controller for adaptive stepsize");

  // 7. Allocate memory 
  //mpY = new C_FLOAT64[mdimension];
  //mpYd = new C_FLOAT64[mdimension];
  //pinterpolated = new C_FLOAT64[mdimension];
  interpolated.resize(mdimension); 
  pinterpolated = interpolated.array(); 
  mYd.resize(mdimension);
  mpYd = mYd.array();
  mY.resize(mdimension);
  mpY = mY.array();

  // set outputime (to zero)
  outputTime = *mpContainerStateTime;

  // 8. Copy current state into local state vector
  memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));

  //========Initialize Roots Related Arguments========
  mNumRoot = mpContainer->getRoots().size();

  if (mRootsFound.array() != NULL)
    {
      delete [] mRootsFound.array();
    }

  mRootsFound.initialize(mNumRoot, new C_INT[mNumRoot]);
  mRootsA.resize(mNumRoot);
  mRootsB.resize(mNumRoot);
  mpRootValueNew = &mRootsA;
  mpRootValueOld = &mRootsB;
  mRootsNonZero.resize(mNumRoot);
  mRootsNonZero = 0.0;
  mLastRootTime = -std::numeric_limits< C_FLOAT64 >::infinity();
  *mpRootValueOld = mpContainer->getRoots();

  errorold = errorold = -std::numeric_limits< C_FLOAT64 >::infinity();
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
  CVector< C_FLOAT64 > yTemp(mdimension);
  memcpy(yTemp.array(), mpContainerStateTime, mdimension * sizeof(C_FLOAT64));

  if (y != mpContainerStateTime)
    memcpy(mpContainerStateTime, y, mdimension * sizeof(C_FLOAT64));
  //this does the actual evaluation and puts it into ydot 
  //mpContainer->updateSimulatedValues(*mpReducedModel);
  mpContainer->updateSimulatedValues(false);
  memcpy(ydot, mpYdot, mdimension * sizeof(C_FLOAT64));
//this is old debugging 
#ifdef DEBUG_NUMERICS
  std::cout << "State:     " << mpContainer->getState(false) << std::endl;
  std::cout << "Rate:      " << mpContainer->getRate(false) << std::endl;
#endif // DEBUG_NUMERICS
//retreave the old value 
  memcpy(mpContainerStateTime, yTemp.array(), mdimension * sizeof(C_FLOAT64));

  return;
}


CTrajectoryMethod::Status CEulerMethod::step(const double & deltaT, const bool & /* final */)
{
  outputTime = *mpContainerStateTime + deltaT;
  memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));
  mpContainer->updateSimulatedValues(false);
  mpContainer->updateRootValues(false);
  *mpRootValueOld = mpContainer->getRoots();
  int internalsteps = 0; 

  //Saving of the initial Math-Container in the mHistory - for interpolation
  TimeStatePair initial;
  initial.time = *mpContainerStateTime;
  initial.state.assign(mpY, mpY + mdimension);
  mHistoryinter.push_back(initial);

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
    if (mStatus == ROOT)//||
          //(mNumRoot > 0 && checkRoots()))
        {
          //clearing the mHistory for next steps 
          mHistoryinter.clear();
          return ROOT;
        }
  } 
  
  //Interpolation -> get the Trajectory Problem defined state at the requested time 
  std::vector<C_FLOAT64> pinterpolatedState = interpolateAttime(outputTime);
  //Update everything to the pinterpolated output 
  memcpy(mpContainerStateTime, pinterpolatedState.data(), mdimension * sizeof(C_FLOAT64));
  memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));
  *mpContainerStateTime = outputTime;
  //updating the state and the roots -> for root finding 
  mpContainer->updateSimulatedValues(false);
  mpContainer->updateRootValues(false);
  
  //clearing the mHistory for next steps 
  mHistoryinter.clear();
  return NORMAL;
}


C_FLOAT64 CEulerMethod::doOneStep(C_FLOAT64 startTime)
{
  //precaution of the step
  mpContainer->updateSimulatedValues(false);
  mpContainer->updateRootValues(false);
  *mpRootValueOld = mpContainer->getRoots();
  bool stepreject = false; 
  bool accepted = false;

  //step calculation 
  while (!accepted)
  {
    C_FLOAT64 current_stepsize = mStepsize;
    //making copies of the original state + vector allocation
    std::vector<C_FLOAT64> y_original(mpY, mpY + mdimension);
    C_FLOAT64 t_old = *mpContainerStateTime;

    std::vector<C_FLOAT64> fullstep(mdimension);
    std::vector<C_FLOAT64> halfstep(mdimension);
    std::vector<C_FLOAT64> deltaerror(mdimension);
    std::vector<C_FLOAT64> scale(mdimension);

    // == 1. calculate rates ==
    evalF(mpContainerStateTime, mpY, mpYd);
    

    // == 2. full step ==
    for (int i = 0; i < mdimension; ++i)
      fullstep[i] = y_original[i] + current_stepsize * mpYd[i];

    // == 3. first half step ==
    for (int i = 0; i < mdimension; ++i)
      mpY[i] = y_original[i] + (current_stepsize*0.5) * mpYd[i];

    *mpContainerStateTime = t_old + (current_stepsize*0.5);

    evalF(mpContainerStateTime, mpY, mpYd);
    // memcpy(mpContainerStateTime, mpY, mdimension * sizeof(C_FLOAT64));
    // mpContainer->updateSimulatedValues(false);
    // memcpy(yd_temp.data(), mpYdot, mdimension * sizeof(C_FLOAT64));

    // == 4. second half step ==
    for (int i = 0; i < mdimension; ++i)
      halfstep[i] = mpY[i] + (current_stepsize*0.5) * mpYd[i];
    
    //*mpContainerStateTime = t_old + mStepsize / 2.0;

    // == 5. error estimation ==
    C_FLOAT64 localerror = 0.0;
    C_FLOAT64 hscale; 
    C_FLOAT64 k =2; 
    C_FLOAT64 beta; 
    if(PI==true)
    {
      beta = 0.4/k;
    }
    else
    {
      beta = 0.0;
    }
    C_FLOAT64 alpha = (1/k) - (0.75 * beta); 
    C_FLOAT64 safe = 0.9; 
    C_FLOAT64 minhscale = 0.2; 
    C_FLOAT64 maxhscale = 10.0; 
    C_FLOAT64 sum_sq = 0.0;
    C_FLOAT64 term = 0.0;
    if(mdimension>1)
    {
    for (int i = 1; i < mdimension; ++i)
    {
      deltaerror[i] = std::abs(halfstep[i] - fullstep[i]);
      scale[i] = euler_atolerance + std::max(std::abs(y_original[i]), std::abs(fullstep[i])) * euler_rtolerance;
      term = deltaerror[i] / scale[i];
      sum_sq += term * term;
    }
    localerror = sqrt(sum_sq / (mdimension-1));
    }
    else
    {
    for (int i = 0; i < mdimension; ++i)
    {
      deltaerror[i] = std::abs(halfstep[i] - fullstep[i]);
      scale[i] = euler_atolerance + std::max(std::abs(y_original[i]), std::abs(fullstep[i])) * euler_rtolerance;
      term = deltaerror[i] / scale[i];
      sum_sq += term * term;
    }
    localerror = sqrt(sum_sq / mdimension);
    }


    // == 6. decition if step can be accepted ==
  if (localerror <= 1.0)
    {
      for (int i = 0; i < mdimension; ++i)
        mpY[i] = fullstep[i];

      *mpContainerStateTime = t_old + current_stepsize;

      //this updates the math container to the new state 
      //*mpContainerStateTime = t_old + mStepsize;
      memcpy(mpContainerStateTime, mpY, mdimension * sizeof(C_FLOAT64));
      mpContainer->updateSimulatedValues(false);
      mpContainer->updateRootValues(false);

      //saving the intermediate state for interpolation
      TimeStatePair ts;
      ts.time = *mpContainerStateTime;
      ts.state.assign(mpY, mpY + mdimension);
      mHistoryinter.push_back(ts);

      accepted = true;

      if(localerror == 0.0)
      {
        hscale = 1; 
      }
      else
      {
        hscale = safe * pow(localerror,-alpha)*pow(errorold, beta);
        if(hscale<minhscale){
          hscale = minhscale;
        }
        if(hscale>maxhscale){
          hscale = maxhscale; 
        }
      }

      if(stepreject) //previous step was not accepted
      {
        hscale = std::max(hscale, 1.0); 
        mStepsize *= hscale;
      }
        mStepsize *= hscale; 
        errorold = std::max(localerror, 1.0e-04); 
        stepreject =false; 

      if(*mpContainerStateTime > outputTime)
      {
        std::vector<C_FLOAT64> outputstate = interpolateAttime(outputTime);
        memcpy(mpContainerStateTime, outputstate.data(), mdimension * sizeof(C_FLOAT64));
        memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));
        *mpContainerStateTime = outputTime;
        mpContainer->updateSimulatedValues(false); 
        mpContainer->updateRootValues(false);
      }
    
    C_FLOAT64 Tolerance = 100.0 * (fabs(*mpContainerStateTime) * std::numeric_limits< C_FLOAT64 >::epsilon() + std::numeric_limits<C_FLOAT64>::min());

      // == EVENTS == 
    if(mNumRoot>0)
    {
      if (checkRoots())
      {
        // C_FLOAT64 t; 
        // C_FLOAT64 f;
        // findRoot(startTime, *mpContainerStateTime, t, f); 
        // C_FLOAT64 RootValue = f; 
        // C_FLOAT64 RootTime = t; 
        C_FLOAT64 RootTime; 
        C_FLOAT64 RootValue; 
        CBrent::findRoot(startTime, *mpContainerStateTime, mpRootValueCalculator, &RootTime, &RootValue, 1e-9);

        Tolerance = 100.0 * (fabs(*mpContainerStateTime) * std::numeric_limits< C_FLOAT64 >::epsilon() + std::numeric_limits< C_FLOAT64 >::min());
        //Precaution if the Root is not the wanted root 
        if (RootTime > outputTime)
        {
          mStatus = NORMAL;
        }

        else if (fabs(RootTime -mLastRootTime) < Tolerance)
        {
          mStatus = NORMAL;
        }
        //now get back to the Event time 
        else if (mLastRootTime < RootTime)
        {
          mLastRootTime = RootTime;
          // if(CBrent::findRoot(startTime, RootTime, mpRootValueCalculator, &RootTime, &RootValue, 1e-9))
          // {
          // CCopasiMessage(CCopasiMessage::ERROR, "Discontinous events");
          // }
          //Update everything to the pinterpolated output 
          std::vector<C_FLOAT64> pinterpolatedRootstate = interpolateAttime(RootTime);
          memcpy(mpContainerStateTime, pinterpolatedRootstate.data(), mdimension * sizeof(C_FLOAT64));
          memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));
          *mpContainerStateTime = RootTime;
          mpContainer->updateSimulatedValues(false); 
          mpContainer->updateRootValues(false);
          *mpRootValueNew = mpContainer->getRoots();

          // Mark the appropriate root
          C_INT * pRootFound = mRootsFound.array();
          C_INT * pRootFoundEnd = pRootFound + mNumRoot;
          C_FLOAT64 * pRootValue = mpRootValueNew->array();

          Tolerance = 100.0 * (fabs(*mpContainerStateTime) * std::numeric_limits< C_FLOAT64 >::epsilon() + std::numeric_limits< C_FLOAT64 >::min());
    
          for (; pRootFound != pRootFoundEnd; ++pRootFound, ++pRootValue)
          {
          // Added a numerical Tolerance just to make sure 
          //if (*pRootValue == RootValue || *pRootValue == -RootValue)
            if (std::fabs(*pRootValue - RootValue) < Tolerance || std::fabs(*pRootValue + RootValue) < Tolerance || std::fabs(*pRootValue-0) < Tolerance )
              {
                *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleBoth);
              }
            else
              {
                *pRootFound = static_cast< C_INT >(CMath::RootToggleType::NoToggle);
              }
          }
        
          //most important thing to fire Events 
          mStatus = ROOT;
          return RootTime - startTime;
        }
      }
    else 
      {
        mStatus = NORMAL; 
        return *mpContainerStateTime;
        mHistoryinter.clear();
      }
    }
  }
  
  else
    {
      // step is not accepted -> stepsize will get reduced 
      //mStepsize = mStepsize * std::sqrt(euler_rtolerance / localerror); 
      //mStepsize = mStepsize * std::sqrt(0.9 / localerror); 
      hscale = std::max(safe*pow(localerror, -alpha), minhscale);
      mStepsize *= hscale;
      //bring back original state
      memcpy(mpY, y_original.data(), mdimension * sizeof(C_FLOAT64));
      *mpContainerStateTime = t_old;
      stepreject = true;
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
  C_FLOAT64 Tolerance = 100.0 * (fabs(*mpContainerStateTime) * std::numeric_limits< C_FLOAT64 >::epsilon() + std::numeric_limits< C_FLOAT64 >::min());
  if((t - mHistoryinter.front().time)<Tolerance)
  {
    for (int i = 0; i < mHistoryinter.size(); ++i)
    {
      const TimeStatePair& p0 = mHistoryinter[i];
      const TimeStatePair& p1 = mHistoryinter[i+1];
      //linear interpolation 
      for (int j = 0; j < mdimension; ++j)
      {
        pinterpolated[j] = p0.state[j]+ ((t - p0.time) / (p1.time - p0.time)) * (p1.state[j] - p0.state[j]);
      }

        // retun: the vector of the pinterpolated state
        return std::vector<C_FLOAT64>(pinterpolated, pinterpolated + mdimension);
      }
    }
  
  else
  {
  //actual interpolation
  for (int i = 0; i < mHistoryinter.size(); ++i)
  {
    C_FLOAT64 timecheck = mHistoryinter[i].time;
    if (timecheck >= t)
    {
      const TimeStatePair& p0 = mHistoryinter[i - 1];
      const TimeStatePair& p1 = mHistoryinter[i];
      //linear interpolation 
      for (int j = 0; j < mdimension; ++j)
      {
        pinterpolated[j] = p0.state[j] + ((t - p0.time) / (p1.time - p0.time)) * (p1.state[j] - p0.state[j]);
      }

      // retun: the vector of the pinterpolated state
      return std::vector<C_FLOAT64>(pinterpolated, pinterpolated + mdimension);
    }
  }
  }
  // should never happen 
  CCopasiMessage(CCopasiMessage::ERROR, MCTrajectoryMethod + 32);
  return {};
}


bool CEulerMethod::checkRoots()
{
  bool hasRoots = false;

  // Swap old and new root values
  // CVector< C_FLOAT64 > * pTmp = mpRootValueOld;
  // mpRootValueOld = mpRootValueNew;
  // mpRootValueNew = pTmp;


  *mpRootValueNew = mpContainer->getRoots();

  //declare the important pointers 
  C_FLOAT64 *pRootValueOld = mpRootValueOld->array();
  C_FLOAT64 *pRootValueNew = mpRootValueNew->array();
  C_FLOAT64 *pRootNonZero = mRootsNonZero.array();

  C_INT *pRootFound = mRootsFound.array();
  C_INT *pRootFoundEnd = pRootFound + mRootsFound.size();

  //for every root - compare the old and new values 
  for (; pRootFound != pRootFoundEnd; 
       pRootValueOld++, pRootValueNew++, pRootFound++, pRootNonZero++)
  {
    if (*pRootValueOld * *pRootValueNew < 0.0 || //sign must have changed 
        (*pRootValueNew == 0.0 && *pRootValueOld != 0.0)) //if new value is exactly 0 -> now event 
    {
      hasRoots = true;
      *pRootFound = static_cast<C_INT>(CMath::RootToggleType::ToggleBoth);
      *pRootNonZero = *pRootValueOld;
    }
    else if (*pRootValueNew == 0.0 &&
               *pRootValueOld != 0.0)
        {
          hasRoots = true;
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleEquality); // toggle only equality
          *pRootNonZero = *pRootValueOld;
        }
      else if (*pRootValueNew != 0.0 &&
               *pRootValueOld == 0.0 &&
               *pRootValueNew **pRootNonZero < 0.0)
        {
          hasRoots = true;
          *pRootFound = static_cast< C_INT >(CMath::RootToggleType::ToggleInequality); // toggle only inequality
        }
    else
    {
      *pRootFound = static_cast<C_INT>(CMath::RootToggleType::NoToggle);
    }
  }
  return hasRoots;
}

void CEulerMethod::findRoot(C_FLOAT64 startTime, C_FLOAT64 endTime, C_FLOAT64 &t, C_FLOAT64 &f)
{
  //initlializing the things we need 
  C_FLOAT64 *pRootValueOld = mpRootValueOld->array();
  C_FLOAT64 *pRootValueNew = mpRootValueNew->array();
  C_FLOAT64 oldtime = startTime;
  C_FLOAT64 newtime = endTime;
  CVector<C_FLOAT64> time (mNumRoot);  
  CVector<C_FLOAT64> oldroot (mNumRoot);  
  CVector<C_FLOAT64> newroot (mNumRoot);  
  time =std::numeric_limits<C_FLOAT64>::infinity();
  oldroot = std::numeric_limits<C_FLOAT64>::infinity();
  newroot =std::numeric_limits<C_FLOAT64>::infinity();

  //calculating the time of the event (root = 0)
  for (int i = 0; i < mNumRoot; ++i)
  {
    C_FLOAT64 fOld = pRootValueOld[i];
    C_FLOAT64 fNew = pRootValueNew[i];
    if(fNew*fOld <= 0) // only for those, which have a sign change 
    {
      C_FLOAT64 t2 = oldtime - fOld * ((newtime - oldtime) / (fNew - fOld));
      //we have to make sure the time is the one we want - time dependent roots are annoying otherways 
      if (t2 > mLastRootTime) // just precaution - should not happen becuase we checked for sign change
      time[i] = t2;
      //next to the time we are saving the root state values for f calculation
      oldroot[i] = fOld;  
      newroot[i] = fNew;  
    }
  }
  //Now retrieving the time of the first event happening in the interval 
  C_FLOAT64 *minIt = std::min_element(time.begin(), time.end());
  int rootIndex = std::distance(time.begin(), minIt);
  t = *minIt;
  //we also would like the root state at that time (should be ~0)
  C_FLOAT64 fOldAtRoot = oldroot[rootIndex];
  C_FLOAT64 fNewAtRoot = newroot[rootIndex];
  // via linear interpolation
  f = fOldAtRoot + ((t - oldtime) / (newtime - oldtime)) * (fNewAtRoot - fOldAtRoot);
}

C_FLOAT64 CEulerMethod::rootValue(const C_FLOAT64 & time)
{
  std::vector<C_FLOAT64> interstate = interpolateAttime(time);
  memcpy(mpContainerStateTime, interstate.data(), mdimension * sizeof(C_FLOAT64));
  memcpy(mpY, mpContainerStateTime, mdimension * sizeof(C_FLOAT64));
  *mpContainerStateTime = time;
  mpContainer->updateSimulatedValues(false); 
  mpContainer->updateRootValues(false);
  *mpRootValueNew = mpContainer->getRoots();

  // *mpContainerStateTime = time;
  // mpContainer->applyUpdateSequence(mUpdateTimeDependentRoots);

  const C_FLOAT64 * pRoot = mpContainer->getRoots().array();
  const C_FLOAT64 * pRootEnd = pRoot + mNumRoot;
  const C_FLOAT64 * pRootOld = mpRootValueOld->array();
  const C_FLOAT64 * pRootNew = mpRootValueNew->array();

  C_FLOAT64 MaxRootValue = - std::numeric_limits< C_FLOAT64 >::infinity();
  C_FLOAT64 RootValue;

  for (; pRoot != pRootEnd; ++pRoot, ++pRootOld, ++pRootNew)
    {
      // We are only looking for roots which change sign in [pOld, pNew]
      if (*pRootOld **pRootNew < 0 || *pRootNew == 0)
        {
          // Assure that the RootValue is increasing between old and new for each
          // candidate root.
          RootValue = (*pRootNew >= *pRootOld) ? *pRoot : -*pRoot;

          if (RootValue > MaxRootValue)
            {
              MaxRootValue = RootValue;
            }
        }
    }
  return MaxRootValue;
}

