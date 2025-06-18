#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/math/CMathContainer.h"
#include <vector>

class CEulerMethod : public CTrajectoryMethod
{
public:
  struct Data
  {
    size_t dim;
    CEulerMethod * pMethod;
  };
  /**
   *  Struct to safe each intermediate step of the euler method
   */
  struct TimeStatePair
  {
    C_FLOAT64 time;
    std::vector<C_FLOAT64> state;
  };
  /**
   *  Construtor 
   */
  CEulerMethod(const CDataContainer * pParent,
               const CTaskEnum::Method & methodType,
               const CTaskEnum::Task & taskType);
/**
   *  Copy-Construtor 
   */
  CEulerMethod(const CEulerMethod & src,
            const CDataContainer * pParent);

/**
   *  Deconstructor 
   * => deletes the used pointers (mpYd and mpY)
   */
  ~CEulerMethod(); 

  /**
   * This methods must be called to elevate subgroups to
   * derived objects. The default implementation does nothing.
   * @return bool success
   */
  virtual bool elevateChildren();

/**
   * 1. calls the start function of TrajectoryMethod to initialize the mContainerstate and time
   * 2. initializing of different parameters: 
   * - mData.dim (dimension of the ode problem)
   * - mpYdot (array of rates)
   * - mTargetTime to get the value of the duration (from TrajectoryProblem)
   * - mStepsize to get the value of the step sizes given in the initializeParamter-function
   * - mpY and mpYd as an empy pointer array of size mData.dim  
   * - copies the state vector into mpY 
   * - Saving of the initial Math-Container in the mHistory - for interpolation
   */
  virtual void start();
/**
   *  Does the Euler integration until the first time interval is reached (Trajectory Problem)
   *  - does euler steps
   *  - interpolates
   */
  virtual CTrajectoryMethod::Status step(const double & deltaT, const bool & final); //added virtual to that 

  /**
   *  Euler step (1 step)
   * -> updates the mathcontainer
   * -> also adapts the stepsize according to the local error 
   * @param t  The time at which the euler step should start 
   */
  bool doOneStep(C_FLOAT64 startTime);

   /**
  * Check if the method is suitable for this problem
  * @return bool suitability of the method
  */
  virtual bool isValidProblem(const CCopasiProblem * pProblem);




  /**
    * Returns a linearly interpolated state for the given time.
    *
    * @param t  The time at which the interpolated state is requested.
    * @return   A vector containing the interpolated state.
    */
   std::vector<C_FLOAT64> interpolateAttime(C_FLOAT64 t) const;

/**
   *  This evaluates the derivatives 
   * => uses the state vector to calculate a new rate vector 
   */
  virtual void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

  //virtual void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);

/**
   *  This evaluates the events 
   */
  virtual void evalR(const C_FLOAT64 * t, const C_FLOAT64 * y, const C_INT * nr, C_FLOAT64 * r);




protected:
  /**
   * Function to initialize parameters. 
   * Just determines the step size used for the Euler Method. 
   * The default is 0.01.
   */
  void initializeParameter();

  /**
   * The status of the integrator
   */
  Status mStatus;

  /**
   * An integer, referring number of roots
   */
  size_t mNumRoot;

   /**
   * 2 Vector for storing root value
   */
  CVector< C_FLOAT64 > mRootsA;
  CVector< C_FLOAT64 > mRootsB;
  CVector< C_FLOAT64 > mRootsNonZero;

  /**
   * Pointer to the vector holding the previously calculated roots
   */
  CVector< C_FLOAT64 > *mpRootValueOld;

  /**
   * Pointer to the vector holding the newly calculated roots
   */
  CVector< C_FLOAT64 > *mpRootValueNew;

  /**
   * The last time dependent root time
   */
  C_FLOAT64 mLastRootTime;

private:

  /**
   * Check whether finds a root
   */
  bool checkRoots();
/**
   * Float, which stores the step size from the initialized parameter "Step size"
   */
  C_FLOAT64 mStepsize;

/**
   * Float, which stores the absolute error tolerance from the initialized parameter "absolute tolerance"
   */
  C_FLOAT64 euler_atolerance;

/**
   * Float, which stores the relative error tolerance from the initialized parameter "relative tolerance"
   */
  C_FLOAT64 euler_rtolerance;

/**
   * value to record the predefined maximal internal steps - used in step 
   */
  int steplimit; 
  
  /**
   * mData.dim is the dimension of the ODE system. 
   * It is used to determine the size of the state vector
   */
  Data mData;

    /**
   * Pointer to the array with left hand side values of the math container (=state vector).
   * => time and concentration of the species 
   */
  C_FLOAT64 * mpY;

  /**
   * Pointer to the array with right hand side values of the math container (=rate vector).
   * => only used for the evalF function 
   */
  const C_FLOAT64 * mpYdot;

  /**
   * Pointer to the array with the rates of the reactions evaluated (=rate vector - intermediate)
   * => used to determine the rate of each equation with a given mpY - more used as temporary storage
   * 
   */
  C_FLOAT64 * mpYd;

  /**
   * Pointer to the array of the interpolated state 
   */

  C_FLOAT64 * interpolated; 
/**
   * value to record the outputtime for each interval
   */
  C_FLOAT64 outputTime; 
  /**
    * A history of time and state pairs stored during the integration.
    * This is used for linear interpolation to obtain the state at arbitrary times.
    */
  std::vector<TimeStatePair> mHistoryinter;
  std::vector<TimeStatePair> mHistoryrootsstart;
  std::vector<TimeStatePair> mHistoryrootsinterpolate;
};