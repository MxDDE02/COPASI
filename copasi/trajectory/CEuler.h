#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/math/CMathContainer.h"
#include <vector>
#include "copasi/utilities/CBrent.h"
class CEulerMethod : public CTrajectoryMethod
{
public:
  /**
   *  Construtor 
   */
  CEulerMethod(const CDataContainer * pParent,
               const CTaskEnum::Method & methodType,
               const CTaskEnum::Task & taskType);
/**
   * Copy-Construtor 
   */
  CEulerMethod(const CEulerMethod & src,
            const CDataContainer * pParent);

/**
   *  Deconstructor 
   * => deletes the used pointers 
   */
  ~CEulerMethod(); 

// == Functions ==

/**
   * 1. calls the start function of TrajectoryMethod to initialize the mContainerstate and time
   * 2. initilization of different parameters 
   */
  virtual void start();

/**
   * This methods must be called to elevate subgroups to
   * derived objects. The default implementation does nothing.
   * @return bool success
   */
  virtual bool elevateChildren();


/**
  * Check if the method is suitable for this problem
  * @return bool suitability of the method
  */
  virtual bool isValidProblem(const CCopasiProblem * pProblem);

/**
   *  This instructs the method to calculate a time step of deltaT
   *  starting with the current state, i.e., the result of the previous
   *  step.
   *  The new state (after deltaT) is expected in the current state.
   *  The return value is the actual timestep taken.
   *  @param const double & deltaT
   *  @param const bool & final (default: false)
   *  @return Status status
   */
  virtual CTrajectoryMethod::Status step(const double & deltaT, const bool & final); 



private:

/**
   *  Euler step (1 step)
   * -> updates the mathcontainer
   * -> also adapts the stepsize according to the local error 
   * -> checks for roots in between the interval of the one step 
   * @param startTime The time at which the euler step should start 
   */
  C_FLOAT64 doOneStep(C_FLOAT64 startTime);

/**
   * Function to initialize parameters.
   */
  void initializeParameter();

/**
   * Function to retrieve the roots from a certain interval (endTime - startTime)
   * returns the time of the first root of that innterval such as its state 
   * @param startTime startpoint of the root finder - start of the interval 
   * @param endTime endpoint of the root finder - end of the interval
   * @param t time of the root - return value 
   * @param f state of the root (should be 0) - return value 
   */
  //C_FLOAT64 findRoot(C_FLOAT64 startTime, C_FLOAT64 endTime, C_FLOAT64 &t, C_FLOAT64 &f);

/**
   * Checks if a root is located between the new calculated state and the previous calculated state 
   */
  bool checkRoots();


/**
    * Returns a linearly interpolated state for the given time.
    * @param t  The time at which the interpolated state is requested.
    * @return   A vector containing the interpolated state.
    */
   std::vector<C_FLOAT64> interpolateAttime(C_FLOAT64 t) const;

/**
   *  This evaluates the derivatives 
   * => uses the state vector to calculate a new rate vector 
   * @param t  The time at which the evaluation is supposed to occur 
   * @param y  state vector for the evauation (what should be put into the math container)
   * @param ydot evaluated rate - return value 
   */
  virtual void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);


  // == ROOT RELATED ARGUMENTS == 

  /**
   * An integer, referring number of roots
   */
  int mNumRoot;

   /**
   * Vectors for storing root values
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

// == INTEGRATOR RELATED ARGUMENTS
/**
   * The status of the integrator (failure, normal, root)
   */
  Status mStatus;
/**
   * Float, which stores the step size from the initialized parameter "Step size"
   */
  C_FLOAT64 mStepsize;

/**
   * value to record the outputtime for each interval
   */
  C_FLOAT64 outputTime;   

/**
   * value WHICH Stores the number of variables in the state vector (time + concentrations)
   */
  C_INT mdimension; 

/**
   * boolean expression, which determines beta at the adpative step size calculation 
   * if PI = false, beta will bet set zero -> no PI controller good for efficient adjustment
   * if PI = true, beta will be set 0.4/k -> PI controller for more stable stepsize adjustment 
   */
  bool PI; 
    /**
   * Pointer to the array with left hand side values of the math container (=state vector).
   * => time and concentration of the species 
   */
  CVector< C_FLOAT64 > mY;
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
  //C_FLOAT64 * mpYd;
  CVector< C_FLOAT64 > mYd;

  C_FLOAT64 * mpYd; 

  

// == INTERPOLATION == 

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
   * Pointer to the array of the interpolated state 
   */
  CVector< C_FLOAT64 > interpolated;
  C_FLOAT64 * pinterpolated; 


/**
   *  Struct to safe each intermediate step of the euler method with time and states (caution: the first entry of the state is the time)
   */
  struct TimeStatePair
  {
    C_FLOAT64 time;
    std::vector<C_FLOAT64> state;
  };

  /**
    * A history of time and state pairs stored during the integration.
    * This is used for linear interpolation to obtain the state at arbitrary times.
    */
  std::vector<TimeStatePair> mHistoryinter;

  C_FLOAT64 errorold; 
  CBrent::Eval * mpRootValueCalculator;
  C_FLOAT64 rootValue(const C_FLOAT64 & time); 

  void findRoot(C_FLOAT64 startTime, C_FLOAT64 endTime, C_FLOAT64 &t, C_FLOAT64 &f);

  std::vector<C_FLOAT64> rootvs(const C_FLOAT64 & time);

};