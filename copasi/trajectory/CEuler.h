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
   *  - does euler steps (with the predefinded step size)
   *  - interpolates, if the step size is bigger then the intervall size 
   */
  virtual CTrajectoryMethod::Status step(const double & deltaT, const bool & final); //added virtual to that 

  /**
   *  Euler step (1 step)
   * -> updates the mathcontainer
   * @param t  The time at which the euler step should start 
   */
  bool doOneStep(C_FLOAT64 startTime);



  /**
    * Returns a linearly interpolated state for the given time.
    *
    * @param t  The time at which the interpolated state is requested.
    * @return   A vector containing the interpolated state.
    */
   std::vector<C_FLOAT64> interpolateAt(C_FLOAT64 t) const;

 /**
   *  Function to calculate the error (2nd order) 
   * @param t  The time at which the error calculation is required 
   * - utilized step doubling to calculate the local error 
   * - stepsize is adjusted according to the error 
   */
   C_FLOAT64 estimateError(C_FLOAT64 t);


 /**
   *  This evaluates the derivatives 
   * => uses the state vector to calculate a new rate vector 
   */
  //virtual void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);



protected:
  /**
   * Function to initialize parameters. 
   * Just determines the step size used for the Euler Method. 
   * The default is 0.01.
   */
  void initializeParameter();

private:
/**
   * Float, which stores the step size from the initialized parameter "Step size"
   */
  C_FLOAT64 mStepsize;
/**
   * Float, which stores the epsilon from the initialized parameter "epsilon" - error tolerance
   */
  C_FLOAT64 euler_epsilon;

/**
   * Float, which stores the absolute error tolerance from the initialized parameter "absolute tolerance"
   */
  C_FLOAT64 euler_atolerance;

/**
   * Float, which stores the relative error tolerance from the initialized parameter "relative tolerance"
   */
  C_FLOAT64 euler_rtolerance;
  
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
    * A history of time and state pairs stored during the integration.
    * This is used for linear interpolation to obtain the state at arbitrary times.
    */
  std::vector<TimeStatePair> mHistoryinter;
};