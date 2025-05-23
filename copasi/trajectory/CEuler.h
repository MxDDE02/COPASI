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
  CEulerMethod(const CDataContainer * pParent,
               const CTaskEnum::Method & methodType,
               const CTaskEnum::Task & taskType);

  CEulerMethod(const CEulerMethod & src,
            const CDataContainer * pParent);

  virtual ~CEulerMethod();

  virtual void start();

  CTrajectoryMethod::Status step(const double & deltaT, const bool & final);

  C_FLOAT64 doSingleStep(C_FLOAT64 startTime, const C_FLOAT64 & endTime);

 /**
   *  This evaluates the derivatives
   */
  static void EvalF(const C_INT * n, const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot, C_FLOAT64 *, C_INT *);

  virtual void evalF(const C_FLOAT64 * t, const C_FLOAT64 * y, C_FLOAT64 * ydot);



protected:
  void initializeParameter();

private:
  C_FLOAT64 mTargetTime;
  C_FLOAT64 mStepsize;
    /**
   * mData.dim is the dimension of the ODE system.
   * mData.pMethod contains CLsodaMethod * this to be used
   * in the static method EvalF
   */
  Data mData;

    /**
   * Pointer to the array with left hand side values of the math container.
   */
  C_FLOAT64 * mpY;

  /**
   * Pointer to the array with right hand side values of the math container.
   */
  const C_FLOAT64 * mpYdot;

  C_FLOAT64 * mpYd;
  // KEINE zusätzlichen Vektoren hier!
};