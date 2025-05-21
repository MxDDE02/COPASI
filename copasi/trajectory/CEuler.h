#include "copasi/trajectory/CTrajectoryMethod.h"
#include "copasi/math/CMathContainer.h"
#include <vector>

class CEulerMethod : public CTrajectoryMethod
{
public:
  CEulerMethod(const CDataContainer * pParent,
               const CTaskEnum::Method & methodType,
               const CTaskEnum::Task & taskType);

CEulerMethod(const CEulerMethod & src,
            const CDataContainer * pParent);

  virtual ~CEulerMethod();

  virtual void start();

  CTrajectoryMethod::Status step(const double & deltaT, const bool & final);

  C_FLOAT64 doSingleStep(C_FLOAT64 startTime, const C_FLOAT64 & endTime);

protected:
  void initializeParameter();

private:
  C_FLOAT64 mTargetTime;
  C_FLOAT64 mStepsize;

  // KEINE zusätzlichen Vektoren hier!
};