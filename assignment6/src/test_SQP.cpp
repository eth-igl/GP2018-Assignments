#include <iostream>

#include <FunctionConstraints.h>
#include <SQPFunctionMinimizer.h>

class Ex1Objective : public ObjectiveFunction {
public:
	// Ex 1
};

class Ex1Constraint : public FunctionConstraints {
public:
	// Ex 1
};

int main() {

	ObjectiveFunction* objective;
	FunctionConstraints* constraints;

	// Ex 1: uncomment these lines
//	objective = new Ex1Objective;
//	constraints = new Ex1Constraint;

	SQPFunctionMinimizer minimizer;
	VectorXd x(2); x << 100, -100;
	minimizer.minimize(objective, constraints, x);

	std::cout << "x            = " << x.transpose() << std::endl;
	std::cout << "f(x)         = " << objective->computeValue(x) << std::endl;
	std::cout << "c(x)         = " << constraints->getEqualityConstraintValues(x) << std::endl;
	std::cout << "# iterations = " << minimizer.lastNumberIterations << std::endl;
}
