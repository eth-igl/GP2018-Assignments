#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction {
public:

    RosenbrockFunction() {
		a = 1; b = 100;
    }

    virtual double computeValue(const VectorXd& x) {

		// Ex 1.1
		// return f(x)

	}

    virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1
		// write df/dx in `grad`

    }

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {

		// Ex 1.2
		// write d^2f/dx^2 in `hessianEntries`
    }

    double a, b;
};
