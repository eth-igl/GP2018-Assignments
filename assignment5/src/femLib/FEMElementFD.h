#pragma once

#include "FEMElement.h"

/**
	This class implements Constant Strain Triangles elements in 2D
*/
class FEMElementFD : public FEMElement {
public:
	FEMElementFD(const std::array<int, 3> &nodeIndices, const VectorXd &X)
		: FEMElement(nodeIndices, X) {
	}
	~FEMElementFD(){
	}

	virtual double getEnergy(const VectorXd& x, const VectorXd& X)
	{
		// Ex 2.1
		return 0;
	}
	virtual void addEnergyGradientTo(const VectorXd& x, const VectorXd& X, VectorXd& grad)
	{
		// Ex 2.1
	}
	virtual void addEnergyHessianTo(const VectorXd& x, const VectorXd& X, std::vector<Tripletd>& hesEntries)
	{
		// Ex 2.2
	}

private:
	double h = 1e-8; // step size
};
