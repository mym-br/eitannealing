#include "intcoef.h"
#include "problem.h"

#include "twodim/problem2D.h"

static intCoef *instance = NULL;

intCoef::intCoef(problem &problem):
    numcoefficients(problem.getNumCoefficients()), 
    x(numcoefficients)
{
    // FIXME: Add a 3D implementation!
	try {
		problem2D &p = dynamic_cast<problem2D&>(problem);
		x.fill(0.0);
		for (const problem2D::triangularElement &e : p.getElements()) {
			const problem2D::node &a = p.getNodes()[e.a];
			const problem2D::node &b = p.getNodes()[e.b];
			const problem2D::node &c = p.getNodes()[e.c];
			Eigen::Vector2d vab(b.x - a.x, b.y - a.y);
			Eigen::Vector2d vac(c.x - a.x, c.y - a.y);
			double area = std::abs(vab.x()*vac.y() - vab.y()*vac.x()) / 2;
			x[p.getNode2Coefficient(e.a)] += area / 3;
			x[p.getNode2Coefficient(e.b)] += area / 3;
			x[p.getNode2Coefficient(e.c)] += area / 3;
		}
	}
	catch (const std::bad_cast&) {
		// Cast failed for 3d case (temporary solution)
		x.fill(0.0);
	}
}

double intCoef::getInt(double* coefficients) const
{
    return x.dot(Eigen::VectorXd::Map(coefficients, numcoefficients));
}


