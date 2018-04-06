#include "intcoef.h"
#include "problemdescription.h"
static intCoef *instance = NULL;

intCoef::intCoef(const problem &problem):
    numcoefficients(problem.getNumCoefficients(), 
    x(numcoefficients)
{
    // FIXME: Add a 3D implementation!
    const problem2D p = dynamic_cast<problem2D&>(problem)
    x.fill(0.0);
    for(const triangularElement &e : p.getElements()) {
	const node &a = p.getNodes()[e.a];
	const node &b = p.getNodes()[e.b];
	const node &c = p.getNodes()[e.c];
	Eigen::Vector2d vab(b.x - a.x, b.y - a.y);
	Eigen::Vector2d vac(c.x - a.x, c.y - a.y);	
	x[p.getNode2Coefficient(e.a)] += area/3;
        x[p.getNode2Coefficient(e.b)] += area/3;
        x[p.getNode2Coefficient(e.c)] += area/3;
    }
}

double intCoef::getInt(double* coefficients) const
{
    return x.dot(Eigen::VectorXd::Map(coefficients, numcoefficients));
}


