#include "intcoef.h"

intCoef::intCoef(
      const std::vector<node> &nodes, 	       
      const std::vector<triangularElement> &elements,
      const std::map<int, int> &node2coefficient,
      int numcoefficients,
      float totalheight):numcoefficients(numcoefficients), x(numcoefficients)
{
    for(const triangularElement &e : elements) {
	const node &a = nodes[e.a];
	const node &b = nodes[e.b];
	const node &c = nodes[e.c];
	Eigen::Vector2d vab(b.x - a.x, b.y - a.y);
	Eigen::Vector2d vac(c.x - a.x, c.y - a.y);
	
	double area = fabs(vab.x()*vac.y()-vab.y()*vac.x())/2;
	std::map<int, int>::const_iterator i;
	if((i=node2coefficient.find(e.a))!=node2coefficient.end())
	  x.coeffRef(i->second) += totalheight*area/3;
	if((i=node2coefficient.find(e.b))!=node2coefficient.end())
	  x.coeffRef(i->second) += totalheight*area/3;
	if((i=node2coefficient.find(e.c))!=node2coefficient.end())
	  x.coeffRef(i->second) += totalheight*area/3;
    }
}

double intCoef::getInt(double* coefficients) const
{
    return x.dot(Eigen::VectorXd::Map(coefficients, numcoefficients));
}
