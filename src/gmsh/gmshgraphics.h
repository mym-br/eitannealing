/*
* gmshgraphics.h
*
*  Created on: Mar 16, 2017
*      Author: aksato
*/

#ifndef GMSHGRAPHICS_H_
#define GMSHGRAPHICS_H_

#include <QObject>
#include <memory>
#include "../problem.h"
#include "onelab.h"
#include<QModelIndex>

class gmshviewport : public QObject {
Q_OBJECT

public:
	gmshviewport(const char* name, std::string _outputFname, std::string _propertyName, const char* address, std::shared_ptr<problem> _input);
	~gmshviewport() {};

public slots:
	void solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight);

private:
	void saveTempMesh(double *sol);

	std::shared_ptr<problem> input;
	std::shared_ptr<onelab::remoteNetworkClient> gmeshClient;
	std::vector<onelab::string> ns;
	std::string outputFname, propertyName;

	int iter;
};

#endif // GMSHGRAPHICS_H_