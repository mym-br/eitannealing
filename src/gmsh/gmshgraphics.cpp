#include "gmshgraphics.h"
#include "../solution.h"
#include <QDir>

gmshviewport::gmshviewport(const char* name, std::string _outputFname, const char* address, std::shared_ptr<problem> _input) : input(_input), iter(0), outputFname(_outputFname) {
	gmeshClient = std::shared_ptr<onelab::remoteNetworkClient>(new onelab::remoteNetworkClient(name, address));
}

void gmshviewport::solution_updated(const QModelIndex & topLeft, const QModelIndex & bottomRight)
{
	{
		// Get solution pointer
		double *solution = (double *)topLeft.internalPointer() - topLeft.row();
		// Save solution file
		solution::saveMesh(solution, outputFname.c_str(), input, iter++);
		// Get working path
		QDir binDir(QDir::currentPath());
		// ask Gmsh to refresh
		onelab::string s("Gmsh/Action", "refresh");
		gmeshClient->set(s);
		if (iter % 10 == 0) gmeshClient->sendParseStringRequest("Delete View[0];");
		gmeshClient->sendMergeFileRequest(binDir.absoluteFilePath(outputFname.c_str()).toStdString());
	}
}