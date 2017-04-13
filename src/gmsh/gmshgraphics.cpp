#include "gmshgraphics.h"
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
		saveTempMesh(solution);
		// Get working path
		QDir binDir(QDir::currentPath());
		// ask Gmsh to refresh
		onelab::string s("Gmsh/Action", "refresh");
		gmeshClient->set(s);
		if (iter % 10 == 0) gmeshClient->sendParseStringRequest("Delete View[0];");
		gmeshClient->sendMergeFileRequest(binDir.absoluteFilePath(outputFname.c_str()).toStdString());
	}
}

void gmshviewport::saveTempMesh(double *sol) {
	std::ofstream myfile;
	myfile.open(outputFname.c_str());

	std::ifstream inputfile(input->getMeshFilename());
	for (int i = 0; inputfile.eof() != true; i++) {
		std::string line;
		std::getline(inputfile, line);
		myfile << line << '\n';
	}

	//Salvando os tensoes nos no's em formato para ser utilizado no gmsh
	myfile << "$NodeData\n1\n\"electric potential\"\n1\n0.0\n3\n" << iter++ << "\n1\n" << input->getNodesCount() << "\n";
	for (int j = 0; j < input->getNodesCount(); j++) {
		myfile << (j + 1) << "\t" << sol[input->getNode2Coefficient(j)] << "\n";
	}
	myfile << "$EndNodeData\n";
	myfile.flush();
	myfile.close();
}