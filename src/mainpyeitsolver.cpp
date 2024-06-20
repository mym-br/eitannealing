#include <iostream>
#include <fstream>
#include <map>
#include "pyeitsolver.h"

void savePotentialsToGmsh(std::vector<double> &sol, const char *filename, const char *meshFileName, int nodesCount, int groundNode, int currentsCount);

int main(int argc, char *argv[])
{
    std::cout << "Hello World!" << std::endl;

    std::map<std::string, int> initDetails = init(argv[1], argv[2]);

    std::cout << initDetails["coeffCount"] << " " << initDetails["currentsCount"] << std::endl;

    std::vector<double> conductivities(initDetails["coeffCount"], 0.380);

    std::vector<double> electrodePotentials = solveFullForwardProblem(conductivities);

    for (size_t i = 0; i < electrodePotentials.size(); i++)
    {
        std::cout << electrodePotentials[i] << " ";
    }

    savePotentialsToGmsh(electrodePotentials, "pyeitsol.msh", argv[1], initDetails["coeffCount"], initDetails["groundNode"], initDetails["currentsCount"]);

    return 0;
}

void savePotentialsToGmsh(std::vector<double> &sol, const char *filename, const char *meshFileName, int nodesCount, int groundNode, int currentsCount)
{
    std::ofstream myfile;
    myfile.open(filename);

    std::ifstream inputfile(meshFileName);
    for (int i = 0; inputfile.eof() != true; i++)
    {
        std::string line;
        std::getline(inputfile, line);
        myfile << line << '\n';
    }

    // Salvando os tensoes nos no's em formato para ser utilizado no gmsh
    int k = 0;
    for (int patterno = 0; patterno < currentsCount; patterno++)
    {
        myfile << "$NodeData\n1\n\"Electric Potential\"\n1\n0.0\n3\n"
               << patterno << "\n1\n"
               << nodesCount << "\n";

        // for (; j < groundNode; j++)
        //     myfile << (j + 1) << "\t" << sol[k++] << "\n";
        // myfile << j + 1 << "\t" << 0 << "\n";
        // for (; j < nodesCount - 1; j++)
        //     myfile << (j + 2) << "\t" << sol[k++] << "\n";
        for (int j = 0; j < nodesCount; j++)
            myfile << (j + 1) << "\t" << sol[k++] << "\n";
        myfile << "$EndNodeData\n";
    }

    myfile.flush();
    myfile.close();
}