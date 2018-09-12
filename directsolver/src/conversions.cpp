#include "conversions.h"
#include <string>
#include "problem.h"
#include <memory>
#include <experimental/filesystem>
extern "C" {
	#include "mm/mmio.h"
}

namespace EITFILECONVERISONS {
	/*
	* Get File name without extension from File path or File Name
	*/
	std::string getFilenameOnly(std::string filePath)
	{
		// Find the last position of '.' in given string
		std::size_t pos = filePath.rfind('.');

		// If last '.' is found
		if (pos != std::string::npos) {
			// return the substring
			return filePath.substr(0,pos);
		}
		// In case of no extension return empty string
		return "";
	}


	void saveMtx(matrix *m, FILE *f) {
		MM_typecode matcode;
		int i;

		mm_initialize_typecode(&matcode);
		mm_set_matrix(&matcode);
		mm_set_coordinate(&matcode);
		mm_set_real(&matcode);
		mm_set_symmetric(&matcode);

		mm_write_banner(f, matcode);
		mm_write_mtx_crd_size(f, m->rows(), m->cols(), m->nonZeros());

		/* NOTE: matrix market files use 1-based indices, i.e. first element
		of a vector has index 1, not 0.  */

		for (int k = 0; k < m->outerSize(); ++k)
			for (matrix::InnerIterator it(*m, k); it; ++it)
				if(it.row() >= it.col())
					fprintf(f, "%d %d %10.16g\n", it.row() + 1, it.col() + 1, it.value());
	}

	std::string convertMeshFile(const std::string filename) {
		std::experimental::filesystem::path filenamepath(filename);
		if (!std::experimental::filesystem::exists(filenamepath)) return "";

		bool is2dProblem;
		std::shared_ptr<problem> input = problem::createNewProblem(filename.c_str(), is2dProblem);
		input->initProblem(filename.c_str());
		input->buildNodeCoefficients();
		input->prepareSkeletonMatrix();
		input->createCoef2KMatrix();

		matrix *m;
		Eigen::VectorXd vcond(input->getNumCoefficients());
		for (int i = 0; i < vcond.rows(); i++) vcond[i] = 0.3815;
		input->assembleProblemMatrix(&vcond[0], &m);
		input->postAssembleProblemMatrix(&m);


		FILE *f;
		std::string mtxFname = getFilenameOnly(filename) + ".mtx";
		if ((f = fopen(mtxFname.c_str(), "w")) == NULL) { std::cerr << "Could not open market matrix file to write" << mtxFname << std::endl; return ""; }
		saveMtx(m,f);
		fclose(f);

		return mtxFname;
	}
}