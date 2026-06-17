#pragma once

#include "Utils/GMM_Macros.h"


namespace MatlabGMMDataExchange
{
	int SetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A);
	int SetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A);
	int SetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A);
	int SetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A);

	int GetEngineDenseMatrix(const char* name, GMMDenseComplexColMatrix& A);
	int GetEngineDenseMatrix(const char* name, GMMDenseColMatrix& A);
	int GetEngineSparseMatrix(const char* name, GMMSparseRowMatrix& A);
	int GetEngineSparseMatrix(const char* name, GMMSparseComplexRowMatrix& A);
	int GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0ComplexRowMatrix& A);
	int GetEngineCompressedSparseMatrix(const char* name, GMMCompressed0RowMatrix& A);
}
