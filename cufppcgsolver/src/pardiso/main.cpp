/*******************************************************************************
 * Copyright 2004-2022 Intel Corporation.
 *
 * This software and the related documents are Intel copyrighted  materials,  and
 * your use of  them is  governed by the  express license  under which  they were
 * provided to you (License).  Unless the License provides otherwise, you may not
 * use, modify, copy, publish, distribute,  disclose or transmit this software or
 * the related documents without Intel's prior written permission.
 *
 * This software and the related documents  are provided as  is,  with no express
 * or implied  warranties,  other  than those  that are  expressly stated  in the
 * License.
 *******************************************************************************/

/*
 *   Content : Intel(R) oneAPI Math Kernel Library (oneMKL) PARDISO C example
 *
 ********************************************************************************
 */
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "PARDISO" routine */
/* on symmetric linear systems */
/* -------------------------------------------------------------------- */
/* This program can be downloaded from the following site: */
/* www.pardiso-project.org */
/* */
/* (C) Olaf Schenk, Department of Computer Science, */
/* University of Basel, Switzerland. */
/* Email: olaf.schenk@unibas.ch */
/* -------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "../args/args.hxx"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mmio.h"

template <typename T_ELEM>
MKL_INT loadMMSparseMatrix(char *filename, char elem_type, bool csrFormat, MKL_INT *m,
                           MKL_INT *n, MKL_INT *nnz, T_ELEM **aVal, MKL_INT **aRowInd,
                           MKL_INT **aColInd, bool extendSymMatrix, bool transposeMatrix);

// Define the format to printf MKL_INT values
#if !defined(MKL_ILP64)
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

int main(int argc, char *argv[])
{
    /* Parse arguments */
    args::ArgumentParser parser("This is a performance test program for the pardiso solver.", "No comments.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::CompletionFlag completion(parser, {"complete"});
    args::ValueFlag<std::string> bfname(parser, "filename", "b vector file", {'b'});
    args::Positional<std::string> Afname(parser, "filename", "A matrix file. Supported files: .mtx or .msh (with automatic conversion)");
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Completion e)
    {
        std::cout << e.what();
        return 0;
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl
                  << parser;
        return 1;
    }

    /* Matrix data. */
    MKL_INT m, n, nnz;
    MKL_INT *ia, *ja;
    double *a;
    MKL_INT mtype = -2; /* Real symmetric matrix */
    std::string fileName = args::get(Afname);

    if (loadMMSparseMatrix<double>((char *)fileName.c_str(), 'd', true, &m,
                                   &n, &nnz, &a, &ia,
                                   &ja, false, true))
    {
        exit(EXIT_FAILURE);
    }

    /* RHS and solution vectors. */
    double *b, *x;
    b = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));
    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i;
    double ddum;  /* Double dummy */
    MKL_INT idum; /* Integer dummy. */
                  /* -------------------------------------------------------------------- */
                  /* .. Setup Pardiso control parameters. */
                  /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;   /* No solver default */
    iparm[1] = 2;   /* Fill-in reordering from METIS */
    iparm[3] = 0;   /* No iterative-direct algorithm */
    iparm[4] = 0;   /* No user fill-in reducing permutation */
    iparm[5] = 0;   /* Write solution into x */
    iparm[6] = 0;   /* Not in use */
    iparm[7] = 2;   /* Max numbers of iterative refinement steps */
    iparm[8] = 0;   /* Not in use */
    iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Not in use */
    iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    maxfct = 1;     /* Maximum number of numerical factorizations. */
    mnum = 1;       /* Which factorization to use. */
    msglvl = 1;     /* Print statistical information in file */
    error = 0;      /* Initialize error flag */
                    /* -------------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: " IFORMAT, error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors = " IFORMAT, iparm[17]);
    printf("\nNumber of factorization MFLOPS = " IFORMAT, iparm[18]);
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: " IFORMAT, error);
        exit(2);
    }
    printf("\nFactorization completed ... ");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2; /* Max numbers of iterative refinement steps. */
    /* Set right hand side to one. */
    for (i = 0; i < n; i++)
    {
        b[i] = 1;
    }
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if (error != 0)
    {
        printf("\nERROR during solution: " IFORMAT, error);
        exit(3);
    }
    printf("\nSolve completed ... ");
    printf("\nThe solution of the system is: ");
    for (i = 0; i < n; i++)
    {
        printf("\n x [" IFORMAT "] = % f", i, x[i]);
    }
    printf("\n");
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    return 0;
}

static void compress_index(
    const MKL_INT *Ind,
    MKL_INT nnz,
    MKL_INT m,
    MKL_INT *Ptr,
    MKL_INT base)
{
    int i;

    /* initialize everything to zero */
    for (i = 0; i < m + 1; i++)
    {
        Ptr[i] = 0;
    }
    /* count elements in every row */
    Ptr[0] = base;
    for (i = 0; i < nnz; i++)
    {
        Ptr[Ind[i] + (1 - base)]++;
    }
    /* add all the values */
    for (i = 0; i < m; i++)
    {
        Ptr[i + 1] += Ptr[i];
    }
}

struct cooFormat
{
    MKL_INT i;
    MKL_INT j;
    MKL_INT p; // permutation
};

int cmp_cooFormat_csr(struct cooFormat *s, struct cooFormat *t)
{
    if (s->i < t->i)
    {
        return -1;
    }
    else if (s->i > t->i)
    {
        return 1;
    }
    else
    {
        return (int)(s->j - t->j);
    }
}

int cmp_cooFormat_csc(struct cooFormat *s, struct cooFormat *t)
{
    if (s->j < t->j)
    {
        return -1;
    }
    else if (s->j > t->j)
    {
        return 1;
    }
    else
    {
        return (int)(s->i - t->i);
    }
}

typedef int (*FUNPTR)(const void *, const void *);
typedef int (*FUNPTR2)(struct cooFormat *s, struct cooFormat *t);

static FUNPTR2 fptr_array[2] = {
    cmp_cooFormat_csr,
    cmp_cooFormat_csc,
};

static int verify_pattern(
    MKL_INT m,
    MKL_INT nnz,
    MKL_INT *csrRowPtr,
    MKL_INT *csrColInd)
{
    MKL_INT i, col, start, end, base_index;
    int error_found = 0;

    if (nnz != (csrRowPtr[m] - csrRowPtr[0]))
    {
        fprintf(stderr, "Error (nnz check failed): (csrRowPtr[%d]=%lld - csrRowPtr[%lld]=%lld) != (nnz=%lld)\n", 0, csrRowPtr[0], m, csrRowPtr[m], nnz);
        error_found = 1;
    }

    base_index = csrRowPtr[0];
    if ((0 != base_index) && (1 != base_index))
    {
        fprintf(stderr, "Error (base index check failed): base index = %lld\n", base_index);
        error_found = 1;
    }

    for (i = 0; (!error_found) && (i < m); i++)
    {
        start = csrRowPtr[i] - base_index;
        end = csrRowPtr[i + 1] - base_index;
        if (start > end)
        {
            fprintf(stderr, "Error (corrupted row): csrRowPtr[%lld] (=%lld) > csrRowPtr[%lld] (=%lld)\n", i, start + base_index, i + 1, end + base_index);
            error_found = 1;
        }
        for (col = start; col < end; col++)
        {
            if (csrColInd[col] < base_index)
            {
                fprintf(stderr, "Error (column vs. base index check failed): csrColInd[%lld] < %lld\n", col, base_index);
                error_found = 1;
            }
            if ((col < (end - 1)) && (csrColInd[col] >= csrColInd[col + 1]))
            {
                fprintf(stderr, "Error (sorting of the column indecis check failed): (csrColInd[%lld]=%lld) >= (csrColInd[%lld]=%lld)\n", col, csrColInd[col], col + 1, csrColInd[col + 1]);
                error_found = 1;
            }
        }
    }
    return error_found;
}

template <typename T_ELEM>
MKL_INT loadMMSparseMatrix(
    char *filename,
    char elem_type,
    bool csrFormat,
    MKL_INT *m,
    MKL_INT *n,
    MKL_INT *nnz,
    T_ELEM **aVal,
    MKL_INT **aRowInd,
    MKL_INT **aColInd,
    bool extendSymMatrix, bool transposeMatrix)
{
    MM_typecode matcode;
    double *tempVal;
    MKL_INT *tempRowInd, *tempColInd;
    double *tval;
    int *trow, *tcol;
    MKL_INT *csrRowPtr, *cscColPtr;
    MKL_INT i, j, error, base, count;
    struct cooFormat *work;

    /* read the matrix */
    int i_m, i_n, i_nnz;
    error = mm_read_mtx_crd(filename, &i_m, &i_n, &i_nnz, &trow, &tcol, &tval, &matcode);
    if (error)
    {
        fprintf(stderr, "!!!! can not open file: '%s'\n", filename);
        return 1;
    }
    *m = (MKL_INT)i_m;
    *n = (MKL_INT)i_n;
    *nnz = (MKL_INT)i_nnz;

    /* start error checking */
    if (mm_is_complex(matcode) && ((elem_type != 'z') && (elem_type != 'c')))
    {
        fprintf(stderr, "!!!! complex matrix requires type 'z' or 'c'\n");
        return 1;
    }

    if (mm_is_dense(matcode) || mm_is_array(matcode) || mm_is_pattern(matcode) /*|| mm_is_integer(matcode)*/)
    {
        fprintf(stderr, "!!!! dense, array, pattern and integer matrices are not supported\n");
        return 1;
    }

    /* if necessary symmetrize the pattern (transform from triangular to full) */
    if ((extendSymMatrix) && (mm_is_symmetric(matcode) || mm_is_hermitian(matcode) || mm_is_skew(matcode)))
    {
        // count number of non-diagonal elements
        count = 0;
        for (i = 0; i < (*nnz); i++)
        {
            if (trow[i] != tcol[i])
            {
                count++;
            }
        }
        // allocate space for the symmetrized matrix
        tempRowInd = (MKL_INT *)malloc((*nnz + count) * sizeof(MKL_INT));
        tempColInd = (MKL_INT *)malloc((*nnz + count) * sizeof(MKL_INT));
        if (mm_is_real(matcode) || mm_is_integer(matcode))
        {
            tempVal = (double *)malloc((*nnz + count) * sizeof(double));
        }
        else
        {
            tempVal = (double *)malloc(2 * (*nnz + count) * sizeof(double));
        }
        // copy the elements regular and transposed locations
        for (j = 0, i = 0; i < (*nnz); i++)
        {
            tempRowInd[j] = (MKL_INT)trow[i];
            tempColInd[j] = (MKL_INT)tcol[i];
            if (mm_is_real(matcode) || mm_is_integer(matcode))
            {
                tempVal[j] = tval[i];
            }
            else
            {
                tempVal[2 * j] = tval[2 * i];
                tempVal[2 * j + 1] = tval[2 * i + 1];
            }
            j++;
            if (trow[i] != tcol[i])
            {
                tempRowInd[j] = (MKL_INT)tcol[i];
                tempColInd[j] = (MKL_INT)trow[i];
                if (mm_is_real(matcode) || mm_is_integer(matcode))
                {
                    if (mm_is_skew(matcode))
                    {
                        tempVal[j] = -tval[i];
                    }
                    else
                    {
                        tempVal[j] = tval[i];
                    }
                }
                else
                {
                    if (mm_is_hermitian(matcode))
                    {
                        tempVal[2 * j] = tval[2 * i];
                        tempVal[2 * j + 1] = -tval[2 * i + 1];
                    }
                    else
                    {
                        tempVal[2 * j] = tval[2 * i];
                        tempVal[2 * j + 1] = tval[2 * i + 1];
                    }
                }
                j++;
            }
        }
        (*nnz) += count;
        // free temporary storage
        free(trow);
        free(tcol);
        free(tval);
    }
    else
    {
        tempRowInd = (MKL_INT *)malloc(*nnz * sizeof(MKL_INT));
        tempColInd = (MKL_INT *)malloc(*nnz * sizeof(MKL_INT));
        for (i = 0; i < *nnz; i++)
        {
            tempRowInd[i] = (MKL_INT)(transposeMatrix ? tcol[i] : trow[i]);
            tempColInd[i] = (MKL_INT)(transposeMatrix ? trow[i] : tcol[i]);
        }
        tempVal = tval;
    }
    // life time of (trow, tcol, tval) is over.
    // please use COO format (tempRowInd, tempColInd, tempVal)

    // use qsort to sort COO format
    work = (struct cooFormat *)malloc(sizeof(struct cooFormat) * (*nnz));
    if (NULL == work)
    {
        fprintf(stderr, "!!!! allocation error, malloc failed\n");
        return 1;
    }
    for (i = 0; i < (*nnz); i++)
    {
        work[i].i = tempRowInd[i];
        work[i].j = tempColInd[i];
        work[i].p = i; // permutation is identity
    }

    if (csrFormat)
    {
        /* create row-major ordering of indices (sorted by row and within each row by column) */
        qsort(work, *nnz, sizeof(struct cooFormat), (FUNPTR)fptr_array[0]);
    }
    else
    {
        /* create column-major ordering of indices (sorted by column and within each column by row) */
        qsort(work, *nnz, sizeof(struct cooFormat), (FUNPTR)fptr_array[1]);
    }

    // (tempRowInd, tempColInd) is sorted either by row-major or by col-major
    for (i = 0; i < (*nnz); i++)
    {
        tempRowInd[i] = work[i].i;
        tempColInd[i] = work[i].j;
    }

    // setup base
    // check if there is any row/col 0, if so base-0
    // check if there is any row/col equal to matrix dimension m/n, if so base-1
    MKL_INT base0 = 0;
    MKL_INT base1 = 0;
    for (i = 0; i < (*nnz); i++)
    {
        const MKL_INT row = tempRowInd[i];
        const MKL_INT col = tempColInd[i];
        if ((0 == row) || (0 == col))
        {
            base0 = 1;
        }
        if ((*m == row) || (*n == col))
        {
            base1 = 1;
        }
    }
    if (base0 && base1)
    {
        printf("Error: input matrix is base-0 and base-1 \n");
        return 1;
    }

    base = 0;
    if (base1)
    {
        base = 1;
    }

    /* compress the appropriate indices */
    if (csrFormat)
    {
        /* CSR format (assuming row-major format) */
        csrRowPtr = (MKL_INT *)malloc(((*m) + 1) * sizeof(csrRowPtr[0]));
        if (!csrRowPtr)
            return 1;
        compress_index(tempRowInd, *nnz, *m, csrRowPtr, base);

        *aRowInd = csrRowPtr;
        *aColInd = (MKL_INT *)malloc((*nnz) * sizeof(MKL_INT));
    }
    else
    {
        /* CSC format (assuming column-major format) */
        cscColPtr = (MKL_INT *)malloc(((*n) + 1) * sizeof(cscColPtr[0]));
        if (!cscColPtr)
            return 1;
        compress_index(tempColInd, *nnz, *n, cscColPtr, base);

        *aColInd = cscColPtr;
        *aRowInd = (MKL_INT *)malloc((*nnz) * sizeof(MKL_INT));
    }

    /* transfrom the matrix values of type double into one of the cusparse library types */
    *aVal = (T_ELEM *)malloc((*nnz) * sizeof(T_ELEM));

    for (i = 0; i < (*nnz); i++)
    {
        if (csrFormat)
        {
            (*aColInd)[i] = tempColInd[i];
        }
        else
        {
            (*aRowInd)[i] = tempRowInd[i];
        }
        (*aVal)[i] = tempVal[work[i].p];
    }

    /* check for corruption */
    MKL_INT error_found;
    if (csrFormat)
    {
        error_found = verify_pattern(*m, *nnz, *aRowInd, *aColInd);
    }
    else
    {
        error_found = verify_pattern(*n, *nnz, *aColInd, *aRowInd);
    }
    if (error_found)
    {
        fprintf(stderr, "!!!! verify_pattern failed\n");
        return 1;
    }

    /* cleanup and exit */
    free(work);
    free(tempVal);
    free(tempColInd);
    free(tempRowInd);

    return 0;
}

/* specific instantiation */
template MKL_INT loadMMSparseMatrix<double>(
    char *filename,
    char elem_type,
    bool csrFormat,
    MKL_INT *m,
    MKL_INT *n,
    MKL_INT *nnz,
    double **aVal,
    MKL_INT **aRowInd,
    MKL_INT **aColInd,
    bool extendSymMatrix, bool transposeMatrix);