/*
 Basic definitions for mp::NLWriter2, mp::SOLReader2, mp::NLSolver.
 As extern "C" to be used also in the C API.

 Copyright (C) 2024 AMPL Optimization, Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc. disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */
#ifndef NLW2_BASIC_DEFS_C_H
#define NLW2_BASIC_DEFS_C_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/// Enum NLW2_ObjSense
typedef enum NLW2_ObjSense {
  NLW2_ObjSenseMinimize = 0,
  NLW2_ObjSenseMaximize = 1
} NLW2_ObjSense;

/// Enum NLW2_VarType
typedef enum NLW2_VarType {
  NLW2_VarTypeContinuous = 0,
  NLW2_VarTypeInteger = 1
//const int NLW2_VarTypeSemiContinuous = 2;
//const int NLW2_VarTypeSemiInteger = 3;
//const int NLW2_VarTypeImplicitInteger = 4;
} NLW2_VarType;

/// Enum NLW2_MatrixFormat
typedef enum NLW2_MatrixFormat {
  NLW2_MatrixFormatIrrelevant = -1,
  NLW2_MatrixFormatUnset = 0,
//const int NLW2_MatrixFormatColwise = 1;
  NLW2_MatrixFormatRowwise = 2
} NLW2_MatrixFormat;

/// Enum NLW2_HessianFormat
typedef enum NLW2_HessianFormat {
  NLW2_HessianFormatTriangular = 1,
  NLW2_HessianFormatSquare = 2
} NLW2_HessianFormat;

/// Variables' data by pointers
typedef struct NLW2_ColData_C {
  /// Num vars
  int num_col_;
  /// lower bounds
  const double *lower_;
  /// upper bounds
  const double *upper_;
  /// type: see enum NLW2_VarType.
  /// Set to NULL if all continuous.
  const int *type_;
} NLW2_ColData_C;

/// Sparse matrix.
typedef struct NLW2_SparseMatrix_C {
  /// Size of the start_ array:
  /// N cols (for colwise) / N rows (for rowwise),
  /// depending on format_.
  int num_colrow_;
  /// Format.
  /// Only rowwise supported.
  NLW2_MatrixFormat format_;
  /// Nonzeros
  size_t num_nz_;
  /// Row / col starts
  const size_t *start_;
  /// Entry index
  const int *index_;
  /// Entry value
  const double *value_;
} NLW2_SparseMatrix_C;

/// Sparse vector.
typedef struct NLW2_SparseVector_C {
  /// Number of entries
  int num_;
  /// Entry index
  const int* index_;
  /// Entry value
  const double* value_;
} NLW2_SparseVector_C;


/// Basic NL options for NLModel.
/// Prefer to create by NLW2_MakeNLOptionsBasic_Default().
typedef struct NLW2_NLOptionsBasic_C {
  /// NL text mode?
  int n_text_mode_;
  /// NL comments in text mode?
  int want_nl_comments_;
  /// Flags (1== want output suffixes)
  int flags_;
} NLW2_NLOptionsBasic_C;

/// Use this to create default NL options for NLModel.
NLW2_NLOptionsBasic_C NLW2_MakeNLOptionsBasic_C_Default(void);

/// NLW2_WriteNLResultCode enum.
enum NLW2_WriteNLResultCode {
  NLW2_WriteNL_Unset = 0,
  NLW2_WriteNL_OK = 1,
  NLW2_WriteNL_CantOpen,
  NLW2_WriteNL_Failed
};

/// SOL read result code
enum NLW2_SOLReadResultCode {
  NLW2_SOLRead_Result_Not_Set = -1,
  NLW2_SOLRead_OK = 0,
  NLW2_SOLRead_Fail_Open,
  NLW2_SOLRead_Early_EOF,
  NLW2_SOLRead_Bad_Format,
  NLW2_SOLRead_Bad_Line,
  NLW2_SOLRead_Bad_Options,
  NLW2_SOLRead_Vector_Not_Finished,
  NLW2_SOLRead_Bad_Suffix
};


#ifdef __cplusplus
}  // extern "C"
#endif

#endif // NLW2_BASIC_DEFS_C_H
