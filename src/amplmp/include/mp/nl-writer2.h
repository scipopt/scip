/*
 NL writer2

 NL is a format for representing optimization problems such as linear,
 quadratic, nonlinear, complementarity and constraint programming problems
 in discrete or continuous variables. It is described in the technical report
 "Writing .nl Files" (http://www.cs.sandia.gov/~dmgay/nlwrite.pdf).

 This is a complete reusable C++ implementation of an NL writer.

 Usage:
   /// Write an NL file:
   auto result = mp::WriteNLFile(filenamebase, feeder, utils);
   if (mp::WriteNL_OK != result.first) {
     ...
   }

 where feeder is an object that provides information on model
 components. See NLFeeder for an interface of a feeder class.

 See also NLReader and NLHandler classes.

 Copyright (C) 2023 AMPL Optimization Inc.

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Gleb Belov
 */

#ifndef MP_NLWriter22_H
#define MP_NLWriter22_H

#include <string>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <cstdarg>
#include <cassert>

#include "mp/nl-solver-basics-c.h"
#include "mp/nl-writer2-misc.h"
/// Just for checking compilation
#include "mp/nl-feeder.h"


namespace mp {

/// Typedef WriteNLResult:
/// result code and error message
using WriteNLResult = std::pair<NLW2_WriteNLResultCode, std::string>;

/// Write NL file and any necessary auxiliary files.
/// @param namebase: name without extension,
/// to be augmented by .nl, .col, ...etc.
/// @return Write status.
/// Note that warnings are reported via \a utl.
template <class NLFeeder>
inline WriteNLResult WriteNLFile(
    const std::string& namebase,
    NLFeeder& nlf, NLUtils& utl);


/// NLWriter2.
///
/// mp::NLWriter2 is a zero-overhead NL model writer
/// implemented with inline template code. In particular,
/// it does not store any intermediate model representation.
///
/// Usage: recommended via mp::NLSOL class.
/// If you need the NL writing part only, proceed as follows:
///
/// Write an NL file:
///    WriteNLFile(filenamebase, feeder, utils);
///
/// where feeder is an object that provides information on model
/// components. See NLFeeder for an interface of a feeder class.
///
/// See also NLReader and NLHandler classes.
///
/// @param Params: a specialization of class NLWriter2Params<>.
template <typename Params>
class NLWriter2 :
    public Params::FormatterType {
public:
  using FormatterType = typename Params::FormatterType;
  using FeederType = typename Params::FeederType;

  /// Constructor.
  NLWriter2(FeederType& f, NLHeader h, NLUtils& utl) :
    FormatterType(utl, f.WantNLComments(), f.OutputPrecision()),
    feeder_(f), header_(h) { }

  /// Write to files with the given namebase.
  /// @param namebase will be used for filenames
  /// namebase.nl, namebase.row, etc.
  /// @return Write status.
  WriteNLResult WriteFiles(const std::string& namebase);

  /// Retrieve WriteNLResultCode and error message
  WriteNLResult GetResult() const
  { return result_; }


protected:
  /// Write auxiliary files .row, .slc, .col, etc,
  /// if data provided by the corresp. feeders
  /// (otherwise, any existing files are deleted).
  void WriteAuxFiles(const std::string& namebase);

  /// Write NL file, if any variables
  void WriteNL(const std::string& namebase);


  void WriteRowObjNames(const std::string& namebase);
  void WriteRowDelNames(const std::string& namebase);
  void WriteColNames(const std::string& namebase);
  void WriteUnusedVarNames(const std::string& namebase);
  void WriteFixedVars(const std::string& namebase);
  void WriteObjAdj(const std::string& namebase);


  void WriteNLHeader();
  void WriteFunctions();        // function defs
  void WriteSuffixes();
  void WritePLSOSConstraints();

  void WriteVarBounds();
  void WriteInitialGuesses();
  void WriteConBounds();
  void WriteDualInitialGuesses();

  /// i=0: used in several places,
  /// i>0 increasing: used in 1 con,
  /// i<0 decreasing: used in 1 obj
  void WriteDefinedVariables(int i);

  void WriteConObjExpressions();

  void WriteColumnSizes();

  void WriteLinearConExpr();
  void WriteObjGradients();

  void WriteRandomVariables();


protected:
  ////////////////////////////////////////////////////////
  /// Writer callbacks.
  /// The callbacks are employed by user-specialized feeder
  /// when asked to feed something.

  /// Sparse vector writer.
  /// Is constructed knowing the number
  /// of nozeros to be written.
  template <class Index, class Value>
  class SparseVectorWriter {
  public:
    /// Index type
    using index_type = Index;
    /// Value type
    using value_type = Value;
    /// Construct
    SparseVectorWriter() { }
    /// Not construct(const&)
    SparseVectorWriter(const SparseVectorWriter& ) = delete;
    /// Construct(&&)
    SparseVectorWriter(SparseVectorWriter&& other) {
      SparseVectorWriter svw;
      svw = std::move(other);
    }
    /// Constructor
    SparseVectorWriter(NLWriter2& nlw, size_t n);
    /// Destructor
    ~SparseVectorWriter() { assert(0 == n_entries_); }
    /// No operator=(const&)
    SparseVectorWriter& operator=(
        const SparseVectorWriter& vw) = delete;
    /// operator=(&&)
    SparseVectorWriter& operator=(
        SparseVectorWriter&& other) {
      if (this!=&other) {
        std::swap(p_nlw_, other.p_nlw_);
        std::swap(n_entries_, other.n_entries_);
      }
      return *this;
    }
    /// Write next entry
    void Write(Index , Value );
    /// Number of outstanding elements
    int NLeft() const { return n_entries_; }

  private:
    NLWriter2* p_nlw_ = nullptr;
    size_t n_entries_ = 0;
  };

  /// Typedef SparseDblVecWriter
  using SparseDblVecWriter = SparseVectorWriter<int, double>;
  /// Typedef SparseIntVecWriter
  using SparseIntVecWriter = SparseVectorWriter<int, int>;


  /// Single sparse vector writer factory.
  /// Can be used once to produce a SparseVectorWriter.
  /// At first it writes the number of nonzeros
  /// using the given format string.
  template <class Index, class Value>
  class SingleSparseVecWrtFactory {
  public:
    /// Index type
    using index_type = Index;
    /// Value type
    using value_type = Value;
    /// Vector writer type
    using writer_type = SparseVectorWriter<Index, Value>;
    /// Construct.
    /// @param fmt: format string for printf()
    ///  containing %d for the number of sparse elements
    ///  and finishing with '\n',
    ///  e.g., "x%d\t# initial guess\n"
    SingleSparseVecWrtFactory(NLWriter2& nlw, const char* fmt)
      : nlw_(nlw), fmt_(fmt) { }
    /// Construct. More general.
    /// @param hdr_prn: a lambda printing something
    /// given N nonzeros.
    SingleSparseVecWrtFactory(NLWriter2& nlw,
                              std::function<void(int)> hdr_prn)
      : nlw_(nlw), hdr_prn_(hdr_prn) { }
    /// Create a vector writer
    SparseVectorWriter<Index, Value> MakeVectorWriter(size_t nnz) {
      assert(0 == nInst_++);
      if (fmt_)
        nlw_.apr(nlw_.nm, fmt_, (int)nnz);
      else if (hdr_prn_)
        hdr_prn_(nnz);
      else {
        assert(0);
      }
      return {nlw_, nnz};
    }

  private:
    NLWriter2& nlw_;
    const char* fmt_ = nullptr;
    std::function<void(int)> hdr_prn_;
    int nInst_ = 0;
  };

  /// Typedef SingleSparseDblVecWrtFactory
  using SingleSparseDblVecWrtFactory
  = SingleSparseVecWrtFactory<int, double>;
  /// Typedef SingleSparseIntVecWrtFactory
  using SingleSparseIntVecWrtFactory
  = SingleSparseVecWrtFactory<int, int>;


  /** Declare actual expression writer. */
  class ExprWriter;

  /** Expression argument writer.
   *  An object of this class is returned by FuncPut()
   *  and Oput..() of the parent expression's ExprWriter.
   *  This class acts as a proxy for ExprWriter's
   *  of the arguments.
   *
   *  ExprArgWriter(na=1) is passed to write the root of a tree,
   *  such as of the tree for a certain constraint or
   *  a defined variable, or recursively, for
   *  a user expression type, when user calls EPut().
   */
  class ExprArgWriter {
  public:
    /// Construct for \a na arguments.
    ExprArgWriter(NLWriter2& nlw, int na);
    /// Destruct. Check that all arguments are written.
    ~ExprArgWriter();

    /// Don't copy-construct
    ExprArgWriter(const ExprArgWriter& ) = delete;
    /// Don't assign
    void operator=(const ExprArgWriter& ) = delete;

    /// Do move-construct
    ExprArgWriter(ExprArgWriter&& eaw)
        : nlw_(eaw.nlw_), nargs_(eaw.nargs_) { eaw.nargs_=0; }
    /// Do move-assign
    ExprArgWriter& operator=(ExprArgWriter&& eaw)
    { assert(&nlw_==&eaw.nlw_); std::swap(nargs_, eaw.nargs_); return *this; }

    /// Write the next arg as Feeder's native expression.
    /// This recursively calls Feeder::FeedExpr().
    void EPut(typename FeederType::Expr e);

    /** Write the next arg as 'variable reference'.
     *  0 <= index < num_vars is a solver variable;
     *  index >= num_vars is a defined variable. */
    void VPut(int v, const char* descr="");

    /** Write numeric constant expression. */
    void NPut(double x);

    /** Write string constant expression. */
    void StrPut(const char* );

    /** Write the next arg as function call expression. */
    ExprArgWriter FuncPut(
        int index, int nArgs, const char* descr="");

    /// Write the next arg as AMPL opcode for a unary op.
    /// @return 1-arg writer.
    ExprArgWriter OPut1(int opcode, const char* descr="");
    /// Write AMPL opcode for a binary op.
    ExprArgWriter OPut2(int opcode, const char* descr="");
    /// Write AMPL opcode for a 3-arg op.
    ExprArgWriter OPut3(int opcode, const char* descr="");

    /// Write AMPL opcode for an iterated op (min, exists, sum, etc).
    /// For a piecewise-linear expression, \a nArgs should be
    /// 2*(N slopes) and the arguments are:
    /// break points, slopes, argument variable.
    ExprArgWriter OPutN(
        int opcode, int nArgs, const char* descr="");

    /// Shortcut: OPut1( struct Opcode )
    template <class Opcode> ExprArgWriter OPut1(Opcode oc)
    { return OPut1(oc.code, oc.name); }
    /// Shortcut: OPut2( struct Opcode )
    template <class Opcode> ExprArgWriter OPut2(Opcode oc)
    { return OPut2(oc.code, oc.name); }
    /// Shortcut: OPut3( struct Opcode )
    template <class Opcode> ExprArgWriter OPut3(Opcode oc)
    { return OPut3(oc.code, oc.name); }
    /// Shortcut: OPutN( struct Opcode, int nArgs )
    template <class Opcode> ExprArgWriter OPutN(
        Opcode oc, int nArgs)
    { return OPutN(oc.code, nArgs, oc.name); }

  protected:
    /// Add an argument explicitly.
    /// @return the expr writer for the argument.
    /// Note that the ..Put. methods transfer their call
    /// to the result of an AddArg().
    ExprWriter AddArg();

  private:
    NLWriter2& nlw_;
    int nargs_;
  };


protected:
  /** Actual expression writer. */
  class ExprWriter {
  public:
    /// Construct
    ExprWriter(NLWriter2& nlw) : nlw_(nlw) { }

    /// Don't copy-construct
    ExprWriter(const ExprWriter& ) = delete;
    /// Don't copy-assign
    void operator=(const ExprWriter& ) = delete;

    /// Do move-construct
    ExprWriter(ExprWriter&& ) = default;
    /// Do move-assign
    ExprWriter& operator=(ExprWriter&& ) = default;

    /// Write Feeder's native expression.
    /// This recursively calls Feeder::FeedExpr().
    void EPut(typename FeederType::Expr e);

    /** Write expression 'variable reference'.
     *  0 <= index < num_vars is a solver variable;
     *  index >= num_vars is a defined variable. */
    void VPut(int v, const char* descr="");

    /** Write numeric constant expression. */
    void NPut(double x);

    /** Write string constant expression. */
    void StrPut(const char* );

    /** Write the next arg as function call expression. */
    ExprArgWriter FuncPut(
        int index, int nArgs, const char* descr="");

    /// Write the next arg as AMPL opcode for a unary op.
    /// @return 1-arg writer.
    ExprArgWriter OPut1(int opcode, const char* descr="");
    /// Write AMPL opcode for a binary op.
    ExprArgWriter OPut2(int opcode, const char* descr="");
    /// Write AMPL opcode for a 3-arg op.
    ExprArgWriter OPut3(int opcode, const char* descr="");

    /// Write AMPL opcode for an iterated op (min, exists, sum, etc).
    /// For a piecewise-linear expression, \a nArgs should be
    /// 2*(N slopes) and the arguments are:
    /// break points, slopes, argument variable.
    ExprArgWriter OPutN(
        int opcode, int nArgs, const char* descr="");

    /// Shortcut: OPut1( struct Opcode )
    template <class Opcode> ExprArgWriter OPut1(Opcode oc)
    { return OPut1(oc.code, oc.name); }
    /// Shortcut: OPut2( struct Opcode )
    template <class Opcode> ExprArgWriter OPut2(Opcode oc)
    { return OPut2(oc.code, oc.name); }
    /// Shortcut: OPut3( struct Opcode )
    template <class Opcode> ExprArgWriter OPut3(Opcode oc)
    { return OPut3(oc.code, oc.name); }
    /// Shortcut: OPutN( struct Opcode, int nArgs )
    template <class Opcode> ExprArgWriter OPutN(
        Opcode oc, int nArgs)
    { return OPutN(oc.code, nArgs, oc.name); }

  private:
    NLWriter2& nlw_;
  };


  /** Writer of a defined variable. */
  class DefVarWriter {
  public:
    /// Construct
    DefVarWriter(NLWriter2& nlw, int nnzlin);

    /// Write entries c*var[v] to the linear part
    /// of the defining expression.
    /// All nnz entries should be written before
    /// the nonlinear expression.
    /// This method can be called only once.
    SparseDblVecWriter GetLinExprWriter();

    /// Retrieve the nonlinear expression writer.
    /// Should be used after the linear expression.
    ExprArgWriter GetExprWriter();

  private:
    int nnzlin_;
    NLWriter2& nlw_;
  };


  /** Write a sequence of defined variables,
   *  for example, all such used in a certain constraint. */
  class DefVarWriterFactory {
  public:
    /// Construct.
    /// For the meaning of \a k, see NLFeeder.
    DefVarWriterFactory(NLWriter2& nlw, int k)
      : nlw_(nlw), k_(k) { }

    /// Start writing a defined variable.
    ///
    /// @param index: defined variable index, used to
    /// reference it in subsequent expression graphs.
    /// Thus, the index should be >= num_vars,
    ///   < num_vars+num_common_exprs.
    /// Providing the index explicitly, because classical NL
    /// likes special order of defined variables.
    ///
    /// @param nnz: number of nonzeros in the linear part.
    ///
    /// @param descr: meta-information - what is this variable,
    /// for example, "nl(t[2])".
    /// Providing it here because it's not
    /// included in ColNames().
    ///
    /// @return A callback object writing
    /// a single defined variable.
    DefVarWriter StartDefVar(
        int index, int nnz, const char* descr="");

  private:
    NLWriter2& nlw_;
    int k_;            // for which constraint / obj
  };


  /** Write \a num_vars variable bounds
   *  (all except defined variables). */
  class VarBndWriter {
  public:
    /// Construct
    VarBndWriter(NLWriter2& nlw) : nlw_(nlw) { }
    /// Write range for the next variable.
    void WriteLbUb(double lb, double ub);
    /// Get N written
    int GetNWritten() const { return nWrt_; }

  private:
    NLWriter2& nlw_;
    int nWrt_ = 0;
  };

  /** Write \a num_algebraic_cons constraint bounds. */
  class ConBndWriter {
  public:
    /// Construct
    ConBndWriter(NLWriter2& nlw) : nlw_(nlw) { }
    /// Write range/complementarity for the next constraint.
    template <class AlgConRange>
    void WriteAlgConRange(AlgConRange );
    /// Get N written
    int GetNWritten() const { return nWrt_; }

  private:
    NLWriter2& nlw_;
    int nWrt_ = 0;
  };

  /** Write num_vars+num_rand_vars-1
   *  column sizes */
  class ColSizeWriter {
  public:
    /// Construct
    ColSizeWriter(NLWriter2& nlw, int k)
      : nlw_(nlw), kind_(k) { }
    /// Write next col's size
    void Write(int s) {
      switch(kind_) {
      case 1:
        sum_ += s;
        nlw_.apr(nlw_.nm, "%z\n", sum_);
        break;
      case 2:
        nlw_.apr(nlw_.nm, "%d\n", s);
        break;
      default:
        assert(0 && "why am I here?");
      }
      ++nWrt_;
    }
    /// Report N written
    int GetNWritten() const { return nWrt_; }

  private:
    NLWriter2& nlw_;
    const int kind_;
    std::size_t sum_ = 0;
    int nWrt_ = 0;
  };

  /** Write a sequence of random variables. */
  class RandVarWriterFactory {
  public:
    /// Construct
    RandVarWriterFactory(NLWriter2& nlw) : nlw_(nlw) { }

    /// Start writing a random variable.
    /// (To be finished by feeding the defining expression.)
    ///
    /// @return A callback writing the defining expression.
    ExprArgWriter StartRandVar(
        int index, const char* descr);

  private:
    NLWriter2& nlw_;
  };

  using SuffixIntWriter = SparseVectorWriter<int, int>;
  using SuffixDblWriter = SparseVectorWriter<int, double>;

  /** Suffixes. */
  class SuffixWriterFactory {
  public:
    /// Construct
    SuffixWriterFactory(NLWriter2& nlw)
      : nlw_(nlw) { }
    /// Start writing an int-valued suffix.
    SuffixIntWriter StartIntSuffix(
        const char* name, int kind, int nnz);
    /// Start writing a dbl-valued suffix.
    SuffixDblWriter StartDblSuffix(
        const char* name, int kind, int nnz);

  private:
    NLWriter2& nlw_;
  };

  /** PL-SOS constraints. */
  class PLSOSWriter {
  public:
    /// Construct.
    PLSOSWriter(NLWriter2& nlw) : swf_(nlw) { }
    /// .sos for variables
    SuffixIntWriter StartSOSVars(int nnz);
    /// .sos for constraints
    SuffixIntWriter StartSOSCons(int nnz);
    /// .sosref for variables
    SuffixDblWriter StartSOSREFVars(int nnz);

  private:
    SuffixWriterFactory swf_;
  };

protected:
  /// Write string vector calling the given feeder.
  /// @return max string length.
  int WriteStringVec2File(
      const std::string& name, std::function<void(StringFileWriter&)> );

  /// Writes nothing if !sz
  template <class SparseVecFeeder>
  void WriteSuffix(File& nm,
                   int kind, const char* name,
                   SparseVecFeeder&& sv);

  void WriteSparseEntry(File& nm, int i, int v);
  void WriteSparseEntry(File& nm, int i, double x);

  void WriteBndRangeOrCompl(File& nm,
                            double L, double U, int k=0, int cvar=0);

  template <class ExprArgsFeeder>
  void WriteExprArgs(File& nm, ExprArgsFeeder&& args);

  /// Reuse Writer's facilities
  using FormatterType::apr;
  using FormatterType::nput;
  using FormatterType::Utils;

  /// Retrieve Writer
  FormatterType& Formatter() { return *this; }
  /// Retrieve feeder
  FeederType& Feeder() const { return feeder_; }
  /// Retrieve header
  const NLHeader& Hdr() const { return header_; }

  /// Infinity
  double Infty() const;
  /// Negative infinity
  double NegInfty() const;

private:
  FeederType &feeder_;
  NLHeader header_;

  int num_vars_and_exprs_;  // Number of variables and common expressions.

  int maxLen_ColName_ {0};
  int maxLen_UnvName_ {0};
  int maxLen_FixName_ {0};

  File nm;

  WriteNLResult result_ {NLW2_WriteNL_Unset, ""};
};

}  // namespace mp

#endif // MP_NLWriter22_H
