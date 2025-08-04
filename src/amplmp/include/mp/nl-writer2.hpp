/*
 Copyright (C) 2024 AMPL Optimization Inc.

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
#ifndef NLWriter22_HPP
#define NLWriter22_HPP

#include <limits>
#include <cstring>

#include "mp/nl-writer2.h"

namespace mp {

template <class NLFeeder>
inline WriteNLResult WriteNLFile(
    const std::string& namebase,
    NLFeeder& nlf, NLUtils& utl) {
  auto h = nlf.Header();
  if (NLHeader::TEXT == h.format) {
    mp::NLWriter2<
        mp::NLWriter2Params<
        mp::TextFormatter, NLFeeder > >
        writer(nlf, h, utl);
    return writer.WriteFiles(namebase);
  } else {
    mp::NLWriter2<
        mp::NLWriter2Params<
        mp::BinaryFormatter, NLFeeder > >
        writer(nlf, h, utl);
    return writer.WriteFiles(namebase);
  }
}

template <typename Params>
WriteNLResult NLWriter2<Params>::WriteFiles(
    const std::string& namebase) {
  WriteAuxFiles(namebase);
  WriteNL(namebase);
  return GetResult();
}

template <typename Params>
void NLWriter2<Params>::WriteAuxFiles(
    const std::string& namebase) {
  WriteRowObjNames(namebase);
  WriteRowDelNames(namebase);
  WriteColNames(namebase);
  WriteUnusedVarNames(namebase);
  WriteFixedVars(namebase);
  header_.max_var_name_len =     // maxLen_FixName_ can be longer
      maxLen_ColName_ + maxLen_UnvName_ + maxLen_FixName_;
  WriteObjAdj(namebase);
}

template <typename Params>
void NLWriter2<Params>::WriteNL(
    const std::string& namebase) {
  auto fln = namebase + ".nl";
  nm = Utils().openf(fln,
                     Hdr().num_vars<=0,     // Delete file if true
                     "wb");
  if (nm) {
    WriteNLHeader();
    WriteFunctions();
    WriteSuffixes();
    WritePLSOSConstraints();
    if (Feeder().WantBoundsFirst()) {
      WriteVarBounds();
      WriteInitialGuesses();
      WriteConBounds();
      WriteDualInitialGuesses();
    }
    WriteDefinedVariables(0);
    WriteConObjExpressions();
    if (!Feeder().WantBoundsFirst()) {
      WriteDualInitialGuesses();
      WriteInitialGuesses();
      WriteConBounds();
      WriteVarBounds();
    }
    if (Feeder().WantColumnSizes())
      WriteColumnSizes();
    WriteLinearConExpr();
    WriteObjGradients();
    WriteRandomVariables();

    if (std::ferror(nm.GetHandle())) {
      result_.first = NLW2_WriteNL_Failed;
      result_.second = fln + ": " + std::strerror(errno);
    } else
      result_.first = NLW2_WriteNL_OK;
  } else {
    result_.first = NLW2_WriteNL_CantOpen;
    result_.second = fln + ": " + std::strerror(errno);
  }
}


/// Write strings to a file.
/// The file is only opened when the user calls operator bool()
/// or upon 1st write.
/// If nothing written, the file is removed (if exists).
class StringFileWriter {
public:
  /// bool: true if should remove the file
  using TOpener = std::function<File(bool)>;
  /// Construct
  StringFileWriter(int& lmax, TOpener opener)
    : len_max_(lmax), opener_(opener) { }
  /// No copy
  StringFileWriter(const StringFileWriter& ) = delete;
  /// Destroy
  ~StringFileWriter() {
    if (!cnt_ && !fTriedOpen_)  // remove file if not opened
      opener_(true);            // or nothing written
  }
  /// No assign
  void operator=(const StringFileWriter& ) = delete;
  /// oper bool()
  operator bool() {
    return !fTriedOpen_ ? Open() : nm;
  }
  /// Write const char*, StrStrValue, StrDblValue
  template <class Value>
  void Write(Value s) {
    if (!fTriedOpen_)
      Open();
    ++cnt_;
    auto l = Print_ReturnLen1(s);
    if (l>len_max_)
      len_max_ = l;
  }
  /// Typedef StrStrValue
  using StrStrValue = std::pair<const char*, const char*>;
  /// Typedef StrDblValue
  using StrDblValue = std::pair<const char*, double>;

protected:
  bool Open() {
    fTriedOpen_ = true;
    nm = opener_(false);
    return nm;
  }
  int Print_ReturnLen1(const char* s) {
    nm.Printf("%s\n", s);
    return std::strlen(s);
  }
  int Print_ReturnLen1(StrStrValue p) {
    nm.Printf("%s\t%s\n", p.first, p.second);
    return std::strlen(p.first);
  }
  int Print_ReturnLen1(StrDblValue p) {
    nm.Printf("%s\t%.17g\n", p.first, p.second);
    return std::strlen(p.first);
  }

private:
  int& len_max_;
  TOpener opener_;
  File nm;
  bool fTriedOpen_{false};
  int cnt_{0};
};

/// operator<< (StringFileWriter&, Value)
template <class Value>
StringFileWriter&
operator << (StringFileWriter& wrt, Value s) {
  wrt.Write(s);
  return wrt;
}

template <typename Params>
void NLWriter2<Params>::WriteRowObjNames(
    const std::string& namebase) {
  header_.max_con_name_len =
      WriteStringVec2File(namebase + ".row",
                          [this](StringFileWriter& w){
      feeder_.FeedRowAndObjNames(w); });
}

template <typename Params>
void NLWriter2<Params>::WriteRowDelNames(
    const std::string& namebase) {
  WriteStringVec2File(namebase + ".slc",
                      [this](StringFileWriter& w){
    feeder_.FeedDelRowNames(w); });
}

template <typename Params>
void NLWriter2<Params>::WriteColNames(
    const std::string& namebase) {
  maxLen_ColName_ =
      WriteStringVec2File(namebase + ".col",
                          [this](StringFileWriter& w){
      feeder_.FeedColNames(w); });
}

template <typename Params>
void NLWriter2<Params>::WriteUnusedVarNames(
    const std::string& namebase) {
  maxLen_UnvName_ =
      WriteStringVec2File(namebase + ".unv",
                          [this](StringFileWriter& w){
      feeder_.FeedUnusedVarNames(w);
});
}

template <typename Params>
void NLWriter2<Params>::WriteFixedVars(
    const std::string& namebase) {
  maxLen_FixName_ =
      WriteStringVec2File(namebase + ".fix",
                          [this](StringFileWriter& w){
      feeder_.FeedFixedVarNames(w);
});
}

template <typename Params>
void NLWriter2<Params>::WriteObjAdj(
    const std::string& namebase) {
  WriteStringVec2File(namebase + ".adj",
                      [this](StringFileWriter& w){
    feeder_.FeedObjAdj(w);
  });
}

template <typename Params>
int NLWriter2<Params>::WriteStringVec2File(
    const std::string& name, std::function<void(StringFileWriter&)> swf) {
  int lMax=0;
  StringFileWriter sw(lMax,
                      [name, this](bool fClose){
    return Utils().openf(name.c_str(),
                         fClose, "w");
  });
  swf(sw);
  return lMax;
}


/* song and dance because of the interaction of cfront's treatment
   of adjacent string literals, lcc's treatment of multi-line strings,
   and MSDOS's stupid \r\n conventions....
 */
#ifdef X_AMPL
#if defined(_WIN32) || defined(_WIN64)
#define EOL "\r\n"
#else
#define EOL "\n"
#endif /* _WIN32 && _WIN64 */
#else /* X_AMPL */
#ifdef MSDOS
#define EOL "\r\n"
#else
#define EOL "\n"
#endif
#endif /* X_AMPL */

#ifdef NL_LIB2_ORIG_HDR
static constexpr char
gl_1__short[]	= "%c%d",
gl_2a[]	= " %d %d %d %d %d",
gl_1a[]	= "\t# problem %s" EOL,
gl_2[]	= "\t# vars, constraints, objectives, ranges, eqns%s" EOL,
gl_3[]	= " %d %d\t# nonlinear constraints, objectives" EOL,
gl_3c[]	= "\t# nonlinear constrs, objs; ccons: lin, nonlin, nd, nzlb%s" EOL,
gl_4[]	= "\t# network constraints: nonlinear, linear" EOL,
gl_4r[]	= " %d\t# network constraints: nonlinear, linear; stages" EOL,
gl_5[]	= " %d %d %d\t# nonlinear vars in constraints, objectives, both" EOL,
gl_6[]	= " %d %d\t# linear network variables; functions" EOL,
gl_6x[]	= " %d %d %d %d\t# linear network variables; functions; arith, flags" EOL,
gl_6y[]	= " %d %d %d %d %d\t# linear network variables; functions; arith, flags; randcalls" EOL,
gl_7[]	= " %d %d %d %d %d\t# discrete variables: binary, integer, nonlinear (b,c,o)" EOL,
gl_8[]	= " %zd %zd\t# nonzeros in Jacobian, gradients" EOL,
gl_9[]	= " %d %d\t# max name lengths: constraints, variables" EOL,
gl_10[]	= "\t# common exprs: b,c,o,c1,o1%s" EOL;
#else
static constexpr char
gl_1__short[]	= "%c%d",
gl_2a[]	= " %d %d %d %d %d",
gl_1a[]	= "\t# problem %s" EOL,
gl_2[]	= "\t# vars, algcons, objs, ranges, eqns%s" EOL,
gl_3[]	= " %d %d\t# nonlinear cons, objs" EOL,
gl_3c[]	= "\t# nonlinear cons, objs; compl: lin, nonlin, range, nzlb%s" EOL,
gl_4[]	= "\t# network cons: nonlinear, linear" EOL,
gl_4r[]	= " %d\t# network cons: nonlinear, linear; stages" EOL,
gl_5[]	= " %d %d %d\t# nonlinear vars in cons, objs, both" EOL,
gl_6[]	= " %d %d\t# linear network vars; funcs" EOL,
gl_6x[]	= " %d %d %d %d\t# linear network vars; funcs; arith, flags" EOL,
gl_6y[]	= " %d %d %d %d %d\t# linear network vars; funcs; arith, flags; randcalls" EOL,
gl_7[]	= " %d %d %d %d %d\t# discrete vars: binary, integer, nonlinear (b,c,o)" EOL,
gl_8[]	= " %zd %zd\t# nonzeros in Jacobian, gradients" EOL,
gl_9[]	= " %d %d\t# max name lengths: cons/objs, vars" EOL,
gl_10[]	= "\t# common exprs: b,c,o,c1,o1%s" EOL;
#endif
#undef EOL


template <typename Params>
void NLWriter2<Params>::WriteNLHeader() {
  assert(Formatter().Mode() == Hdr().format);

  /// But the header is always text.
  nm.Printf(gl_1__short,
            NLHeader::TEXT==Hdr().format ? 'g' : 'b',
            Hdr().num_ampl_options);
  for (int i = 0; i < Hdr().num_ampl_options; ++i)
    nm.Printf(" %ld", Hdr().ampl_options[i]);
  if (Hdr().ampl_options[VBTOL_OPTION_INDEX] == USE_VBTOL_FLAG)
    nm.Printf(" %.g", Hdr().ampl_vbtol);
  nm.Printf(gl_1a, Hdr().prob_name);

  /// Num variables, constraints, obj, ...
  nm.Printf(gl_2a,
            Hdr().num_vars,
            Hdr().num_algebraic_cons,
            Hdr().num_objs,
            Hdr().num_ranges,
            Hdr().num_eqns);
  const char* s = "";
  if (Hdr().num_rand_vars) {            /*SNL2006*/
    s = ", lcons, randvars";
    nm.Printf(" %d %d",
              Hdr().num_logical_cons,
              Hdr().num_rand_vars);
  }
  else if (Hdr().num_logical_cons) {
    s = ", lcons";
    nm.Printf(" %d", Hdr().num_logical_cons);
  }
  nm.Printf(gl_2, s);

  /// Nonlinear / random cons/obj
  if (Hdr().num_compl_conds |
      Hdr().num_rand_cons | Hdr().num_rand_objs) {
    nm.Printf(" %d %d %d %d %d %d",
              Hdr().num_nl_cons, Hdr().num_nl_objs,
              Hdr().num_compl_conds - Hdr().num_nl_compl_conds,
              Hdr().num_nl_compl_conds, Hdr().num_compl_dbl_ineqs,
              Hdr().num_compl_vars_with_nz_lb);
    s = "";
    if (Hdr().num_rand_cons | Hdr().num_rand_objs) {
      nm.Printf(" %d %d",
                Hdr().num_rand_cons, Hdr().num_rand_objs);
      s = "; rand constrs, objs";
    }
    nm.Printf(gl_3c, s);
  }
  else
    nm.Printf(gl_3, Hdr().num_nl_cons, Hdr().num_nl_objs);

  /// Network cons, stages
  nm.Printf(" %d %d",
            Hdr().num_nl_net_cons, Hdr().num_linear_net_cons);
  nm.Printf(Hdr().num_stages > 1 ? gl_4r : gl_4,
            Hdr().num_stages);

  /// NL vars in cons, objs, both
  nm.Printf(gl_5, Hdr().num_nl_vars_in_cons,
            Hdr().num_nl_vars_in_objs, Hdr().num_nl_vars_in_both);

  /// Linear network vars, functions, ak, flags, rnd calls
  const char* fmt =
      Hdr().num_rand_vars ? gl_6y :
                            Hdr().flags | Hdr().arith_kind ?
                              gl_6x : gl_6;
  nm.Printf(fmt, Hdr().num_linear_net_vars, Hdr().num_funcs,
            NLHeader::TEXT==Hdr().format ? 0 : Hdr().arith_kind,
            Hdr().flags, Hdr().num_rand_calls);

  /// Integer vars
  nm.Printf(gl_7, Hdr().num_linear_binary_vars,
            Hdr().num_linear_integer_vars,
            Hdr().num_nl_integer_vars_in_both,
            Hdr().num_nl_integer_vars_in_cons,
            Hdr().num_nl_integer_vars_in_objs);

  /// Nonzeros
  nm.Printf(gl_8, Hdr().num_con_nonzeros,
            Hdr().num_obj_nonzeros);

  /// Max name lengths
  nm.Printf(gl_9, Hdr().max_con_name_len,
            Hdr().max_var_name_len);

  /// Common exprs
  nm.Printf(" %d %d %d %d %d",
            Hdr().num_common_exprs_in_both,
            Hdr().num_common_exprs_in_cons,
            Hdr().num_common_exprs_in_objs,
            Hdr().num_common_exprs_in_single_cons,
            Hdr().num_common_exprs_in_single_objs);
  s = "";
  if (Hdr().num_rand_common_exprs) {       /*SNL2006*/
    nm.Printf(" %d", Hdr().num_rand_common_exprs);
    s = ",rand";
  }
  nm.Printf(gl_10, s);
}

template <typename Params>
void NLWriter2<Params>::WriteFunctions() {
  for (int i=0; i<Hdr().num_funcs; ++i) {
    auto func = Feeder().Function(i);
    apr(nm, "F%d %d %d %s\n", i,
        (int)func.Type(),
        (int)func.NumArgs(),
        func.Name());
  }
}

template <typename Params>
void NLWriter2<Params>::
WriteSparseEntry(File& nm_, int i, int v) {
  apr(nm_, "%d %d\n", i, v);
}

template <typename Params>
void NLWriter2<Params>::
WriteSparseEntry(File& nm_, int i, double x) {
  apr(nm_, "%d %g\n", i, x);
}

template <typename Params>
void NLWriter2<Params>::WriteSuffixes() {
  SuffixWriterFactory swf(*this);
  Feeder().FeedSuffixes(swf);
}

template <typename Params>
void NLWriter2<Params>::WritePLSOSConstraints() {
  PLSOSWriter plsos(*this);
  Feeder().FeedPLSOS(plsos);
}

template <typename Params>
void NLWriter2<Params>::WriteVarBounds() {
	apr(nm, "b\t#%d bounds (on variables)\n",
      (int)Hdr().num_vars);
  VarBndWriter vbw(*this);
  Feeder().FeedVarBounds(vbw);
  assert((int)Hdr().num_vars == vbw.GetNWritten());
}

template <typename Params>
void NLWriter2<Params>::WriteInitialGuesses() {
  SingleSparseDblVecWrtFactory
      vwf(*this, "x%d\t# initial guess\n");
  Feeder().FeedInitialGuesses(vwf);
}

template <typename Params>
void NLWriter2<Params>::WriteConBounds() {
  if (Hdr().num_algebraic_cons) {
		apr(nm, "r\t#%d ranges (rhs's)\n",  // and complementarity
        (int)Hdr().num_algebraic_cons);
    ConBndWriter cbw(*this);
    Feeder().FeedConBounds(cbw);
    assert((int)Hdr().num_algebraic_cons == cbw.GetNWritten());
  }
}

template <typename Params>
void NLWriter2<Params>::WriteDualInitialGuesses() {
  SingleSparseDblVecWrtFactory
      vwf(*this, "d%d\t# initial dual guess\n");
  Feeder().FeedInitialDualGuesses(vwf);
}

template <typename Params>
void NLWriter2<Params>::WriteBndRangeOrCompl(
    File& nm_,
    double L, double U, int k, int cvar) {
  if (k<=0) {          // normal algebraic constraint
    if (L <= NegInfty())
      apr(nm_, U >= Infty()
          ? "3\n" : "1 %.16g\n", U);
    else
      apr(nm_, U >= Infty() ? "2 %.16g\n"
                           : L == U ? "4 %.16g\n"
                                    : "0 %.16g %.16g\n",
          L, U);
  } else {             // complementarity
    apr(nm_, "5 %d %d\n", k, cvar + 1);   // add +1 here
  }
}


template <typename Params>
template <class ExprArgsFeeder>
void NLWriter2<Params>::WriteExprArgs(File& nm_,
                                      ExprArgsFeeder&& args) {
  while (args) {
    auto e = args.Next();
    WriteExpression(nm_, e);
  }
}

template <typename Params>
void NLWriter2<Params>::WriteDefinedVariables(int i) {
  DefVarWriterFactory dvwf(*this, i);
  Feeder().FeedDefinedVariables(i, dvwf);
}

template <typename Params>
void NLWriter2<Params>::WriteConObjExpressions() {
  int i=0;
  for (; i<Hdr().num_algebraic_cons; ++i) {
    WriteDefinedVariables(i+1);              // i+1
    apr(nm, "%c%d\t#%s\n", 'C', i,
        Feeder().ConDescription(i));
    ExprArgWriter ew(*this, 1);
    Feeder().FeedConExpression(i, ew);
  }
  for (;
       i<Hdr().num_algebraic_cons + Hdr().num_logical_cons;
       ++i) {
    WriteDefinedVariables(i+1);              // i+1
    apr(nm, "%c%d\t#%s\n", 'L',
        i - Hdr().num_algebraic_cons,
        Feeder().ConDescription(i));
    ExprArgWriter ew(*this, 1);
    Feeder().FeedConExpression(i, ew);
  }
  for (i=0; i<Hdr().num_objs; ++i) {
    WriteDefinedVariables(-i-1);              // -i-1
    apr(nm, "%c%d %d\t#%s\n", 'O', i,
        (int)Feeder().ObjType(i),
        Feeder().ObjDescription(i));
    ExprArgWriter ew(*this, 1);
    Feeder().FeedObjExpression(i, ew);
  }
}

template <typename Params>
void NLWriter2<Params>::WriteColumnSizes() {
  switch(Feeder().WantColumnSizes()) {
  case 1:
    apr(nm,
		#ifdef NL_LIB2_ORIG_HDR
				"k%d\t#intermediate Jacobian column lengths\n",
		#else
				"k%d\t#intermediate Jacobian column lengths (cumulative)\n",
		#endif
        Hdr().num_vars + Hdr().num_rand_vars - 1);
  {
    ColSizeWriter csw(*this, 1);
    Feeder().FeedColumnSizes(csw);
    assert(Hdr().num_vars + Hdr().num_rand_vars - 1
           == csw.GetNWritten());
  }
    break;
  case 2:
    apr(nm, "K%d\t#intermediate Jacobian column lengths\n",
        Hdr().num_vars + Hdr().num_rand_vars - 1);
  {
    ColSizeWriter csw(*this, 2);
    Feeder().FeedColumnSizes(csw);
    assert(Hdr().num_vars + Hdr().num_rand_vars - 1
           == csw.GetNWritten());
  }
    break;
  case 0:
    break;            // SKIP
  default:
    assert(0 && "why am I here?");
  }
}

template <typename Params>
void NLWriter2<Params>::WriteLinearConExpr() {
  for (int i=0; i<Hdr().num_algebraic_cons; ++i) {
    SingleSparseDblVecWrtFactory
        vwf(*this,
            [i, this](int nnz){
      this->apr(this->nm, "J%d %d\n", i, nnz);
    });
    Feeder().FeedLinearConExpr(i, vwf);
  }
}

template <typename Params>
void NLWriter2<Params>::WriteObjGradients() {
  for (int i=0; i<Hdr().num_objs; ++i) {
    SingleSparseDblVecWrtFactory
        vwf(*this,
            [i, this](int nnz){
      this->apr(this->nm, "G%d %d\n", i, nnz);
    });
    Feeder().FeedObjGradient(i, vwf);
  }
}

template <typename Params>
void NLWriter2<Params>::WriteRandomVariables() {
  RandVarWriterFactory rvwf(*this);
  Feeder().FeedRandomVariables(rvwf);
}


///////////////////// SparseVectorWriter //////////////////
template <typename Params>
template <class Index, class Value>
NLWriter2<Params>::SparseVectorWriter<Index, Value>::
SparseVectorWriter(NLWriter2& nlw, size_t n)
  : p_nlw_(&nlw), n_entries_(n) { }

template <typename Params>
template <class Index, class Value>
void NLWriter2<Params>::SparseVectorWriter<Index, Value>::
Write(Index i, Value v) {
  --n_entries_;
  p_nlw_->WriteSparseEntry(p_nlw_->nm, i, v);
}


///////////////////// ExprArgWriter ///////////////////////
template <typename Params>
NLWriter2<Params>::ExprArgWriter::
ExprArgWriter(NLWriter2& nlw, int na)
  : nlw_(nlw), nargs_(na) { assert(nargs_>0); }

template <typename Params>
NLWriter2<Params>::ExprArgWriter::
~ExprArgWriter()
{ assert(0==nargs_); }   // Check all written.

template <typename Params>
typename NLWriter2<Params>::ExprWriter
NLWriter2<Params>::ExprArgWriter::AddArg() {
  --nargs_;
  return {nlw_};
}

template <typename Params>
void NLWriter2<Params>::ExprArgWriter::EPut(
    typename FeederType::Expr e)
{ AddArg().EPut(e); }

template <typename Params>
void NLWriter2<Params>::ExprArgWriter::VPut(
    int v, const char* descr)
{ AddArg().VPut(v, descr); }

template <typename Params>
void NLWriter2<Params>::ExprArgWriter::NPut(
    double x) { AddArg().NPut(x); }

template <typename Params>
void NLWriter2<Params>::ExprArgWriter::StrPut(
    const char* s) { AddArg().StrPut(s); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprArgWriter::FuncPut(
    int index, int nArgs, const char* descr)
{ return AddArg().FuncPut(index, nArgs, descr); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprArgWriter::OPut1(
    int opcode, const char* descr)
{ return AddArg().OPut1(opcode, descr); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprArgWriter::OPut2(
    int opcode, const char* descr)
{ return AddArg().OPut2(opcode, descr); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprArgWriter::OPut3(
    int opcode, const char* descr)
{ return AddArg().OPut3(opcode, descr); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprArgWriter::OPutN(
    int opcode, int nArgs, const char* descr)
{ return AddArg().OPutN(opcode, nArgs, descr); }


///////////////////// ExprWriter ///////////////////////
template <typename Params>
void NLWriter2<Params>::ExprWriter::EPut(
    typename FeederType::Expr e) {
  ExprArgWriter ew(nlw_, 1);
  nlw_.Feeder().FeedExpr(e, ew);
}

template <typename Params>
void NLWriter2<Params>::ExprWriter::VPut(
    int v, const char* descr) {
  nlw_.apr(nlw_.nm,
           "v%d\t#%s\n", v, descr);
}

template <typename Params>
void NLWriter2<Params>::ExprWriter::NPut(
    double x) { nlw_.nput(nlw_.nm, x); }

template <typename Params>
void NLWriter2<Params>::ExprWriter::StrPut(
    const char* s) {
  nlw_.apr(nlw_.nm,
           "h%d:%s\n", (int)std::strlen(s), s);
}

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprWriter::FuncPut(
    int index, int nArgs, const char* descr) {
  nlw_.apr(nlw_.nm,
           "f%d %d\t#%s\n", index, nArgs, descr);
  return ExprArgWriter(nlw_, nArgs);
}

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprWriter::OPut1(
    int opcode, const char* descr) {
  nlw_.apr(nlw_.nm,
           "o%d\t#%s\n", opcode, descr);
  return ExprArgWriter(nlw_, 1);
}

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprWriter::OPut2(
    int opcode, const char* descr) {
  nlw_.apr(nlw_.nm,
           "o%d\t#%s\n", opcode, descr);
  return ExprArgWriter(nlw_, 2);
}

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprWriter::OPut3(
    int opcode, const char* descr) {
  nlw_.apr(nlw_.nm,
           "o%d\t#%s\n", opcode, descr);
  return ExprArgWriter(nlw_, 3);
}

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::ExprWriter::OPutN(
    int opcode, int nArgs, const char* descr) {
  nlw_.apr(nlw_.nm,
           "o%d\t#%s\n", opcode, descr);
  int n2write = nArgs;
  if (64 == opcode) {       // piecewise-linear
    assert(0 == n2write % 2);
    n2write /= 2;
  }
  nlw_.apr(nlw_.nm, "%d\n", n2write);
  return ExprArgWriter(nlw_, nArgs);
}


///////////////////// DEFINED VARIABLES ///////////////////
template <typename Params>
NLWriter2<Params>::DefVarWriter::
DefVarWriter(NLWriter2& nlw, int nnzlin)
  : nnzlin_(nnzlin), nlw_(nlw) { }

template <typename Params>
typename NLWriter2<Params>::SparseDblVecWriter
NLWriter2<Params>::DefVarWriter::GetLinExprWriter()
{ return SparseDblVecWriter(nlw_, nnzlin_); }

template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::DefVarWriter::GetExprWriter() {
  // Check that linw_ used? Still might not be written!
  return ExprArgWriter(nlw_, 1);
}

template <typename Params>
typename NLWriter2<Params>::DefVarWriter
NLWriter2<Params>::DefVarWriterFactory::
StartDefVar(int index, int nnz, const char* descr) {
  nlw_.apr(nlw_.nm,
           "V%d %d %d\t#%s\n",
           index, nnz,
           k_>=0 ? k_
                 : nlw_.Hdr().num_algebraic_cons
                   + nlw_.Hdr().num_logical_cons - k_,
           descr);
  return DefVarWriter(nlw_, nnz);
}


/////////////////////// BOUNDS ////////////////////////////
template <typename Params>
void
NLWriter2<Params>::VarBndWriter::WriteLbUb(
    double lb, double ub) {
  nlw_.WriteBndRangeOrCompl(nlw_.nm, lb, ub);
  ++nWrt_;
}

template <typename Params>
template <typename AlgConRange>
void
NLWriter2<Params>::ConBndWriter::WriteAlgConRange(
    AlgConRange acr) {
  nlw_.WriteBndRangeOrCompl(nlw_.nm,
                            acr.L, acr.U, acr.k, acr.cvar);
  ++nWrt_;
}


////////////////////// SUFFIXES ///////////////////////////
template <typename Params>
typename NLWriter2<Params>::SuffixIntWriter
NLWriter2<Params>::SuffixWriterFactory::StartIntSuffix(
    const char* name, int kind, int nnz) {
  assert(0 == (kind & 4));
  if (nnz)
    nlw_.apr(nlw_.nm, "S%d %d %s\n", kind, nnz, name);
  return SuffixIntWriter(nlw_, nnz);
}

template <typename Params>
typename NLWriter2<Params>::SuffixDblWriter
NLWriter2<Params>::SuffixWriterFactory::StartDblSuffix(
    const char* name, int kind, int nnz) {
  assert(0 != (kind & 4));
  if (nnz)
    nlw_.apr(nlw_.nm, "S%d %d %s\n", kind, nnz, name);
  return SuffixDblWriter(nlw_, nnz);
}


////////////////////// PL-SOS SUFFIXES /////////////
template <typename Params>
typename NLWriter2<Params>::SuffixIntWriter
NLWriter2<Params>::PLSOSWriter::StartSOSVars(int nnz)
{
  return swf_.StartIntSuffix("sos", 0, nnz);
}

template <typename Params>
typename NLWriter2<Params>::SuffixIntWriter
NLWriter2<Params>::PLSOSWriter::StartSOSCons(int nnz)
{
  return swf_.StartIntSuffix("sos", 1, nnz);
}

template <typename Params>
typename NLWriter2<Params>::SuffixDblWriter
NLWriter2<Params>::PLSOSWriter::StartSOSREFVars(int nnz)
{
  return swf_.StartDblSuffix("sosref", 4, nnz);
}


////////////////////// RANDOM VARS /////////////////
template <typename Params>
typename NLWriter2<Params>::ExprArgWriter
NLWriter2<Params>::RandVarWriterFactory::
StartRandVar(int index, const char* descr) {
  nlw_.apr(nlw_.nm,
           "R%d\t# %s\n", index, descr);
  return ExprArgWriter(nlw_, 1);
}


////////////////////// MATH ////////////////////////
template <typename Params>
double NLWriter2<Params>::Infty() const {
  return std::numeric_limits<double>::max();
}

template <typename Params>
double NLWriter2<Params>::NegInfty() const {
  return -Infty();
}

}  // namespace mp

#endif // NLWriter22_HPP
