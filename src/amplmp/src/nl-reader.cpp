/*
 .nl file support.

 Copyright (C) 2013 AMPL Optimization Inc

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

 Author: Victor Zverovich
 */

#include "mp/nl-reader.h"

namespace {
enum {
  USE_VBTOL_OPTION = 1,
  READ_VBTOL       = 3
};
}

void mp::ReadError::init(fmt::CStringRef filename, int line, int column,
                         fmt::CStringRef format_str, fmt::ArgList args) {
  filename_ = filename.c_str();
  line_ = line;
  column_ = column;
  fmt::MemoryWriter w;
  w.write("{}:{}:{}: ", filename, line, column);
  w.write(format_str, args);
  SetMessage(w.c_str());
}

mp::arith::Kind mp::arith::GetKind() {
  // Unlike ASL, we don't try detecting floating-point arithmetic at
  // configuration time because it doesn't work with cross-compiling.
  if (sizeof(double) != 2 * sizeof(uint32_t))
    return UNKNOWN;
  union {
    double d;
    uint32_t i[2];
  } u;
  u.i[0] = u.i[1] = 0;
  u.d = 1e13;
  if (u.i[1] == 0x42A2309C && u.i[0] == 0xE5400000)
    return IEEE_LITTLE_ENDIAN;
  if (u.i[0] == 0x42A2309C && u.i[1] == 0xE5400000)
    return IEEE_BIG_ENDIAN;
  return UNKNOWN;
}

fmt::Writer &mp::operator<<(fmt::Writer &w, const NLHeader &h) {
  w << (h.format == NLHeader::TEXT ? 'g' : 'b') << h.num_ampl_options;
  for (int i = 0; i < h.num_ampl_options; ++i)
    w << ' ' << h.ampl_options[i];
  if (h.ampl_options[USE_VBTOL_OPTION] == READ_VBTOL)
    w << ' ' << h.ampl_vbtol;
  w << '\n';
  w.write(" {} {} {} {} {} {}\n",
      h.num_vars, h.num_algebraic_cons, h.num_objs,
      h.num_ranges, h.num_eqns, h.num_logical_cons);
  w.write(" {} {} {} {} {} {}\n",
      h.num_nl_cons, h.num_nl_objs,
      h.num_compl_conds - h.num_nl_compl_conds,
      h.num_nl_compl_conds, h.num_compl_dbl_ineqs,
      h.num_compl_vars_with_nz_lb);
  w.write(" {} {}\n", h.num_nl_net_cons, h.num_linear_net_cons);
  w.write(" {} {} {}\n",
      h.num_nl_vars_in_cons, h.num_nl_vars_in_objs, h.num_nl_vars_in_both);
  w.write(" {} {} {} {}\n",
      h.num_linear_net_vars, h.num_funcs,
          h.format == NLHeader::TEXT ? 0 : h.arith_kind, h.flags);
  w.write(" {} {} {} {} {}\n",
      h.num_linear_binary_vars, h.num_linear_integer_vars,
      h.num_nl_integer_vars_in_both, h.num_nl_integer_vars_in_cons,
      h.num_nl_integer_vars_in_objs);
  w.write(" {} {}\n", h.num_con_nonzeros, h.num_obj_nonzeros);
  w.write(" {} {}\n", h.max_con_name_len, h.max_var_name_len);
  w.write(" {} {} {} {} {}\n",
      h.num_common_exprs_in_both, h.num_common_exprs_in_cons,
      h.num_common_exprs_in_objs, h.num_common_exprs_in_single_cons,
      h.num_common_exprs_in_single_objs);
  return w;
}

mp::internal::ReaderBase::ReaderBase(NLStringRef data, fmt::CStringRef name)
: ptr_(data.c_str()), start_(ptr_), end_(ptr_ + data.size()),
  token_(ptr_), name_(name.c_str()) {}

template <typename Locale>
mp::internal::TextReader<Locale>::TextReader(
    NLStringRef data, fmt::CStringRef name)
: ReaderBase(data, name), line_start_(ptr_), line_(1) {}

template <typename Locale>
void mp::internal::TextReader<Locale>::DoReportError(
    const char *loc, fmt::CStringRef format_str, const fmt::ArgList &args) {
  int line = line_;
  const char *line_start = line_start_;
  if (loc < line_start) {
    --line;
    // Find the beginning of the previous line.
    line_start = loc;
    if (*line_start == '\n')
      --line_start;
    while (*line_start != '\n' && line_start != start_)
      --line_start;
    if (*line_start == '\n')
      ++line_start;
  }
  int column = static_cast<int>(loc - line_start + 1);
  throw ReadError(name_, line, column, format_str, args);
}

template <typename Locale>
bool mp::internal::TextReader<Locale>::ReadOptionalDouble(double &value) {
  SkipSpace();
  if (*ptr_ == '\n')
    return false;
  char *end = 0;
  value = std::strtod(ptr_, &end);
  bool has_value = ptr_ != end;
  ptr_ = end;
  return has_value;
}

template <typename Locale>
fmt::StringRef mp::internal::TextReader<Locale>::ReadString() {
  int length = ReadUInt();
  if (*ptr_ != ':')
    DoReportError(ptr_, "expected ':'");
  ++ptr_;
  const char *start = ptr_;
  for (int i = 0; i < length; ++i, ++ptr_) {
    char c = *ptr_;
    if (c == '\n') {
      line_start_ = ptr_  + 1;
      ++line_;
    } else if (!c && ptr_ == end_) {
      DoReportError(ptr_, "unexpected end of file in string");
    }
  }
  if (*ptr_ != '\n')
    DoReportError(ptr_, "expected newline");
  ++line_;
  line_start_ = ++ptr_;
  return fmt::StringRef(length != 0 ? start : 0, length);
}

template <typename Locale>
fmt::StringRef mp::internal::TextReader<Locale>::ReadName() {
  SkipSpace();
  const char *start = ptr_;
  if (*ptr_ == '\n' || !*ptr_)
    ReportError("expected name");
  do ++ptr_;
  while (!std::isspace(*ptr_) && *ptr_);
  return fmt::StringRef(start, ptr_ - start);
}

template <typename Locale>
void mp::internal::TextReader<Locale>::ReadHeader(NLHeader &header) {
  // Read the format (text or binary).
  switch (ReadChar()) {
  case 'g':
    header.format = NLHeader::TEXT;
    break;
  case 'b':
    header.format = NLHeader::BINARY;
    break;
  default:
    ReportError("expected format specifier");
    break;
  }

  // Read options.
  ReadOptionalUInt(header.num_ampl_options);
  if (header.num_ampl_options > MAX_AMPL_OPTIONS)
    ReportError("too many options");
  for (int i = 0; i < header.num_ampl_options; ++i) {
    double tmp;
    if (!ReadOptionalDouble(tmp))
      break;
    header.ampl_options[i] = (long)tmp;
    if (header.ampl_options[i] != tmp)
      break;
  }
  if (header.ampl_options[USE_VBTOL_OPTION] == READ_VBTOL)
    ReadOptionalDouble(header.ampl_vbtol);
  ReadTillEndOfLine();

  // Read problem dimensions.
  header.num_vars = ReadUInt();
  header.num_algebraic_cons = ReadUInt();
  header.num_objs = ReadUInt();
  header.num_eqns = -1;
  if (ReadOptionalUInt(header.num_ranges) &&
      ReadOptionalUInt(header.num_eqns)) {
      ReadOptionalUInt(header.num_logical_cons);
  }
  ReadTillEndOfLine();

  // Read the nonlinear and complementarity information.
  header.num_nl_cons = ReadUInt();
  header.num_nl_objs = ReadUInt();
  bool all_compl =
      ReadOptionalUInt(header.num_compl_conds) &&
      ReadOptionalUInt(header.num_nl_compl_conds) &&
      ReadOptionalUInt(header.num_compl_dbl_ineqs) &&
      ReadOptionalUInt(header.num_compl_vars_with_nz_lb);
  header.num_compl_conds += header.num_nl_compl_conds;
  if (header.num_compl_conds > 0 && !all_compl)
    header.num_compl_dbl_ineqs = -1;
  ReadTillEndOfLine();

  // Read the information about network constraints.
  header.num_nl_net_cons = ReadUInt();
  header.num_linear_net_cons = ReadUInt();
  ReadTillEndOfLine();

  // Read the information about nonlinear variables.
  header.num_nl_vars_in_cons = ReadUInt();
  header.num_nl_vars_in_objs = ReadUInt();
  header.num_nl_vars_in_both = -1;
  ReadOptionalUInt(header.num_nl_vars_in_both);
  ReadTillEndOfLine();

  header.num_linear_net_vars = ReadUInt();
  header.num_funcs = ReadUInt();
  int arith_kind = 0;
  if (ReadOptionalUInt(arith_kind)) {
    if (arith_kind > arith::LAST)
      ReportError("unknown floating-point arithmetic kind");
    header.arith_kind = arith_kind;
    ReadOptionalUInt(header.flags);
  }
  ReadTillEndOfLine();

  // Read the information about discrete variables.
  header.num_linear_binary_vars = ReadUInt();
  header.num_linear_integer_vars = ReadUInt();
  if (header.num_nl_vars_in_both >= 0) {  // ampl versions >= 19930630
    header.num_nl_integer_vars_in_both = ReadUInt();
    header.num_nl_integer_vars_in_cons = ReadUInt();
    header.num_nl_integer_vars_in_objs = ReadUInt();
  }
  ReadTillEndOfLine();

  // Read the information about nonzeros.
  header.num_con_nonzeros = ReadUInt<std::size_t>();
  header.num_obj_nonzeros = ReadUInt<std::size_t>();
  ReadTillEndOfLine();

  // Read the information about names.
  header.max_con_name_len = ReadUInt();
  header.max_var_name_len = ReadUInt();
  ReadTillEndOfLine();

  // Read the information about common expressions checking for overflow
  // as the variable indices go from 0 to num_vars + num_common_exprs.
  int max_vars = header.num_vars;
  header.num_common_exprs_in_both = ReadUInt(max_vars);
  header.num_common_exprs_in_cons = ReadUInt(max_vars);
  header.num_common_exprs_in_objs = ReadUInt(max_vars);
  header.num_common_exprs_in_single_cons = ReadUInt(max_vars);
  header.num_common_exprs_in_single_objs = ReadUInt(max_vars);
  ReadTillEndOfLine();
}

void mp::internal::BinaryReaderBase::ReportError(
    fmt::CStringRef format_str, const fmt::ArgList &args) {
  fmt::MemoryWriter w;
  std::size_t offset = token_ - start_;
  w.write("{}:offset {}: ", name_, offset);
  w.write(format_str, args);
  throw BinaryReadError(name_, offset, w.c_str());
}

template <typename File>
void mp::internal::NLFileReader<File>::Open(fmt::CStringRef filename) {
  file_ = File(filename, fmt::File::RDONLY | fmt::File::BINARY);
  size_ = ConvertFileToMmapSize(file_.size(), filename);
  // Round size up to a multiple of page_size. The remainded of the last
  // partial page is zero-filled both on POSIX and Windows so the resulting
  // memory buffer is zero terminated.
  std::size_t page_size = fmt::getpagesize();
  std::size_t remainder = size_ % page_size;
  rounded_size_ = remainder != 0 ? (size_ + page_size - remainder) : size_;
}

template <typename File>
void mp::internal::NLFileReader<File>::Read(
    fmt::internal::MemoryBuffer<char, 1> &array) {
  array.resize(size_ + 1);
  std::size_t offset = 0;
  while (offset < size_)
    offset += file_.read(&array[offset], size_ - offset);
  array[size_] = 0;
}


namespace internal {

class NameHandler {
 private:
  std::vector<const char *> &names_;
  fmt::StringRef name_;

 public:
  explicit NameHandler(std::vector<const char *> &names)
    : names_(names), name_("") {}

  void OnName(fmt::StringRef name) {
    name_ = name;
    names_.push_back(name.data());
  }

  fmt::StringRef name() const { return name_; }
};

}  // namespace internal


mp::NameProvider::NameProvider(
    fmt::CStringRef filename, fmt::CStringRef gen_name, std::size_t num_items)
  : gen_name_(gen_name.c_str()) {
  ReadNames(filename, num_items);
}

mp::NameProvider::NameProvider(
    fmt::CStringRef gen_name, fmt::CStringRef n2)
  : gen_name_(gen_name.c_str()), gen_name_2_(n2.c_str()) { }

void mp::NameProvider::ReadNames(
    fmt::CStringRef filename, std::size_t num_items) {
  ::internal::NameHandler handler(names_);
  names_.reserve(num_items + 1);
  try {
    reader_.Read(filename, handler);
  } catch (const fmt::SystemError &) {
    return; // System error, ignore the file and use generated names.
  }
  fmt::StringRef last_name = handler.name();
  names_.push_back(last_name.data() + last_name.size() + 1);
}

fmt::StringRef mp::NameProvider::name(
    std::size_t index, std::size_t i2) {
  if (index + 1 < names_.size()) {
    const char *name = names_[index];
    const auto* pos1past = names_[index + 1] - 1;
    assert( ('\n' == *pos1past) || ('\r' == *pos1past));
    if (pos1past>name && '\r' == *(pos1past-1))            // Windows
      --pos1past;
    return fmt::StringRef(name, pos1past - name);
  }
  writer_.clear();
  if (index>=i2)
    writer_ << gen_name_2_ << '[' << (index - i2 + 1) << ']';
  else
    writer_ << gen_name_ << '[' << (index + 1) << ']';
  return fmt::StringRef(writer_.c_str(), writer_.size());
}

size_t mp::NameProvider::number_read() const {
  return (names_.size() ? names_.size()-1 : 0);
}

std::vector<std::string> mp::NameProvider::get_names(
    size_t n, size_t i2) {
  std::vector<std::string> result;
  result.reserve(n);
  for (size_t i=0; i<n; ++i)
    result.push_back(name(i, i2));
  return result;
}

// Instantiate
template class mp::internal::TextReader<>;
template class mp::internal::NLFileReader<>;
