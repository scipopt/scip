# This Makefile enables the automatic generation of macro definition
# headers, --i options and size options for Lint based on command-line
# switches passed to GCC.
#
# Usage:
#
# 	make -f co-gcc.mak \
#		GCC_BIN='name of the gcc binary' \
#		GXX_BIN='name of the g++ binary' \
# 		CFLAGS='[usual C compile switches here]' \
# 		CXXFLAGS='[usual C++ compile switches here]' \
# 		CPPFLAGS='[usual common preprocessor switches here]' \
# 		COMMON_FLAGS='[usual C & C++ compile switches here]' \
#
# ... where 'make' is the name of the GNU Make program on your system.
# That invocation should generate the following files:
#
#       lint_cmac.h
#       lint_cppmac.h
#       gcc-include-path.lnt
#       size-options.lnt
#
# Note, if you do not supply the options that you actually compile with,
# you may see undesired results.  Examples:
#
# 1) If you usually compile with -m64 but do not pass this in the
# COMMON_FLAGS variable when you run `make -f co-gcc.mak` then Lint may
# see the wrong size options (so it may think e.g. that sizeof(void*) is
# 4, which of course is inappropriate if you compile your code in 64-bit
# mode).
#
# 2) The set of compile switches (even non-preprocessor switches like -O3)
# can affect the configuration of GCC's preprocessor, which means it can
# affect how the preprocessor views the contents of system headers (and
# hence the token sequence it generates).  So if we don't see the right
# set of compile switches here then the program you Lint might not be the
# program you compile (even though the same .c and .cpp files are
# involved).
#
# See also the file gcc-readme.txt (supplied with the Lint distribution).

COMMON_GCC_OPTS:= $(COMMON_FLAGS) $(CPPFLAGS) -Wno-long-long
# We want to enable 'long long' for the purpose of extracting the value of
# 'sizeof(long long)'; see the 'sizes' target below.

C_OPTS:=  $(CFLAGS) $(COMMON_GCC_OPTS)
CXX_OPTS:=$(CXXFLAGS) $(COMMON_GCC_OPTS)
# Note, we're not *yet* able to handle some of the header contents when
# -std=c++0x is given.

GCC_BIN:=gcc
GXX_BIN:=g++

GCC:=$(GCC_BIN) $(C_OPTS)
GXX:=$(GXX_BIN) $(CXX_OPTS)

TEMP_FILE_PREFIX:=co-gcc.mak.temp

E:=$(TEMP_FILE_PREFIX)-empty
SIZE_GEN:=$(TEMP_FILE_PREFIX)-generate-size-options

ECHO:=echo
TOUCH:=touch
AWK:=awk

.PHONY = clean clean_temps

config: clean preprocessor sizes clean_temps

preprocessor: empty_files macros include_path

empty_files:
	$(RM) $(E)*
	$(TOUCH) $(E).cpp $(E).c

macros:
	$(GCC) -E -dM $(E).c   -o lint_cmac.h
	$(GXX) -E -dM $(E).cpp -o lint_cppmac.h

include_path:
	@# Here we make options for the #include search path.
	@# Note, frameworks (a feature of Apple's GCC) are not supported
	@# yet so for now we filter them out.  Each remaining search
	@# directory 'foo' is transformed into '--i"foo"' after
	@# superfluous directory separators  are removed (as well as each
	@# CR character appearing immediately before a newline):
	$(GXX) -v -c $(E).cpp 2>&1 \
	| $(AWK) '					\
	    BEGIN  {S=0}				\
	    /search starts here:/  {S=1;next;}		\
	    S && /Library\/Frameworks/ {next;}		\
	    S && /^ /  {				\
		sub("^ ","");				\
		gsub("//*","/");			\
		sub("\xd$$","");			\
		sub("/$$","");				\
		printf("--i\"%s\"\n", $$0);		\
		next;					\
	    }						\
	    S  {exit;}					\
	    ' >gcc-include-path.lnt
	@# Note, we deliberately use '--i' instead of '-i' here; the effect
	@# is that the directories named with the double-dash form are
	@# searched after directories named with the single-dash form.
	@# (See also the entry for '--i' in section 5.7 of the Lint
	@# manual.)
	@#
	@# We typically use '--i' when we want to name a system include
	@# directory, which GCC searches only after it searches all
	@# directories named in a '-I' option.  The upshot is that the
	@# correct search order (i.e., project includes before system
	@# includes) is preserved even when double-dash-i options are given
	@# before single-dash-i options.
	@#
	@# Also note, no harm is done if '-I' options are passed to GCC
	@# here:  directories named with '-I' will appear before the
	@# sys-include-dirs in GCC's output, so even though Lint might then
	@# see a project-include directory named with a '--i' option, that
	@# directory will still be searched before the sys-includes because
	@# of the ordering of '--i' options.  (Just make sure you don't use
	@# the double-dash form with project include dirs outside of this
	@# limited & generated sub-sequence of options because this is the
	@# only place where we are certain that project directories always
	@# come before system directories.)
	@#
	@# XXX:  We need to do something for people without an AWK
	@# implementation---hopefully not here, but perhaps we could
	@# provide an 'awk' with the Lint distro.


sizes:
	$(RM) $(SIZE_GEN)*
	@# 'echo' seems to vary in behavior with respect to its handling
	@# of '\n'.  (Is it a newline, or a literal backslash followed by
	@# a literal 'n'?  It seems to depend on your platform.)  So we
	@# deliberately avoid the use of explicit newline characters here.
	@$(ECHO) '\
extern  "C" int printf(const char*, ...);\
int main() {\
printf( "-ss%lu  ", sizeof(short) );\
printf( "-si%lu  ", sizeof(int) );\
printf( "-sl%lu  ", sizeof(long) );\
printf( "-sll%lu  ", sizeof(long long) );\
printf( "-sf%lu  ", sizeof(float) );\
printf( "-sd%lu  ", sizeof(double) );\
printf( "-sld%lu  ", sizeof(long double) );\
printf( "-sp%lu  ", sizeof(void*) );\
printf( "-sw%lu  ", sizeof(wchar_t) );\
}' >$(SIZE_GEN).cc
	$(GXX) $(SIZE_GEN).cc -o $(SIZE_GEN)
	./$(SIZE_GEN) >size-options.lnt
	@# ... and make it newline-terminated:
	@$(ECHO) ""  >>size-options.lnt

clean_temps:
	$(RM) $(TEMP_FILE_PREFIX)*

clean:
	$(RM) \
	    lint_cppmac.h \
	    lint_cmac.h \
	    gcc-include-path.lnt \
	    size-options.lnt

# vim:ts=8
