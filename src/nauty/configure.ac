# Process this file with autoconf to produce a configure script.
AC_INIT([nauty],[2.8.6])
AC_CONFIG_SRCDIR([nauty-h.in])

# ======================================================================
# First we define some functions

dnl AC_TRY_CFLAGS(CFLAGS, [ACTION-IF-WORKS], [ACTION-IF-FAILS])
dnl Check if $CC supports a given set of cflags
AC_DEFUN([AC_TRY_CFLAGS],
    [AC_MSG_CHECKING([if $CC supports $1 flags])
    SAVE_CFLAGS="$CFLAGS"
    CFLAGS="$1"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], \
                [[]])],[ac_cv_try_cflags_ok=yes],[ac_cv_try_cflags_ok=no])
    CFLAGS="$SAVE_CFLAGS"
    AC_MSG_RESULT([$ac_cv_try_cflags_ok])
    if test "x$ac_cv_try_cflags_ok" = xyes; then
        ifelse([$2],[],[:],[$2])
    else
        ifelse([$3],[],[:],[$3])
    fi])

dnl CVT_YESNO([YES-NO VARIABLE],[0-1 VARIABLE])
dnl Make a 0-1 output variable from a yes/no shell variable
AC_DEFUN([CVT_YESNO],[if test "x$$1" = xyes; then
  $2=1
else
  $2=0
fi
AC_SUBST($2)])

dnl Check for thread-local attribute
AC_DEFUN([AX_TLS], [
  AC_MSG_CHECKING([for thread local storage (TLS) class])
  AC_CACHE_VAL([ac_cv_tls],
   [for ax_tls_keyword in thread_local _Thread_local __thread '__declspec(thread)' none; do
       AS_CASE([$ax_tls_keyword],
          [none], [ac_cv_tls=none ; break],
          [AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
              [#include <stdlib.h>],
              [static  $ax_tls_keyword  int bar;]
            )],
            [ac_cv_tls=$ax_tls_keyword ; break],
            [ac_cv_tls=none]
          )]
        )
    done ]
  )
  AC_MSG_RESULT([$ac_cv_tls])

  AS_IF([test "$ac_cv_tls" != "none"],
    [AC_DEFINE_UNQUOTED([TLS],[$ac_cv_tls],[If the compiler supports a TLS storage class, define it to that here])
     m4_ifnblank([$1],[$1],[[:]])],
    [m4_ifnblank([$2],[$2],[[:]])])
])

# ======================================================================
# First we check for configuration switches

dnl --disable-popcnt disallows popcount instruction, otherwise they are tested for
AC_ARG_ENABLE([popcnt],
   AS_HELP_STRING([--disable-popcnt], [Disable popcnt extensions]),
   [enablearg=given], [enablearg=notgiven])
AS_IF([test $enablearg = notgiven || test "x$enable_popcnt" = xyes],[allow_popcnt=1],[allow_popcnt=0])

dnl --disable-clz disallows clz instructions, otherwise they are tested for
AC_ARG_ENABLE([clz],
   AS_HELP_STRING([--disable-clz], [Disable clz extensions]),
   [enablearg=given], [enablearg=notgiven])
AS_IF([test $enablearg = notgiven || test "x$enable_clz" = xyes],[allow_clz=1],[allow_clz=0])

dnl --enable-generic disallows -march=native, otherwise it is tested for
AC_ARG_ENABLE([generic],
   AS_HELP_STRING([--enable-generic], [Disable -march=native switch]),
   [enablearg=given], [enablearg=notgiven])
AS_IF([test $enablearg = notgiven || test "x$enable_generic" = xno],[allow_native=1],[allow_native=0])

dnl --enable-tls specifies thread-local storage, if possible
AC_ARG_ENABLE([tls],
   AS_HELP_STRING([--enable-tls], [Enable thread-local storage]),
   [use_tls="#define USE_TLS"], [use_tls=])

dnl --enable-ansi specifies ANSI controls are allowed in output
AC_ARG_ENABLE([ansi],
   AS_HELP_STRING([--enable-ansi], [Enable ANSI controls in output]),
   [enablearg=given], [enablearg=notgiven])
AS_IF([test $enablearg = notgiven || test "x$enable_ansi" = xno],[allow_ansi=0],[allow_ansi=1])

dnl --disable-interrupt prevents nauty or traces catching control-C, otherwise allow it
AC_ARG_ENABLE([interrupt],
   AS_HELP_STRING([--disable-interrupt], [Disable interrupt catching in nauty&traces]),
   [enablearg=given], [enablearg=notgiven])
AS_IF([test $enablearg = notgiven || test "x$enable_interrupt" = xyes],[allow_interrupt=1],[allow_interrupt=0])

dnl --enable-wordsize=value overrides default rule for choosing WORDSIZE
AC_ARG_ENABLE([wordsize],
   AS_HELP_STRING([--enable-wordsize=value], [Override default WORDSIZE choice]))
AS_IF([test "x$enable_wordsize" = x16],[default_wordsize=16],[default_wordsize=0])
AS_IF([test "x$enable_wordsize" = x32],[default_wordsize=32])
AS_IF([test "x$enable_wordsize" = x64],[default_wordsize=64])

# ======================================================================
# Next we get to work

dnl Checks for system features
AC_CANONICAL_HOST

case "$host" in
  *powerpc*)
     is_powerpc=1 ;;
  *)
     is_powerpc=0 ;;
esac

case "$host" in
  aarch64*)
     is_ARM64=1 ;;
  *)
     is_ARM64=0 ;;
esac
AC_SUBST(is_ARM64)

case "$build" in
  *cygwin*)
     is_cygwin=1 ;;
  *)
     is_cygwin=0 ;;
esac

dnl Checks for C compiler and sets CFLAGS if not set by user
user_cflags="$CFLAGS"
AC_PROG_CC
CFLAGS=$user_cflags
MORECFLAGS=""
dnl  we need  AC_SYS_LARGEFILE and AC_FUNC_FSEEKO
AC_SYS_LARGEFILE
AS_IF([test "x$ac_cv_sys_file_offset_bits" = xno],[ac_cv_sys_file_offset_bits=0])
AC_SUBST(ac_cv_sys_file_offset_bits)

AS_IF([test "x$user_cflags" = x || test "x$user_cflags" = x-m32],
[
  AS_IF([test "$CC" = "icc"],
    [AC_TRY_CFLAGS([-O3],[CFLAGS="$CFLAGS -O3"]);
     AS_IF([test "x$ac_cv_try_cflags_ok" = xno],
        [AC_TRY_CFLAGS([-O2],[CFLAGS="$CFLAGS -O2"],[CFLAGS="$CFLAGS -O"])])],
    [AC_TRY_CFLAGS([-O4 -Werror],[CFLAGS="$CFLAGS -O4"]);
     AS_IF([test "x$ac_cv_try_cflags_ok" = xno],
       [AC_TRY_CFLAGS([-O3],[CFLAGS="$CFLAGS -O3"])
        AS_IF([test "x$ac_cv_try_cflags_ok" = xno],
          [AC_TRY_CFLAGS([-O2],[CFLAGS="$CFLAGS -O2"],[CFLAGS="$CFLAGS -O"])])
       ])
    ])
])
 
case "$CC" in
  gcc*)
     gccver=`"$CC" -dumpversion` ;;
  *)
     gccver=0.0.0 ;;
esac

dnl Checks relevant to Cygwin and other Windows things
AC_EXEEXT

dnl Checks for header files.

AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS(signal.h stddef.h unistd.h sys/types.h sys/wait.h stdlib.h string.h errno.h limits.h stdint.h)
CVT_YESNO(ac_cv_header_stddef_h,header_stddef_h)
CVT_YESNO(ac_cv_header_unistd_h,header_unistd_h)
CVT_YESNO(ac_cv_header_stdlib_h,header_stdlib_h)
CVT_YESNO(ac_cv_header_string_h,header_string_h)
CVT_YESNO(ac_cv_header_sys_wait_h,header_sys_wait_h)
CVT_YESNO(ac_cv_header_sys_types_h,header_sys_types_h)
CVT_YESNO(ac_cv_header_errno_h,header_errno_h)
CVT_YESNO(ac_cv_header_signal_h,header_signal_h)
CVT_YESNO(ac_cv_header_limits_h,header_limits_h)
CVT_YESNO(ac_cv_header_stdint_h,header_stdint_h)

AC_TYPE_PID_T
CVT_YESNO(ac_cv_type_pid_t,have_pid_t)

AC_CHECK_DECL([ftell],ftell_dec=1,ftell_dec=0)
AC_SUBST(ftell_dec)
AC_CHECK_DECL([popen],popen_dec=1,popen_dec=0)
AC_SUBST(popen_dec)
AC_CHECK_DECL([fdopen],fdopen_dec=1,fdopen_dec=0)
AC_SUBST(fdopen_dec)
AC_CHECK_DECL([putenv],putenv_dec=1,putenv_dec=0)
AC_SUBST(putenv_dec)
AC_CHECK_DECL([setenv],setenv_dec=1,setenv_dec=0)
AC_SUBST(setenv_dec)

AC_CHECK_DECL([malloc],malloc_dec=1,malloc_dec=0)
AS_IF([test $malloc_dec -eq 0],
      [AC_EGREP_HEADER([(^|[^a-zA-Z_0-9])malloc[^a-zA-Z_0-9]],malloc.h,
       malloc_dec=2,malloc_dec=0)])
AC_SUBST(malloc_dec)

AC_MSG_CHECKING(if compiler accepts flexible array members)
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[]], \
  [[void foo(void*);struct { int a; int b[]; } ab; foo(&ab); ]])],
[flex_array_ok=1],[flex_array_ok=0])
AC_SUBST(flex_array_ok)
AC_MSG_RESULT($flex_array_ok)

AC_MSG_CHECKING(if INFINITY is declared in math.h)
AC_EGREP_CPP(yes,
[#include <stdio.h>
#include <math.h>
#ifdef INFINITY
  yes
#endif
], has_math_inf=1, has_math_inf=0)
AC_SUBST(has_math_inf)
AC_MSG_RESULT($has_math_inf)

dnl Checks for sizes of integer types; avoid 64-bit if necessary
AC_CHECK_SIZEOF([int])
AC_SUBST(ac_cv_sizeof_int)
AC_CHECK_SIZEOF([long])
AC_SUBST(ac_cv_sizeof_long)
AC_CHECK_SIZEOF([long long])
AC_SUBST(ac_cv_sizeof_long_long)
AC_CHECK_SIZEOF([void*])
AC_SUBST(ac_cv_sizeof_voidp)
AC_CHECK_SIZEOF([__int128])
AC_SUBST(ac_cv_sizeof___int128)
AC_CHECK_SIZEOF([__int128_t])
AC_SUBST(ac_cv_sizeof___int128_t)
lok=0
testprogs="dreadtest dreadtestS dreadtestS1 dreadtest4K"
testprogs="$testprogs dreadtest1 dreadtestW dreadtestW1"

AS_IF([test $ac_cv_sizeof_long_long -eq 8],
      [lok=1
       testprogs="$testprogs dreadtestL1 dreadtestL"],
      [AS_IF([test $ac_cv_sizeof_long -eq 8],
             [lok=1])
    ])
AC_SUBST(lok)
AC_SUBST(testprogs)

dnl Checks for library functions.
AC_CHECK_FUNC(isatty,have_isatty=1,have_isatty=0)
AC_SUBST(have_isatty)
AC_CHECK_FUNC(time,have_time=1,have_time=0)
AC_SUBST(have_time)
AC_CHECK_FUNC(gettimeofday,have_gettimeofday=1,have_gettimeofday=0)
AC_SUBST(have_gettimeofday)
AC_CHECK_FUNC(clock_gettime,have_clock_gettime=1,have_clock_gettime=0)
AC_SUBST(have_clock_gettime)
AC_CHECK_FUNC(clock,have_clock=1,have_clock=0)
AC_SUBST(have_clock)
AC_CHECK_FUNC(times,have_times=1,have_times=0)
AC_SUBST(have_times)
AC_CHECK_FUNC(getrusage,have_getrusage=1,have_getrusage=0)
AC_SUBST(have_getrusage)
AC_CHECK_FUNC(perror,have_perror=1,have_perror=0)
AC_SUBST(have_perror)
AC_CHECK_FUNC(pipe,have_pipe=1,have_pipe=0)
AC_SUBST(have_pipe)
AC_CHECK_FUNC(wait,have_wait=1,have_wait=0)
AC_SUBST(have_wait)
AC_CHECK_FUNC(popen,have_popen=1,have_popen=0)
AC_SUBST(have_popen)
AC_CHECK_FUNC(putenv,have_putenv=1,have_putenv=0)
AC_SUBST(have_putenv)
AC_CHECK_FUNC(setenv,have_setenv=1,have_setenv=0)
AC_SUBST(have_setenv)
AC_FUNC_FORK
CVT_YESNO(ac_cv_func_fork,have_fork)
AC_SUBST(have_fork)
AC_CHECK_FUNC(fseeko,have_fseeko=1,have_fseeko=0)
AC_SUBST(have_fseeko)
AC_CHECK_FUNC(sigaction,have_sigaction=1,have_sigaction=0)
AC_SUBST(have_sigaction)
AC_CHECK_FUNC(sigprocmask,have_sigprocmask=1,have_sigprocmask=0)
AC_SUBST(have_sigprocmask)

dnl --check if popcnt instruction is available and desired
AC_MSG_CHECKING(if popcnt instruction is available and requested)
AC_RUN_IFELSE([AC_LANG_PROGRAM([],[[if (__builtin_cpu_supports("popcnt")) return 0; else return 1;]])],
  popsup=1,popsup=0)

AS_IF([test "$allow_popcnt" -eq 1],
      [have_hwpopcnt=$popsup
test -f /proc/cpuinfo && grep -i popcnt /proc/cpuinfo >/dev/null && have_hwpopcnt=1 
test $have_hwpopcnt -eq 0 && sysctl machdep.cpu.features 2>/dev/null | grep -i popcnt >/dev/null && have_hwpopcnt=1
test $have_hwpopcnt -eq 0 && isainfo -x -v 2>/dev/null | grep -i popc >/dev/null && have_hwpopcnt=1
test $have_hwpopcnt -eq 0 && test $is_ARM64 -eq 1 && have_hwpopcnt=1
],
      [have_hwpopcnt=0
])
# Some problems with popcnt on powerpc are not resolved yet
test $is_powerpc -eq 1 && have_hwpopcnt=0
AC_MSG_RESULT($have_hwpopcnt)
AC_SUBST(have_hwpopcnt)

SAVE_CFLAGS="$CFLAGS"
have_popcnt=0
have_popcntl=0
have_popcntll=0
have_mmpop32=0
have_mmpop64=0

AS_IF([test $have_hwpopcnt -eq 1],
      [AS_IF([test "$CC" = "icc"],
             [AC_MSG_CHECKING(if _mm_popcnt_u32() is supported)
              AC_LINK_IFELSE([AC_LANG_SOURCE([#include <nmmintrin.h>
                   int main(){unsigned int x; x = _mm_popcnt_u32(x);}])],
                             [have_mmpop32=1], [have_mmpop32=0])
              AC_MSG_RESULT($have_mmpop32)
              AC_MSG_CHECKING(if _mm_popcnt_u64() is supported)
              AC_LINK_IFELSE([AC_LANG_SOURCE([#include <nmmintrin.h>
                   int main(){unsigned int x; x = _mm_popcnt_u64(x);}])],
                             [have_mmpop64=1], [have_mmpop64=0])
              AC_MSG_RESULT($have_mmpop64)],
             [AS_IF([test "$is_ARM64" -eq 0], [CFLAGS="$CFLAGS -mpopcnt"])
              AC_MSG_CHECKING(if __builtin_popcount() is supported)
              AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned int x; x = __builtin_popcount(x);}])],
                             [have_popcnt=1], [have_popcnt=0])
              AC_MSG_RESULT($have_popcnt)
              AC_MSG_CHECKING(if __builtin_popcountl() is supported)
              AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned long x; x = __builtin_popcountl(x);}])],
                             [have_popcntl=1], [have_popcntl=0])
              AC_MSG_RESULT($have_popcntl)
              AC_MSG_CHECKING(if __builtin_popcountll() is supported)
              AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned long long x; x = __builtin_popcountll(x);}])],
                             [have_popcntll=1], [have_popcntll=0])
              AC_MSG_RESULT($have_popcntll)
              AS_IF([test $is_ARM64 -eq 1], [],
                    [test $have_popcnt -eq 1 || test $have_popcntl -eq 1 || test $have_popcntll -eq 1],
                    [MORECFLAGS="$MORECFLAGS -mpopcnt"])
])])
CFLAGS="$SAVE_CFLAGS"
AC_SUBST(have_popcnt)
AC_SUBST(have_popcntl)
AC_SUBST(have_popcntll)
AC_SUBST(have_mmpop32)
AC_SUBST(have_mmpop64)

# Unless --enable-generic is given to configure, we try to use the -march=native compiler
#   switch.  However, it breaks some versions of gcc on MacOSX due to a compiler bug, so
#   we will also add --mno-avx if that is necessary.
# But don't add -march-native if the user provided an architecture via CFLAGS or CC

case "$CFLAGS" in
 *-march=* ) allow_native=0 ;;
esac
case "$CC" in
 *-march=* ) allow_native=0 ;;
esac
case "$build" in
 *musl ) allow_native=0 ;;
esac

AS_IF([test "$allow_native" -eq 1],
[
    SAVE_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS -march=native"
    AC_MSG_CHECKING([if $CC supports -march=native])
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#include <stdio.h>
                                        int main() { double t; long int q; scanf("%ld",&q);
                                        t = q; printf("t = %f\n",t); return 0; }]]
                                     )
                      ],
                      [ac_cv_native_ok=yes], [ac_cv_native_ok=no])
    AC_MSG_RESULT([$ac_cv_native_ok])
    CFLAGS="$SAVE_CFLAGS"
    AS_IF([test "x$ac_cv_native_ok" = xyes],
          [MORECFLAGS="$MORECFLAGS -march=native"],
          [AC_TRY_CFLAGS([-march=native -mno-avx],
                         [MORECFLAGS="$MORECFLAGS -march=native -mno-avx"])]
         )
])

AC_SUBST(MORECFLAGS)

echo CC="$CC"
echo CFLAGS="$CFLAGS"
echo MORECFLAGS="$MORECFLAGS"

dnl --check if lzcnt instruction is available and requested
AC_MSG_CHECKING(if lzcnt instruction is available and requested)
AS_IF([test "$allow_clz" -eq 1],
      [have_hwlzcnt=0
test -f /proc/cpuinfo && egrep -i "abm|bmi1" /proc/cpuinfo >/dev/null && have_hwlzcnt=1 
test $have_hwlzcnt -eq 0 && sysctl machdep.cpu.features 2>/dev/null | egrep -i "bmi1|abm" >/dev/null && have_hwlzcnt=1
test $have_hwlzcnt -eq 0 && sysctl machdep.cpu.extfeatures 2>/dev/null | egrep -i "lzcnt" >/dev/null && have_hwlzcnt=1
test $have_hwlzcnt -eq 0 && isainfo -x -v 2>/dev/null | egrep -i "bmi1|abm" >/dev/null && have_hwlzcnt=1
test $have_hwlzcnt -eq 0 && test $is_ARM64 -eq 1 && have_hwlzcnt=1 ],
      [have_hwlzcnt=0])
AC_MSG_RESULT($have_hwlzcnt)
AC_SUBST(have_hwlzcnt)

dnl Check for existence of gcc clz extensions
AS_IF([test "$allow_clz" -eq 1],
      [ AC_MSG_CHECKING(if __builtin_clz() is supported)
        AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned int x; x = __builtin_clz(x);}])],
                       [have_clz=1], [have_clz=0])
        AC_MSG_RESULT($have_clz)
        AC_MSG_CHECKING(if __builtin_clzl() is supported)
        AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned long x; x = __builtin_clzl(x);}])],
                       [have_clzl=1], [have_clzl=0])
        AC_MSG_RESULT($have_clzl)
        AC_MSG_CHECKING(if __builtin_clzll() is supported)
        AC_LINK_IFELSE([AC_LANG_SOURCE([int main(){unsigned long long x; x = __builtin_clzll(x);}])],
                       [have_clzll=1], [have_clzll=0])
        AC_MSG_RESULT($have_clzll) ],
      [
        have_clz=0
        have_clzl=0
        have_clzll=0
])
AC_SUBST(have_clz)
AC_SUBST(have_clzl)
AC_SUBST(have_clzll)

AC_SUBST(default_wordsize)

AS_IF([test $have_pipe -eq 1 && test $have_wait -eq 1 && test $have_fork -eq 1],
    [shortg=1
     shortg_or_null=shortg],
    [shortg=0
     shortg_or_null=])
AC_SUBST(shortg_or_null)

AC_CHECK_FUNC(getc_unlocked,have_getc_unlocked=1,have_getc_unlocked=0)
AC_CHECK_FUNC(flockfile,have_flockfile=1,have_flockfile=0)
AC_SUBST(have_flockfile)
AS_IF([test $have_getc_unlocked -eq 1 && test $have_flockfile -eq 1],
      [stdio_nolock=1],[stdio_nolock=0])
AC_SUBST(stdio_nolock)

dnl Set thread-local variables
AX_TLS([tls_supported=1],[tls_supported=0])
AS_IF([test "$tls_supported" -eq 0],[ac_cv_tls=])
AC_SUBST(ac_cv_tls)
AC_SUBST(tls_supported)
AC_SUBST(use_tls)

dnl ANSI terminal control sequences
have_ansicontrols=$allow_ansi
AC_SUBST(have_ansicontrols)

dnl SIGINT catching in dreadnaut
AC_SUBST(allow_interrupt)

dnl Existence and attributes of sort program
AS_IF([test $is_cygwin -eq 1],
      [ AC_PATH_PROGS(sort_prog,[gsort sort],no_sort_found,
             [/usr/bin$PATH_SEPARATOR/bin$PATH_SEPARATOR$PATH])],
      [ AC_CHECK_PROGS(sort_prog,[gsort sort],no_sort_found)])

AC_MSG_CHECKING(if sort supports the -k switch)
AS_IF([{ $sort_prog -k "2,3" </dev/null >/dev/null 2>&1 ; }],
      [sort_newparams_ok=1], [sort_newparams_ok=0])
AC_SUBST(sort_newparams_ok)
AC_MSG_RESULT($sort_newparams_ok)

edit_msg="++++++ This file is automatically generated, don't edit it by hand! ++++++"
AC_SUBST(edit_msg)

AC_PROG_MKDIR_P
AC_PROG_INSTALL

dnl Existence of gunzip program
AC_CHECK_PROGS(gunzip_prog,[gunzip],none)
AS_IF([test "$gunzip_prog" = "gunzip"],[have_gunzip=1],[have_gunzip=0])
AC_SUBST(have_gunzip)

AC_CONFIG_FILES([makefile:makefile.in nauty.h:nauty-h.in
    naututil.h:naututil-h.in gtools.h:gtools-h.in
    nauty.pc:nauty-pc.in])
AC_OUTPUT
