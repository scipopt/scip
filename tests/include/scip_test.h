/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_test.h
 * @brief  SCIP unit test framework header (Unity-based)
 *
 * This header provides Unity-style assertion macros (the TEST_ASSERT* family and
 * the SOFT_ASSERT* extension that records failures and continues execution),
 * together with helpers for writing and running SCIP unit tests.
 */

#ifndef SCIP_TEST_H
#define SCIP_TEST_H

/* Suppress warnings that are expected in test code:
 * - unused-variable/function: not all TUs use every static helper in this header
 * - missing-prototypes: test functions are defined without prior declarations
 * - format-extra-args: the assertion macros always append an empty sentinel
 *   argument (see SCIP_EXPAND below), which no format string consumes. Only
 *   this one format check is disabled; mismatches between conversions and
 *   argument types are still reported.
 * Intentionally no diagnostic push/pop so these apply to the entire test TU.
 */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wmissing-prototypes"
#pragma GCC diagnostic ignored "-Wformat-extra-args"
#endif

/* Define CR_API to disable assertions in SCIP source files that are included
 * directly by tests (e.g., #include "scip/cons_nonlinear.c"). These assertions
 * fail because including .c files results in duplicate static functions at
 * different addresses. See #3543.
 */
#ifndef CR_API
#define CR_API
#endif

/* MSVC only defines M_PI and the other <math.h> constants when _USE_MATH_DEFINES
 * is set before <math.h> is first included (pulled in transitively by scip.h). */
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "scip/scip.h"
#include <string.h>
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#ifdef _WIN32
#include <process.h>
#include <io.h>
#include <stdlib.h>
#define getpid() _getpid()
/* MSVC has no POSIX setenv; _putenv_s always overwrites, so the overwrite flag
 * (last argument) is dropped. Assigning an empty value is how _putenv_s removes
 * a variable, which gives us unsetenv as well. */
#define setenv(name, value, overwrite) _putenv_s((name), (value))
#define unsetenv(name) _putenv_s((name), "")
/* MSVC spells the POSIX fd helpers with a leading underscore */
#define scip_test_dup    _dup
#define scip_test_dup2   _dup2
#define scip_test_fileno _fileno
#define scip_test_close  _close
#else
#include <unistd.h>
#define scip_test_dup    dup
#define scip_test_dup2   dup2
#define scip_test_fileno fileno
#define scip_test_close  close
#endif

/* Directory prefix for capture temp files: /tmp on POSIX (avoids littering the
 * working tree), the current directory on Windows (which has no /tmp). */
#ifdef _WIN32
#define SCIP_TEST_TMP_PREFIX ""
#else
#define SCIP_TEST_TMP_PREFIX "/tmp/"
#endif

/* the method TESTsetSCIPStage(scip, stage) can be called in SCIP_STAGE_PROBLEM and can get to
 *  SCIP_STAGE_TRANSFORMED
 *  SCIP_STAGE_PRESOLVING
 *  SCIP_STAGE_PRESOLVED
 *  SCIP_STAGE_SOLVING
 *  SCIP_STAGE_SOLVED
 *
 *  If stage == SCIP_STAGE_SOLVING and enableNLP is true, then SCIP will build its NLP
 */
SCIP_RETCODE TESTscipSetStage(SCIP* scip, SCIP_STAGE stage, SCIP_Bool enableNLP);

/** assembles path of testfile using another files directory as directory name */
void TESTsetTestfilename(
   char*                 filename,           /**< buffer to write to, assumed to have length at least SCIP_MAXSTRLEN */
   const char*           file,               /**< name of file, usually including full path, from which to take directory name */
   const char*           testfile            /**< name of file to append, assumed to be in same directory as file */
);

/* Include the .c file here because the plugins implemented in scip_test.c should
 * be available every test.
 * */
#include "scip_test.c"

/*
 * Unity Test Framework
 */
#define UNITY_INCLUDE_DOUBLE
#include "unity/unity.h"

/*
 * Unity Assertion Wrappers (C99-compatible)
 *
 * These wrappers accept Criterion-style messages, including printf-style ones
 * such as ("expected %d, got %d", a, b), and forward them to Unity.
 *
 * Two-layer macro pattern: the outer macro takes (...) and appends a sentinel
 * so the inner macro's "..." always receives at least one argument. This avoids
 * the C99 pedantic warning about zero variadic arguments. The sentinel also
 * doubles as the "no message given" marker (see scip_test_msg below).
 *
 * SCIP_EXPAND works around the MSVC traditional preprocessor bug where
 * __VA_ARGS__ is not properly split when passed to another macro.
 */
#define SCIP_EXPAND(...) __VA_ARGS__

/* Criterion's assertion messages are printf-style, but Unity's message
 * parameter is a plain string, so the message has to be formatted before it is
 * handed over. Returns NULL for the empty sentinel, which makes Unity fall back
 * to its own default description.
 *
 * A single shared buffer is enough: Unity consumes the string before the next
 * assertion can run, and the unit tests are single-threaded.
 */
static char scip_test_msgbuf[1024];

/* Sentinel appended by the assertion wrappers when no user message is given.
 * It must be a non-empty string literal: GCC's -Wformat-zero-length rejects an
 * empty one, and using a literal (rather than NULL) keeps -Wformat-security
 * quiet on distros where -Wformat implies -Wformat=2. Recognized by content in
 * the helpers below, which then fall back to the default wording.
 */
#define SCIP_TEST_NOMSG " "

#ifdef __GNUC__
__attribute__((format(printf, 1, 2)))
#endif
static const char* scip_test_msg(const char* fmt, ...)
{
   va_list ap;

   if( fmt == NULL || fmt[0] == '\0' || strcmp(fmt, SCIP_TEST_NOMSG) == 0 )
      return NULL;

   va_start(ap, fmt);
   (void) vsnprintf(scip_test_msgbuf, sizeof(scip_test_msgbuf), fmt, ap);
   va_end(ap);

   return scip_test_msgbuf;
}

/* Like scip_test_msg, but when no user message was given (the empty sentinel)
 * returns a caller-supplied default instead of NULL. A failing assertion with no
 * explicit message would otherwise print a bare "FAIL" with no detail at all, so
 * the boolean and comparison macros pass the same defaults Unity itself uses.
 */
#ifdef __GNUC__
__attribute__((format(printf, 2, 3)))
#endif
static const char* scip_test_msg_default(const char* def, const char* fmt, ...)
{
   va_list ap;

   if( fmt == NULL || fmt[0] == '\0' || strcmp(fmt, SCIP_TEST_NOMSG) == 0 )
      return def;

   va_start(ap, fmt);
   (void) vsnprintf(scip_test_msgbuf, sizeof(scip_test_msgbuf), fmt, ap);
   va_end(ap);

   return scip_test_msgbuf;
}

#undef TEST_ASSERT_TRUE
#define TEST_ASSERT_TRUE(...) SCIP_EXPAND(TEST_ASSERT_TRUE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_TRUE_(cond, ...) UNITY_TEST_ASSERT((cond), __LINE__, scip_test_msg_default(" Expected TRUE Was FALSE", __VA_ARGS__))

#undef TEST_ASSERT_FALSE
#define TEST_ASSERT_FALSE(...) SCIP_EXPAND(TEST_ASSERT_FALSE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_FALSE_(cond, ...) UNITY_TEST_ASSERT(!(cond), __LINE__, scip_test_msg_default(" Expected FALSE Was TRUE", __VA_ARGS__))

/* TEST_ASSERT and TEST_ASSERT_NOT are the canonical condition checks and mirror
 * Criterion's cr_assert and cr_assert_not. Unity also spells the same two checks
 * TEST_ASSERT_TRUE/TEST_ASSERT_FALSE and TEST_ASSERT/TEST_ASSERT_UNLESS; keep all
 * of those as working spellings and let every one accept the optional message.
 */
#undef TEST_ASSERT
#define TEST_ASSERT(...) SCIP_EXPAND(TEST_ASSERT_TRUE_(__VA_ARGS__, SCIP_TEST_NOMSG))

#undef TEST_ASSERT_NOT
#define TEST_ASSERT_NOT(...) SCIP_EXPAND(TEST_ASSERT_FALSE_(__VA_ARGS__, SCIP_TEST_NOMSG))

#undef TEST_ASSERT_UNLESS
#define TEST_ASSERT_UNLESS(...) SCIP_EXPAND(TEST_ASSERT_FALSE_(__VA_ARGS__, SCIP_TEST_NOMSG))

#undef TEST_ASSERT_EQUAL
#define TEST_ASSERT_EQUAL(...) SCIP_EXPAND(TEST_ASSERT_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
/* Unlike UNITY_TEST_ASSERT, the UnityAssert* functions take the message as an
 * ordinary argument, so it would be formatted on every *passing* assertion too.
 * The macros below therefore pre-check the condition and only format when it
 * looks like a failure. Unity still performs the real comparison and decides
 * pass/fail; a pre-check that disagrees can at worst drop or waste a message,
 * never change the verdict. Operands go through locals so that they are
 * evaluated exactly once, as before.
 */
#define TEST_ASSERT_EQUAL_(expected, actual, ...) do { \
    UNITY_INT scip_test_e_ = (UNITY_INT)(expected); \
    UNITY_INT scip_test_a_ = (UNITY_INT)(actual); \
    UNITY_TEST_ASSERT_EQUAL_INT(scip_test_e_, scip_test_a_, __LINE__, \
        scip_test_e_ == scip_test_a_ ? NULL : scip_test_msg(__VA_ARGS__)); \
} while(0)

#undef TEST_ASSERT_NOT_EQUAL
#define TEST_ASSERT_NOT_EQUAL(...) SCIP_EXPAND(TEST_ASSERT_NOT_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_NOT_EQUAL_(expected, actual, ...) UNITY_TEST_ASSERT(((expected) != (actual)), __LINE__, scip_test_msg_default(" Expected Not-Equal", __VA_ARGS__))

#undef TEST_ASSERT_NULL
#define TEST_ASSERT_NULL(...) SCIP_EXPAND(TEST_ASSERT_NULL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_NULL_(ptr, ...) UNITY_TEST_ASSERT_NULL((ptr), __LINE__, scip_test_msg(__VA_ARGS__))

#undef TEST_ASSERT_NOT_NULL
#define TEST_ASSERT_NOT_NULL(...) SCIP_EXPAND(TEST_ASSERT_NOT_NULL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_NOT_NULL_(ptr, ...) UNITY_TEST_ASSERT_NOT_NULL((ptr), __LINE__, scip_test_msg(__VA_ARGS__))

#undef TEST_ASSERT_EQUAL_STRING
#define TEST_ASSERT_EQUAL_STRING(...) SCIP_EXPAND(TEST_ASSERT_EQUAL_STRING_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_EQUAL_STRING_(expected, actual, ...) do { \
    const char* scip_test_e_ = (const char*)(expected); \
    const char* scip_test_a_ = (const char*)(actual); \
    UNITY_TEST_ASSERT_EQUAL_STRING(scip_test_e_, scip_test_a_, __LINE__, \
        (scip_test_e_ != NULL && scip_test_a_ != NULL && strcmp(scip_test_e_, scip_test_a_) == 0) \
        ? NULL : scip_test_msg(__VA_ARGS__)); \
} while(0)

#undef TEST_ASSERT_EQUAL_MEMORY
#define TEST_ASSERT_EQUAL_MEMORY(...) SCIP_EXPAND(TEST_ASSERT_EQUAL_MEMORY_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_EQUAL_MEMORY_(expected, actual, size, ...) do { \
    const void* scip_test_e_ = (const void*)(expected); \
    const void* scip_test_a_ = (const void*)(actual); \
    size_t scip_test_n_ = (size_t)(size); \
    UNITY_TEST_ASSERT_EQUAL_MEMORY(scip_test_e_, scip_test_a_, scip_test_n_, __LINE__, \
        (scip_test_e_ != NULL && scip_test_a_ != NULL && memcmp(scip_test_e_, scip_test_a_, scip_test_n_) == 0) \
        ? NULL : scip_test_msg(__VA_ARGS__)); \
} while(0)

#undef TEST_ASSERT_DOUBLE_WITHIN
#define TEST_ASSERT_DOUBLE_WITHIN(...) SCIP_EXPAND(TEST_ASSERT_DOUBLE_WITHIN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_DOUBLE_WITHIN_(actual, expected, delta, ...) do { \
    UNITY_DOUBLE scip_test_a_ = (UNITY_DOUBLE)(actual); \
    UNITY_DOUBLE scip_test_e_ = (UNITY_DOUBLE)(expected); \
    UNITY_DOUBLE scip_test_d_ = (UNITY_DOUBLE)(delta); \
    UNITY_TEST_ASSERT_DOUBLE_WITHIN(scip_test_d_, scip_test_e_, scip_test_a_, __LINE__, \
        fabs(scip_test_a_ - scip_test_e_) <= scip_test_d_ ? NULL : scip_test_msg(__VA_ARGS__)); \
} while(0)

/* NOTE: These match Criterion semantics: first arg comparison second arg
 * e.g., TEST_ASSERT_LESS_THAN(a, b) asserts a < b
 */
#undef TEST_ASSERT_LESS_THAN
#define TEST_ASSERT_LESS_THAN(...) SCIP_EXPAND(TEST_ASSERT_LESS_THAN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_LESS_THAN_(a, b, ...) UNITY_TEST_ASSERT(((a) < (b)), __LINE__, scip_test_msg_default(" Expected LESS THAN", __VA_ARGS__))

#undef TEST_ASSERT_GREATER_THAN
#define TEST_ASSERT_GREATER_THAN(...) SCIP_EXPAND(TEST_ASSERT_GREATER_THAN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_GREATER_THAN_(a, b, ...) UNITY_TEST_ASSERT(((a) > (b)), __LINE__, scip_test_msg_default(" Expected GREATER THAN", __VA_ARGS__))

#undef TEST_ASSERT_LESS_OR_EQUAL
#define TEST_ASSERT_LESS_OR_EQUAL(...) SCIP_EXPAND(TEST_ASSERT_LESS_OR_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_LESS_OR_EQUAL_(a, b, ...) UNITY_TEST_ASSERT(((a) <= (b)), __LINE__, scip_test_msg_default(" Expected LESS_OR_EQUAL", __VA_ARGS__))

#undef TEST_ASSERT_GREATER_OR_EQUAL
#define TEST_ASSERT_GREATER_OR_EQUAL(...) SCIP_EXPAND(TEST_ASSERT_GREATER_OR_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_GREATER_OR_EQUAL_(a, b, ...) UNITY_TEST_ASSERT(((a) >= (b)), __LINE__, scip_test_msg_default(" Expected GREATER_OR_EQUAL", __VA_ARGS__))

#undef TEST_ASSERT_EQUAL_DOUBLE
#define TEST_ASSERT_EQUAL_DOUBLE(...) SCIP_EXPAND(TEST_ASSERT_EQUAL_DOUBLE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_EQUAL_DOUBLE_(expected, actual, ...) do { \
    UNITY_DOUBLE scip_test_e_ = (UNITY_DOUBLE)(expected); \
    UNITY_DOUBLE scip_test_a_ = (UNITY_DOUBLE)(actual); \
    UNITY_TEST_ASSERT_EQUAL_DOUBLE(scip_test_e_, scip_test_a_, __LINE__, \
        scip_test_e_ == scip_test_a_ ? NULL : scip_test_msg(__VA_ARGS__)); \
} while(0)

/*
 * Soft Assertion Support
 *
 * Criterion's cr_expect macros continue execution even on failure and report
 * all failures at the end. We emulate this with manual tracking.
 */
static int scip_test_soft_failure_count = 0;
static char scip_test_soft_failure_messages[8192];

#define SOFT_ASSERT_RESET() do { \
    scip_test_soft_failure_count = 0; \
    scip_test_soft_failure_messages[0] = '\0'; \
} while(0)

#define SOFT_ASSERT_CHECK() do { \
    if( scip_test_soft_failure_count > 0 ) { \
        TEST_FAIL_MESSAGE(scip_test_soft_failure_messages); \
    } \
} while(0)

/*
 * Signal Handling Support
 *
 * For tests that expect signals (e.g., SIGABRT).
 */
static jmp_buf scip_test_signal_jmp;
static volatile sig_atomic_t scip_test_signal_caught = 0;

static void scip_test_signal_handler(int sig)
{
   scip_test_signal_caught = sig;
   longjmp(scip_test_signal_jmp, 1);
}

#define TEST_EXPECT_SIGNAL(sig, code) do { \
    signal(sig, scip_test_signal_handler); \
    scip_test_signal_caught = 0; \
    if( setjmp(scip_test_signal_jmp) == 0 ) { \
        code; \
        TEST_FAIL_MESSAGE("Expected signal " #sig " but none was raised"); \
    } else { \
        TEST_ASSERT_EQUAL_INT(sig, (int)scip_test_signal_caught); \
    } \
    signal(sig, SIG_DFL); \
} while(0)

/*
 * Soft Assertion Support (SOFT_ASSERT_*)
 *
 * Soft assertions continue execution even on failure and report all failures
 * at the end of the test. Unity doesn't have this built-in.
 *
 * Usage: Call SOFT_ASSERT_RESET() in setUp(), SOFT_ASSERT_CHECK() in tearDown().
 */
/* Record one soft-assert failure. "detail" describes the failed condition and is
 * assembled from stringified macro arguments; "msg" is the optional user message
 * (NULL when none was given), appended so that soft assertions report the same
 * information as the hard ones.
 */
static void scip_test_soft_fail(
   const char*           detail,             /**< description of the failed condition */
   const char*           file,               /**< file the assertion appears in */
   int                   line,               /**< line the assertion appears in */
   const char*           msg                 /**< user message, or NULL */
   )
{
   size_t len = strlen(scip_test_soft_failure_messages);

   scip_test_soft_failure_count++;
   (void) snprintf(scip_test_soft_failure_messages + len, sizeof(scip_test_soft_failure_messages) - len,
      "  %s at %s:%d%s%s\n", detail, file, line, msg != NULL ? " - " : "", msg != NULL ? msg : "");
}

#define SOFT_ASSERT(...) SCIP_EXPAND(SOFT_ASSERT_TRUE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_TRUE(...) SCIP_EXPAND(SOFT_ASSERT_TRUE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_TRUE_(cond, ...) do { \
    if( !(cond) ) \
        scip_test_soft_fail("SOFT_ASSERT FAILED: " #cond, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_NOT(...) SCIP_EXPAND(SOFT_ASSERT_FALSE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_FALSE(...) SCIP_EXPAND(SOFT_ASSERT_FALSE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_FALSE_(cond, ...) do { \
    if( (cond) ) \
        scip_test_soft_fail("SOFT_ASSERT_NOT FAILED: " #cond, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_EQUAL(...) SCIP_EXPAND(SOFT_ASSERT_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_EQUAL_(expected, actual, ...) do { \
    if( (actual) != (expected) ) \
        scip_test_soft_fail("SOFT_ASSERT_EQUAL FAILED: " #actual " != " #expected, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_NOT_EQUAL(...) SCIP_EXPAND(SOFT_ASSERT_NOT_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_NOT_EQUAL_(expected, actual, ...) do { \
    if( (actual) == (expected) ) \
        scip_test_soft_fail("SOFT_ASSERT_NOT_EQUAL FAILED: " #actual " == " #expected, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

/* NOTE: These match Criterion semantics: SOFT_ASSERT_LESS_THAN(a, b) asserts a < b */
#define SOFT_ASSERT_LESS_THAN(...) SCIP_EXPAND(SOFT_ASSERT_LESS_THAN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_LESS_THAN_(a, b, ...) do { \
    if( !((a) < (b)) ) \
        scip_test_soft_fail("SOFT_ASSERT_LESS_THAN FAILED: " #a " >= " #b, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_GREATER_THAN(...) SCIP_EXPAND(SOFT_ASSERT_GREATER_THAN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_GREATER_THAN_(a, b, ...) do { \
    if( !((a) > (b)) ) \
        scip_test_soft_fail("SOFT_ASSERT_GREATER_THAN FAILED: " #a " <= " #b, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_LESS_OR_EQUAL(...) SCIP_EXPAND(SOFT_ASSERT_LESS_OR_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_LESS_OR_EQUAL_(a, b, ...) do { \
    if( !((a) <= (b)) ) \
        scip_test_soft_fail("SOFT_ASSERT_LESS_OR_EQUAL FAILED: " #a " > " #b, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_GREATER_OR_EQUAL(...) SCIP_EXPAND(SOFT_ASSERT_GREATER_OR_EQUAL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_GREATER_OR_EQUAL_(a, b, ...) do { \
    if( !((a) >= (b)) ) \
        scip_test_soft_fail("SOFT_ASSERT_GREATER_OR_EQUAL FAILED: " #a " < " #b, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_NULL(...) SCIP_EXPAND(SOFT_ASSERT_NULL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_NULL_(ptr, ...) do { \
    if( (ptr) != NULL ) \
        scip_test_soft_fail("SOFT_ASSERT_NULL FAILED: " #ptr, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_NOT_NULL(...) SCIP_EXPAND(SOFT_ASSERT_NOT_NULL_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_NOT_NULL_(ptr, ...) do { \
    if( (ptr) == NULL ) \
        scip_test_soft_fail("SOFT_ASSERT_NOT_NULL FAILED: " #ptr, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_DOUBLE_WITHIN(...) SCIP_EXPAND(SOFT_ASSERT_DOUBLE_WITHIN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_DOUBLE_WITHIN_(actual, expected, delta, ...) do { \
    double _a = (actual), _e = (expected), _d = (delta); \
    if( (_a - _e) > _d || (_e - _a) > _d ) \
        scip_test_soft_fail("SOFT_ASSERT_DOUBLE_WITHIN FAILED: " #actual " != " #expected " (delta=" #delta ")", __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

/* EXPECTFEQ - soft floating-point equality check within a tolerance. Reports a
 * failure (at the end of the test) if |a - b| exceeds the tolerance. The default
 * is 1e-6; EXPECTFEQ_WITHIN lets callers pick a different one (the interval
 * arithmetic tests use 1e-12). Centralized here so every test area shares the
 * same definition, including the MSVC-safe ABS((a) - (b)) expansion: a bare
 * a-b with a negative b would otherwise become '--' under MSVC's preprocessor.
 */
#define EXPECTFEQ_WITHIN(tol, a, b) SOFT_ASSERT_DOUBLE_WITHIN((a), (b), (tol), "%s = %g != %g (dif %g)", #a, (a), (b), ABS((a) - (b)))
#define EXPECTFEQ(a, b) EXPECTFEQ_WITHIN(1e-6, (a), (b))

#define SOFT_ASSERT_EQUAL_STRING(...) SCIP_EXPAND(SOFT_ASSERT_EQUAL_STRING_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_EQUAL_STRING_(expected, actual, ...) do { \
    if( strcmp((actual), (expected)) != 0 ) \
        scip_test_soft_fail("SOFT_ASSERT_EQUAL_STRING FAILED: " #actual " != " #expected, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

/* Soft assert that a and b are NOT within delta of each other (for floats/doubles) */
/* NOTE: as for TEST_ASSERT_EQUAL_MEMORY, size is a number of bytes */
#define SOFT_ASSERT_EQUAL_MEMORY(...) SCIP_EXPAND(SOFT_ASSERT_EQUAL_MEMORY_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_EQUAL_MEMORY_(expected, actual, size, ...) do { \
    if( memcmp((actual), (expected), (size)) != 0 ) \
        scip_test_soft_fail("SOFT_ASSERT_EQUAL_MEMORY FAILED: " #actual " != " #expected, __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

#define SOFT_ASSERT_DOUBLE_NOT_WITHIN(...) SCIP_EXPAND(SOFT_ASSERT_DOUBLE_NOT_WITHIN_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define SOFT_ASSERT_DOUBLE_NOT_WITHIN_(actual, expected, delta, ...) do { \
    double _a = (actual), _e = (expected), _d = (delta); \
    if( !((_a - _e) > _d || (_e - _a) > _d) ) \
        scip_test_soft_fail("SOFT_ASSERT_DOUBLE_NOT_WITHIN FAILED: " #actual " == " #expected " (delta=" #delta ")", __FILE__, __LINE__, scip_test_msg(__VA_ARGS__)); \
} while(0)

/* TEST_LOG_INFO - informational logging during a test. Criterion writes this to its
 * own log stream rather than the test's stdout, so route it to stderr to avoid
 * polluting TEST_CAPTURE_STDOUT() comparisons. */
#define TEST_LOG_INFO(...) (void)fprintf(stderr, __VA_ARGS__)

/* Float eq with infinity handling - compares floats, treating infinities specially */
#define TEST_ASSERT_DOUBLE_WITHIN_INF(...) SCIP_EXPAND(TEST_ASSERT_DOUBLE_WITHIN_INF_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_DOUBLE_WITHIN_INF_(actual, expected, eps, ...) do { \
    double _a = (actual), _e = (expected), _eps = (eps); \
    if( _a == _e ) { /* handles infinities */ } \
    else if( (_a - _e) <= _eps && (_e - _a) <= _eps ) { } \
    else { TEST_FAIL_MESSAGE("Float assertion failed"); } \
} while(0)

/*
 * Stdout/Stderr Capture Support
 *
 * For tests that verify output. Using temporary files since fmemopen
 * may not be available on all platforms.
 */
/* Saved file descriptors of the real stdout/stderr while redirected. We capture
 * via freopen()+dup() rather than `stdout = fopen(...)` because stdout/stderr are
 * not assignable lvalues on MSVC. -1 means "not currently redirected". */
static int scip_test_stdout_savedfd = -1;
static int scip_test_stderr_savedfd = -1;
static char scip_test_stdout_buffer[65536];
static char scip_test_stderr_buffer[65536];
static char scip_test_stdout_tmpfile[256] = "";
static char scip_test_stderr_tmpfile[256] = "";

/* Helper to read captured output into buffer */
static void scip_test_read_captured_file(const char* filename, char* buffer, size_t bufsize)
{
   FILE* f = fopen(filename, "r");
   if( f != NULL )
   {
      size_t n = fread(buffer, 1, bufsize - 1, f);
      buffer[n] = '\0';
      fclose(f);
   }
   else
   {
      buffer[0] = '\0';
   }
}

/* Restore a redirected standard stream (stdout/stderr) from its saved fd.
 * Uses dup2() rather than `stream = saved` because stdout/stderr are not
 * assignable lvalues on MSVC. After this the captured temp file is closed
 * (its fd was replaced), so it can be reopened for reading. */
static void scip_test_restore_stream(FILE* stream, int* savedfd)
{
   if( *savedfd != -1 )
   {
      fflush(stream);
      scip_test_dup2(*savedfd, scip_test_fileno(stream));
      scip_test_close(*savedfd);
      *savedfd = -1;
   }
}

/* Redirect a standard stream to a temp file so that a test can inspect what was
 * written to it. The compare macros below restore the stream again.
 *
 * Capturing a stream that is still captured (a test that redirects but never
 * compares, of which there are several) would leak the saved descriptor and
 * leave the stream pointing at the temp file for the rest of the process, so
 * restore an active capture first. Criterion forked per test and did not need
 * this; all tests of a Unity binary share one process.
 */
static void scip_test_capture_stream(
   FILE*                 stream,             /**< stream to redirect */
   int*                  savedfd,            /**< buffer to store the stream's original descriptor */
   char*                 tmpfile,            /**< buffer to store the name of the temp file */
   size_t                tmpfilesize,        /**< size of the temp file name buffer */
   const char*           name                /**< stream name, used to build the temp file name */
   )
{
   scip_test_restore_stream(stream, savedfd);

   fflush(stream);
   *savedfd = scip_test_dup(scip_test_fileno(stream));
   (void) snprintf(tmpfile, tmpfilesize, SCIP_TEST_TMP_PREFIX "scip_test_%s_%d", name, (int)getpid());
   if( freopen(tmpfile, "w+", stream) == NULL ) {}
}

#define TEST_CAPTURE_STDOUT() scip_test_capture_stream(stdout, &scip_test_stdout_savedfd, \
    scip_test_stdout_tmpfile, sizeof(scip_test_stdout_tmpfile), "stdout")

#define TEST_CAPTURE_STDERR() scip_test_capture_stream(stderr, &scip_test_stderr_savedfd, \
    scip_test_stderr_tmpfile, sizeof(scip_test_stderr_tmpfile), "stderr")

/* The SOFT_ASSERT_* macros only collect failures; something has to inspect the
 * result afterwards or they assert nothing at all. Leaving that to each test
 * file is easy to forget, so hook it into RUN_TEST: reset before the test and
 * report after it. The report runs inside Unity's protected block, so
 * TEST_FAIL_MESSAGE behaves as it does in an ordinary assertion.
 *
 * Redirected streams are restored before the report and again before the next
 * test: Unity writes its results to stdout, so a test that captures a stream
 * and does not compare it would otherwise swallow its own failure message.
 */
static void (*scip_test_body)(void);

/* Restore both standard streams at a test boundary and drop the capture temp
 * file. A test that compares a captured stream removes the file itself, but a
 * test that only captures to keep expected error output off the console (as the
 * datatree tests do) leaves it behind, once per process and never cleaned up.
 */
static void scip_test_restore_streams(void)
{
   scip_test_restore_stream(stdout, &scip_test_stdout_savedfd);
   scip_test_restore_stream(stderr, &scip_test_stderr_savedfd);

   if( scip_test_stdout_tmpfile[0] != '\0' )
   {
      (void) remove(scip_test_stdout_tmpfile);
      scip_test_stdout_tmpfile[0] = '\0';
   }

   if( scip_test_stderr_tmpfile[0] != '\0' )
   {
      (void) remove(scip_test_stderr_tmpfile);
      scip_test_stderr_tmpfile[0] = '\0';
   }
}

static void scip_test_run_body(void)
{
   scip_test_body();
   scip_test_restore_streams();
   SOFT_ASSERT_CHECK();
}

#undef RUN_TEST
#define RUN_TEST(func) do { \
    scip_test_body = (func); \
    scip_test_restore_streams(); \
    SOFT_ASSERT_RESET(); \
    UnityDefaultTestRun(scip_test_run_body, #func, __LINE__); \
} while(0)

#define TEST_ASSERT_STDOUT_EQUAL_STRING(...) SCIP_EXPAND(TEST_ASSERT_STDOUT_EQUAL_STRING_(__VA_ARGS__, "", ""))
#define TEST_ASSERT_STDOUT_EQUAL_STRING_(expected, msg, ...) do { \
    if( scip_test_stdout_savedfd != -1 ) { \
        scip_test_restore_stream(stdout, &scip_test_stdout_savedfd); \
        scip_test_read_captured_file(scip_test_stdout_tmpfile, scip_test_stdout_buffer, sizeof(scip_test_stdout_buffer)); \
        remove(scip_test_stdout_tmpfile); \
    } \
    TEST_ASSERT_EQUAL_STRING_MESSAGE(expected, scip_test_stdout_buffer, msg); \
} while(0)

#define TEST_ASSERT_STDERR_EQUAL_STRING(...) SCIP_EXPAND(TEST_ASSERT_STDERR_EQUAL_STRING_(__VA_ARGS__, "", ""))
#define TEST_ASSERT_STDERR_EQUAL_STRING_(expected, msg, ...) do { \
    if( scip_test_stderr_savedfd != -1 ) { \
        scip_test_restore_stream(stderr, &scip_test_stderr_savedfd); \
        scip_test_read_captured_file(scip_test_stderr_tmpfile, scip_test_stderr_buffer, sizeof(scip_test_stderr_buffer)); \
        remove(scip_test_stderr_tmpfile); \
    } \
    TEST_ASSERT_EQUAL_STRING_MESSAGE(expected, scip_test_stderr_buffer, msg); \
} while(0)

/*
 * File Comparison Support
 *
 * For tests that compare file contents (e.g., LP file round-trip tests).
 */

/** Compare two FILE* streams line by line (hard assertion) */
#define TEST_ASSERT_FILES_EQUAL(file1, file2) do { \
    char _buf1[4096], _buf2[4096]; \
    int _line = 0; \
    rewind(file1); \
    rewind(file2); \
    while( 1 ) { \
        char* _r1 = fgets(_buf1, sizeof(_buf1), (file1)); \
        char* _r2 = fgets(_buf2, sizeof(_buf2), (file2)); \
        _line++; \
        if( _r1 == NULL && _r2 == NULL ) break; \
        if( _r1 == NULL || _r2 == NULL ) { \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), "Files differ at line %d: one file ended early", _line); \
            TEST_FAIL_MESSAGE(_msg); \
        } \
        if( strcmp(_buf1, _buf2) != 0 ) { \
            char _msg[256]; \
            snprintf(_msg, sizeof(_msg), "Files differ at line %d", _line); \
            TEST_FAIL_MESSAGE(_msg); \
        } \
    } \
} while(0)

/** Compare two FILE* streams line by line (soft assertion) */
#define SOFT_ASSERT_FILES_EQUAL(file1, file2) do { \
    char _buf1[4096], _buf2[4096]; \
    int _line = 0; \
    rewind(file1); \
    rewind(file2); \
    while( 1 ) { \
        char* _r1 = fgets(_buf1, sizeof(_buf1), (file1)); \
        char* _r2 = fgets(_buf2, sizeof(_buf2), (file2)); \
        _line++; \
        if( _r1 == NULL && _r2 == NULL ) break; \
        if( _r1 == NULL || _r2 == NULL ) { \
            scip_test_soft_failure_count++; \
            snprintf(scip_test_soft_failure_messages + strlen(scip_test_soft_failure_messages), \
                     sizeof(scip_test_soft_failure_messages) - strlen(scip_test_soft_failure_messages), \
                     "  Files differ at line %d: one file ended early at %s:%d\n", _line, __FILE__, __LINE__); \
            break; \
        } \
        if( strcmp(_buf1, _buf2) != 0 ) { \
            scip_test_soft_failure_count++; \
            snprintf(scip_test_soft_failure_messages + strlen(scip_test_soft_failure_messages), \
                     sizeof(scip_test_soft_failure_messages) - strlen(scip_test_soft_failure_messages), \
                     "  Files differ at line %d at %s:%d\n", _line, __FILE__, __LINE__); \
        } \
    } \
} while(0)

/** Compare captured stdout with a reference FILE* (hard assertion).
 *  Closes captured stdout, reads temp file, compares with reffile. */
#define TEST_ASSERT_STDOUT_EQUAL_FILE(...) SCIP_EXPAND(TEST_ASSERT_STDOUT_EQUAL_FILE_(__VA_ARGS__, SCIP_TEST_NOMSG))
#define TEST_ASSERT_STDOUT_EQUAL_FILE_(reffile, ...) do { \
    if( scip_test_stdout_savedfd != -1 ) { \
        FILE* _captured; \
        scip_test_restore_stream(stdout, &scip_test_stdout_savedfd); \
        _captured = fopen(scip_test_stdout_tmpfile, "r"); \
        if( _captured == NULL ) { \
            remove(scip_test_stdout_tmpfile); \
            TEST_FAIL_MESSAGE("Failed to open captured stdout file"); \
        } \
        TEST_ASSERT_FILES_EQUAL(_captured, reffile); \
        fclose(_captured); \
        remove(scip_test_stdout_tmpfile); \
    } else { \
        TEST_FAIL_MESSAGE("stdout was not redirected"); \
    } \
} while(0)

/*
 * SCIP_CALL override for test framework
 *
 * Use Unity assertions for error checking.
 */
#undef SCIP_CALL
#define SCIP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             char _msg_[256];                                                                 \
                             snprintf(_msg_, sizeof(_msg_), "Error <%d> in function call", _restat_);         \
                             TEST_FAIL_MESSAGE(_msg_);                                                        \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/*
 * Suite-level setup/teardown support
 *
 * In Criterion, setup/teardown ran once per TestSuite. In Unity, setUp/tearDown
 * run before/after every test. This causes significant slowdown when setup
 * creates expensive resources like SCIP instances.
 *
 * Use SCIP_SUITE_SETUP/TEARDOWN in your setUp/tearDown functions for
 * suite-level initialization that only runs once:
 *
 *   void setUp(void) { SCIP_SUITE_SETUP(setup); }
 *   void tearDown(void) { SCIP_SUITE_TEARDOWN(teardown); }
 *
 * The setup function runs only on first call; teardown runs via atexit().
 */
static int scip_suite_initialized = 0;
static void (*scip_suite_teardown_fn)(void) = NULL;

static void scip_suite_atexit_handler(void)
{
   if( scip_suite_teardown_fn != NULL )
   {
      scip_suite_teardown_fn();
      scip_suite_teardown_fn = NULL;
   }
}

#define SCIP_SUITE_SETUP(setup_fn) do { \
   if( !scip_suite_initialized ) { \
      setup_fn(); \
      scip_suite_initialized = 1; \
   } \
} while(0)

#define SCIP_SUITE_TEARDOWN(teardown_fn) do { \
   if( scip_suite_initialized && scip_suite_teardown_fn == NULL ) { \
      scip_suite_teardown_fn = teardown_fn; \
      atexit(scip_suite_atexit_handler); \
   } \
} while(0)

#endif /* SCIP_TEST_H */
