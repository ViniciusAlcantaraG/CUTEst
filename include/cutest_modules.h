/* \file cutest_modules.h */

/*
 * assign names for each CUTEst package using the C pre-processor.
 * possibilities are long (i8) and normal (i4, default) integers and
 * half (r2), single (r4), double (r8, default) and quadruple (r16) reals
 *
 * Nick Gould for CUTEst
 * initial version, 2023-11-02
 * this version 2023-11-02
 */

#ifdef CUTEST_LONG
#define CUTEST_KINDS_integer CUTEST_KINDS_long
#else
#define CUTEST_KINDS_integer CUTEST_KINDS_int
#endif

#if defined CUTEST_HALF
#ifdef CUTEST_LONG
#define CUTEST_KINDS_precision CUTEST_KINDS_half_long
#define CUTEST_precision CUTEST_half_long
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_half_long
#define CUTEST_iNTERFACE_precision CUTEST_INTERFACE_half_long
#define CUTEST_LQP_precision CUTEST_LQP_half_long
#else
#define CUTEST_KINDS_precision CUTEST_KINDS_half
#define CUTEST_precision CUTEST_half
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_half
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_half
#define CUTEST_LQP_precision CUTEST_LQP_half
#endif
#elif defined CUTEST_SINGLE
#ifdef CUTEST_LONG
#define CUTEST_KINDS_precision CUTEST_KINDS_single_long
#define CUTEST_precision CUTEST_single_long
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_single_long
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_single_long
#define CUTEST_LQP_precision CUTEST_LQP_single_long
#else
#define CUTEST_KINDS_precision CUTEST_KINDS_single
#define CUTEST_precision CUTEST_single
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_single
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_single
#define CUTEST_LQP_precision CUTEST_LQP_single
#endif
#elif defined CUTEST_QUAD
#ifdef CUTEST_LONG
#define CUTEST_KINDS_precision CUTEST_KINDS_quadruple_long
#define CUTEST_precision CUTEST_quadruple_long
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_quadruple_long
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_quadruple_long
#define CUTEST_LQP_precision CUTEST_LQP_quadruple_long
#else
#define CUTEST_KINDS_precision CUTEST_KINDS_quadruple
#define CUTEST_precision CUTEST_quadruple
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_quadruple
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_quadruple
#define CUTEST_LQP_precision CUTEST_LQP_quadruple
#endif
#else
#ifdef CUTEST_LONG
#define CUTEST_KINDS_precision CUTEST_KINDS_double_long
#define CUTEST_precision CUTEST_double_long
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_double_long
#define CUTEST_iNTERFACE_precision CUTEST_INTERFACE_double_long
#define CUTEST_LQP_precision CUTEST_LQP_double_long
#else
#define CUTEST_KINDS_precision CUTEST_KINDS_double
#define CUTEST_precision CUTEST_double
#define CUTEST_PROBLEM_precision CUTEST_PROBLEM_double
#define CUTEST_INTERFACE_precision CUTEST_INTERFACE_double
#define CUTEST_LQP_precision CUTEST_LQP_double
#endif
#endif
