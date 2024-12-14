!!#if   defined(__MIC__) || defined(__AVX512F__)
!!#define SIMD512
!!#undef  SIMD256
!!#elif defined(__AVX__) || defined(__AVX2__)
!!#define SIMD256
!!#undef  SIMD512
!!#endif

#if   defined(SIMD512)
#define SIMD_WIDTH_64BIT_NATIVE     8
#define SIMD_WIDTH_32BIT_NATIVE     16
#undef  memalign32
#elif defined(SIMD256)
#define SIMD_WIDTH_64BIT_NATIVE     4
#define SIMD_WIDTH_32BIT_NATIVE     8
#undef  memalign64
#else
#define SIMD_WIDTH_64BIT_NATIVE     1
#define SIMD_WIDTH_32BIT_NATIVE     1
#endif

#if defined(SIMD512) || defined(SIMD256)
#define SIMD_OPTIMIZED
#define SIMD_TYPES
#else
#undef  SIMD_OPTIMIZED
#undef  SIMD_TYPES
#endif

#if defined(SIMD512) || defined(SIMD256)
#define SIMD_WIDTH_64BIT_LOGICAL_2X 2 * SIMD_WIDTH_64BIT_NATIVE
#define SIMD_WIDTH_32BIT_LOGICAL_2X 2 * SIMD_WIDTH_32BIT_NATIVE
#define SIMD_WIDTH_64BIT_LOGICAL_4X 4 * SIMD_WIDTH_64BIT_NATIVE
#define SIMD_WIDTH_32BIT_LOGICAL_4X 4 * SIMD_WIDTH_32BIT_NATIVE
#else
#define SIMD_WIDTH_64BIT_LOGICAL_2X 1 * SIMD_WIDTH_64BIT_NATIVE
#define SIMD_WIDTH_32BIT_LOGICAL_2X 1 * SIMD_WIDTH_32BIT_NATIVE
#define SIMD_WIDTH_64BIT_LOGICAL_4X 1 * SIMD_WIDTH_64BIT_NATIVE
#define SIMD_WIDTH_32BIT_LOGICAL_4X 1 * SIMD_WIDTH_32BIT_NATIVE
#endif

#define SIMD_WIDTH_64BIT            SIMD_WIDTH_64BIT_LOGICAL_2X

#ifndef GNU_SIMD_EXPERIMENTAL
#define SIMD_EXP                    EXP
#define SIMD_LOG                    LOG
#define SIMD_POW(X_1, X_2)          ((X_1) ** (X_2))
#endif

#if defined(_OPENMP45)
#define __simd_static               simd:static
#define OMP_SCHEDULE_SIMD_STATIC    simd:static
#define OMP_SIMDLEN(x)              simdlen(x)
#else
#define __simd_static               static
#define OMP_SCHEDULE_SIMD_STATIC    static
#define OMP_SIMDLEN(x)
#endif
