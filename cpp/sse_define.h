#ifndef SSE_DEFINE_H_
#define SSE_DEFINE_H_

#if defined(HAVE_SSE) && !defined(HAVE_AVX)
    #include <xmmintrin.h>
    #include <pmmintrin.h>
    #define MVAR  __m128
    #define MVARi __m128i
    #define MSTEP 4
    #define MEM_ALIGN_SIZE 16 // in bytes
    #define MM_set1_ps      _mm_set1_ps
    #define MM_load_ps      _mm_load_ps
    #define MM_loadu_ps     _mm_loadu_ps
    #define MM_store_ps     _mm_store_ps
    #define MM_storeu_ps    _mm_storeu_ps
    #define MM_stream_ps    _mm_stream_ps
    #define MM_mul_ps       _mm_mul_ps
    #define MM_div_ps       _mm_div_ps
    #define MM_sub_ps       _mm_sub_ps
    #define MM_add_ps       _mm_add_ps
    #define MM_min_ps       _mm_min_ps
    #define MM_max_ps       _mm_max_ps
    #define MM_add_ps       _mm_add_ps
    #define MM_cmpgt_ps     _mm_cmpgt_ps 
    #define MM_cmplt_ps     _mm_cmplt_ps
    #define MM_or_ps        _mm_or_ps
    #define MM_xor_ps       _mm_xor_ps
    #define MM_and_ps       _mm_and_ps
    #define MM_andnot_ps    _mm_andnot_ps
    #define MM_load_si(x)   _mm_load_si128((__m128i*)(x))

#elif defined(HAVE_AVX)
    #include <immintrin.h>
    #define MVAR  __m256
    #define MVARi __m256i
    #define MSTEP 8
    #define MEM_ALIGN_SIZE 32 // in bytes
    #define MM_set1_ps      _mm256_set1_ps
    #define MM_load_ps      _mm256_load_ps
    #define MM_loadu_ps     _mm256_loadu_ps
    #define MM_store_ps     _mm256_store_ps
    #define MM_storeu_ps    _mm256_storeu_ps
    #define MM_stream_ps    _mm256_stream_ps
    #define MM_mul_ps       _mm256_mul_ps
    #define MM_div_ps       _mm256_div_ps
    #define MM_sub_ps       _mm256_sub_ps
    #define MM_add_ps       _mm256_add_ps
    #define MM_min_ps       _mm256_min_ps
    #define MM_max_ps       _mm256_max_ps
    #define MM_add_ps       _mm256_add_ps
    #define MM_cmpgt_ps     _mm256_cmpgt_ps 
    #define MM_cmplt_ps     _mm256_cmplt_ps
    #define MM_or_ps        _mm256_or_ps
    #define MM_xor_ps       _mm256_xor_ps
    #define MM_and_ps       _mm256_and_ps
    #define MM_andnot_ps    _mm256_andnot_ps
    #define MM_load_si(x)   _mm256_load_si256((__m256i*)(x))
    #define SIGN_VEC SIGN_VEC

    static inline __m256  __attribute__((always_inline)) MM_cmplt_ps(__m256 a, __m256 b) { 
        return _mm256_cmp_ps(a, b, _CMP_LT_OQ); }
    static inline __m256  __attribute__((always_inline)) MM_cmpgt_ps(__m256 a, __m256 b) { 
        return _mm256_cmp_ps(a, b, _CMP_GT_OQ); }
    static inline __m256  __attribute__((always_inline)) MM_cmpeq_ps(__m256 a, __m256 b) { 
        return _mm256_cmp_ps(a, b, _CMP_EQ_OQ); }
#endif

#if defined(HAVE_SSE) || defined(HAVE_AVX)
#define HAVE_VECTOR_MATH
#endif 
 
#endif // SSE_define.h