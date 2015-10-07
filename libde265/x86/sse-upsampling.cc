/*
 * H.265 video codec.
 * Copyright (c) 2013 openHEVC contributors
 * Copyright (c) 2013-2014 struktur AG, Dirk Farin <farin@struktur.de>
 *
 * This file is part of libde265.
 *
 * libde265 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * libde265 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libde265.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <emmintrin.h>
#include <tmmintrin.h> // SSSE3
#if HAVE_SSE4_1
#include <smmintrin.h>
#endif

#include "sse-upsampling.h"
#include "libde265/util.h"
#include <assert.h>
#include <stdlib.h>

/*

Possible operations to use:
- _mm_madd_epi16(a,b) : 
r0 := (a0 * b0) + (a1 * b1)
r1 := (a2 * b2) + (a3 * b3)
r2 := (a4 * b4) + (a5 * b5)
r3 := (a6 * b6) + (a7 * b7)

- _mm_insert_epi16(a,b,im) : Inserts the least significant 16 bits of b into the selected 16-bit integer of a
r0 := (imm == 0) ? b : a0;
r1 := (imm == 1) ? b : a1;
...
r7 := (imm == 7) ? b : a7;

- _mm_maddubs_epi1 6_mm_maddubs_epi16: multiplies and adds integers. (16*8bit unsigned, 16*8bit signed) (SSSE3)
r0 := SATURATE_16((a0 * b0) + (a1 * b1))
r1 := SATURATE_16((a2 * b2) + (a3 * b3))
...
r7 := SATURATE_16((a14 * b14) + (a15 * b15))

- _mm_mullo_epi16 : Multiplies the 8 signed or unsigned 16-bit integers from a by the 8 signed or unsigned 16-bit integers from b.
r0 := (a0 * b0)[15:0]
r1 := (a1 * b1)[15:0]
...
r7 := (a7 * b7)[15:0]

- _mm_shuffle_epi8 : This instruction shuffles 16-byte parameters from a 128-bit parameter.
  Maybe for resorting something ...

// horizontal add of four 32 bit partial sums and return result
// Two steps to add 4 32 bit counters
vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));

- _mm_cvtsi128_si32: Moves the least significant 32 bits of a to a 32-bit integer.

- _mm_setr_epi64: Set the 2 __m64 values

- SSE4: _mm_insert_epi8 (Insert byte into 128bit vector)
- SSE4: _mm_extract_epi8 (Extract 8bit)
- SSE4: _mm_blendv_epi8 (Reorder from a,b to output per byte)

- SSE2: _mm_slli_epi64  (Shift left 2 64 bit values)
- SSE2: _mm_unpackhi_epi8 Interleaves 8 bit wise


*/

#ifdef HAVE_SSE4_1
#define MEMORY_PADDING  8
#else
#define MEMORY_PADDING  0
#endif

#define STANDARD_ALIGNMENT 16
#define MAX_CU_SIZE 64

#ifdef HAVE___MINGW_ALIGNED_MALLOC
#define ALLOC_ALIGNED(alignment, size)         __mingw_aligned_malloc((size), (alignment))
#define FREE_ALIGNED(mem)                      __mingw_aligned_free((mem))
#elif _WIN32
#define ALLOC_ALIGNED(alignment, size)         _aligned_malloc((size), (alignment))
#define FREE_ALIGNED(mem)                      _aligned_free((mem))
#elif defined(HAVE_POSIX_MEMALIGN)
static inline void *ALLOC_ALIGNED(size_t alignment, size_t size) {
	void *mem = NULL;
	if (posix_memalign(&mem, alignment, size) != 0) {
		return NULL;
	}
	return mem;
};
#define FREE_ALIGNED(mem)                      free((mem))
#else
#define ALLOC_ALIGNED(alignment, size)      memalign((alignment), (size))
#define FREE_ALIGNED(mem)                   free((mem))
#endif

#define ALLOC_ALIGNED_16(size)              ALLOC_ALIGNED(16, size)

// H.8.1.4.1.3 Derivation process for reference layer sample location in units of 1/16-th sample
// The position_params array contains the precomputed values needed for this. (H 63)
#define MAP_X_TO_REF_16(xP) (((xP - position_params[0]) * position_params[4] + position_params[6] + (1 << 11)) >> 12) + position_params[2]
#define MAP_Y_TO_REF_16(yP) (((yP - position_params[1]) * position_params[5] + position_params[7] + (1 << 11)) >> 12) + position_params[3]
// Get xRef in samples positions and the phase (H 63) plus (H 29) and (H 30).
#define MAP_X_TO_REF_PHASE(xP, xRef, xPhase) int xRef16 = MAP_X_TO_REF_16(xP); xRef = xRef16 >> 4; xPhase = xRef16 % 16;
#define MAP_Y_TO_REF_PHASE(yP, yRef, yPhase) int yRef16 = MAP_Y_TO_REF_16(yP); yRef = yRef16 >> 4; yPhase = yRef16 % 16;

//// Table H.1 – 16-phase luma resampling filter
//int fL[16][8] = { { 0, 0,   0, 64,  0,   0, 0,  0},
//                  { 0, 1,  -3, 63,  4,  -2, 1,  0},
//                  {-1, 2,  -5, 62,  8,  -3, 1,  0},
//                  {-1, 3,  -8, 60, 13,  -4, 1,  0},
//                  {-1, 4, -10, 58, 17,  -5, 1,  0},
//                  {-1, 4, -11, 52, 26,  -8, 3, -1},
//                  {-1, 3,  -9, 47, 31, -10, 4, -1},
//                  {-1, 4, -11, 45, 34, -10, 4, -1},
//                  {-1, 4, -11, 40, 40, -11, 4, -1},
//                  {-1, 4, -10, 34, 45, -11, 4, -1},
//                  {-1, 4, -10, 31, 47,  -9, 3, -1},
//                  {-1, 3,  -8, 26, 52, -11, 4, -1},
//                  { 0, 1,  -5, 17, 58, -10, 4, -1},
//                  { 0, 1,  -4, 13, 60,  -8, 3, -1},
//                  { 0, 1,  -3,  8, 62,  -5, 2, -1},
//                  { 0, 1,  -2,  4, 63,  -3, 1,  0} };

// 8 bit representation of the above filters. Each 8 bits are one coefficient.
//ALIGNED_16(const __m64) fL_8bit[16] = { 0x0000004000000000,
//                                        0x0001fd3f04fe0100,
//                                        0xff02fb3e08fd0100,
//                                        0xff03f83c0dfc0100,
//                                        0xff04f63a11fb0100,
//                                        0xff04f5341af803ff,
//                                        0xff03f72f1ff604ff,
//                                        0xff04f52d22f604ff,
//                                        0xff04f52828f504ff,
//                                        0xff04f6222df504ff,
//                                        0xff04f61f2ff703ff,
//                                        0xff03f81a34f504ff,
//                                        0x0001fb113af604ff,
//                                        0x0001fc0d3cf803ff,
//                                        0x0001fd083efb02ff,
//                                        0x0001fe043ffd0100 };

// The filters in 8 bit fromat. 2 times the same filter to perform 2 multiplications with the filter at once.
ALIGNED_16(const int8_t) fL_8bit_2x[256] = {  0, 0,   0, 64,  0,   0, 0,  0 ,  0, 0,   0, 64,  0,   0, 0,  0 ,
                                              0, 1,  -3, 63,  4,  -2, 1,  0 ,  0, 1,  -3, 63,  4,  -2, 1,  0 ,
                                             -1, 2,  -5, 62,  8,  -3, 1,  0 , -1, 2,  -5, 62,  8,  -3, 1,  0 ,
                                             -1, 3,  -8, 60, 13,  -4, 1,  0 , -1, 3,  -8, 60, 13,  -4, 1,  0 ,
                                             -1, 4, -10, 58, 17,  -5, 1,  0 , -1, 4, -10, 58, 17,  -5, 1,  0 ,
                                             -1, 4, -11, 52, 26,  -8, 3, -1 , -1, 4, -11, 52, 26,  -8, 3, -1 ,
                                             -1, 3,  -9, 47, 31, -10, 4, -1 , -1, 3,  -9, 47, 31, -10, 4, -1 ,
                                             -1, 4, -11, 45, 34, -10, 4, -1 , -1, 4, -11, 45, 34, -10, 4, -1 ,
                                             -1, 4, -11, 40, 40, -11, 4, -1 , -1, 4, -11, 40, 40, -11, 4, -1 ,
                                             -1, 4, -10, 34, 45, -11, 4, -1 , -1, 4, -10, 34, 45, -11, 4, -1 ,
                                             -1, 4, -10, 31, 47,  -9, 3, -1 , -1, 4, -10, 31, 47,  -9, 3, -1 ,
                                             -1, 3,  -8, 26, 52, -11, 4, -1 , -1, 3,  -8, 26, 52, -11, 4, -1 ,
                                              0, 1,  -5, 17, 58, -10, 4, -1 ,  0, 1,  -5, 17, 58, -10, 4, -1 ,
                                              0, 1,  -4, 13, 60,  -8, 3, -1 ,  0, 1,  -4, 13, 60,  -8, 3, -1 ,
                                              0, 1,  -3,  8, 62,  -5, 2, -1 ,  0, 1,  -3,  8, 62,  -5, 2, -1 ,
                                              0, 1,  -2,  4, 63,  -3, 1,  0 ,  0, 1,  -2,  4, 63,  -3, 1,  0 };

// The filters in 16 bit format. Every line is one filter.
ALIGNED_16(const int16_t) fL_16bit[128] = {  0, 0,   0, 64,  0,   0, 0,  0 ,
                                             0, 1,  -3, 63,  4,  -2, 1,  0 ,
                                            -1, 2,  -5, 62,  8,  -3, 1,  0 ,
                                            -1, 3,  -8, 60, 13,  -4, 1,  0 ,
                                            -1, 4, -10, 58, 17,  -5, 1,  0 ,
                                            -1, 4, -11, 52, 26,  -8, 3, -1 ,
                                            -1, 3,  -9, 47, 31, -10, 4, -1 ,
                                            -1, 4, -11, 45, 34, -10, 4, -1 ,
                                            -1, 4, -11, 40, 40, -11, 4, -1 ,
                                            -1, 4, -10, 34, 45, -11, 4, -1 ,
                                            -1, 4, -10, 31, 47,  -9, 3, -1 ,
                                            -1, 3,  -8, 26, 52, -11, 4, -1 ,
                                             0, 1,  -5, 17, 58, -10, 4, -1 ,
                                             0, 1,  -4, 13, 60,  -8, 3, -1 ,
                                             0, 1,  -3,  8, 62,  -5, 2, -1 ,
                                             0, 1,  -2,  4, 63,  -3, 1,  0 };

// Measure the execution time of the upsampling and print it to stdout
#define MEASURE_EXECUTION_TIME 0
// Profile by getting CPU ticks and running the upsampling function multiple times
#define PROFILE_FUNCTION 0

#if MEASURE_EXECUTION_TIME
#define INTMAX_MAX 9223372036854775807LL
#include <chrono>
using namespace std;
using namespace std::chrono;
#endif
#if PROFILE_FUNCTION
#ifdef _WIN32
#include <intrin.h>
#include <inttypes.h>
uint64_t rdtsc(){
    return __rdtsc();
}
#else
// Linux
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#endif
#endif

#if HAVE_SSE4_1
// Upsampling process for 8 bit input ...
void resampling_process_of_luma_sample_values_sse(uint8_t *src, ptrdiff_t srcstride, int src_size[2],
                                                  uint8_t *dst, ptrdiff_t dststride, int dst_size[2],
                                                  int position_params[10])
{
#if MEASURE_EXECUTION_TIME
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif
#if PROFILE_FUNCTION

  // Measure the operations (loading/saving/arithmetic)
  uint64_t lastTicks1;
  uint64_t loadTicks_ver = 0;
  uint64_t loadTicks_hor = 0;
  uint64_t arithmeticTicks_ver = 0;
  uint64_t arithmeticTicks_hor = 0;
  uint64_t saveTicks_ver = 0;
  uint64_t saveTicks_hor = 0;

  // Measure the parts (memAlloc/Hor/Ver)
  uint64_t lastTicksPart;
  uint64_t horizontalTicks = 0;
  uint64_t verticalTicks = 0;
  uint64_t allocBuffer = 0;
    
  uint64_t beforeWholeFunction = rdtsc();
  // Run this function 100 times
  for (int run=0; run<100; run++) {
#endif

  static int16_t* s_tmp = NULL;  // Init static pointer to tempraray array for upsampling
  static int s_tmp_size = -1;

  // Reference layer size
  int PicHeightInSamplesRefLayerY = src_size[1];
  int PicWidthInSamplesRefLayerY  = src_size[0];
  assert(PicHeightInSamplesRefLayerY % 2 == 0);   // Pic height must be dividible by 2 for this function to work. This is a requirement of HEVC so this should always be the case.

  // Current layer size
  int PicHeightInSamplesCurLayerY = dst_size[1];
  int PicWidthInSamplesCurLayerY  = dst_size[0];

  int BitDepthRefLayerY = position_params[8];
  int BitDepthCurrY     = position_params[9];
  int clipMax           = (1 << BitDepthCurrY) - 1;

  // 4. The variables shift1, shift2 and offset are derived as follows:
  int shift2 = 20 - BitDepthCurrY;     // (H 34)
  int offset = 1 << (shift2 - 1);      // (H 35)

  // Perform horizontal / vertical upsampling seperately.
  
  // Allocate temporaray buffer
#if PROFILE_FUNCTION
  lastTicksPart = rdtsc();
#endif
  int alignment = STANDARD_ALIGNMENT;
  int tmpStride = ((PicWidthInSamplesCurLayerY + alignment-1) / alignment * alignment); // In bytes
  int tmpHeight = PicHeightInSamplesRefLayerY;
  if (s_tmp_size == -1 && s_tmp == NULL) {
    // Allocate the temporaray tmp buffer
    s_tmp = (int16_t *)ALLOC_ALIGNED_16(tmpStride * tmpHeight * sizeof(int16_t) + MEMORY_PADDING);
    s_tmp_size = tmpStride * tmpHeight * sizeof(int16_t) + MEMORY_PADDING;
  }
  else if (s_tmp_size < tmpStride * tmpHeight * sizeof(int16_t) + MEMORY_PADDING) {
    // The buffer is not big enough. TODO (Maybe a version of realloc would be better??)
    // Allocate buffer that is big enough
    FREE_ALIGNED(s_tmp);
    s_tmp = (int16_t *)ALLOC_ALIGNED_16(tmpStride * tmpHeight * sizeof(int16_t) + MEMORY_PADDING);
    s_tmp_size = tmpStride * tmpHeight * sizeof(int16_t) + MEMORY_PADDING;
  }
  int16_t* tmp = s_tmp;
#if PROFILE_FUNCTION
  allocBuffer += rdtsc() - lastTicksPart;
#endif

  // -------- Horizontal upsampling ------------ //
#if PROFILE_FUNCTION
  lastTicksPart = rdtsc();
#endif
  __m128i in;           // Contains 16 input sample values that we will perform the multiplication on
  __m128i filter;       // Contains two filters
  __m128i result;       // Result of multiplication

  in = _mm_setzero_si128(); // Set in to 0
  
  // Get the number of samples we have to pad
  int xRef, xPhase;
  MAP_X_TO_REF_PHASE(0, xRef, xPhase);
  xRef -= 3;  // For the filter (3 coefficients on the left)
  assert(xRef > -8 && xRef <= 0); // No shift for more than 8 bytes. (Is this even possible?)
  int shift = abs(xRef);
  
  // Get pointers to y lines
  uint8_t *yLine0 = src;
  uint8_t *yLine1 = src + srcstride;
  int16_t *yLineOut0 = tmp;
  int16_t *yLineOut1 = tmp + tmpStride;
  for (int y = 0; y < PicHeightInSamplesRefLayerY; y += 2) {    // Always process pairs of y lines

    int xReadPos = 0;
    
#if PROFILE_FUNCTION
    lastTicks1 = rdtsc();
#endif
    // Fill the input buffer (in) and apply padding.
    int byte0 = *yLine0;    // Get the left most byte of this and the next line
    int byte1 = *yLine1;
    int iByteCount;
    for (iByteCount = 0; iByteCount < shift; iByteCount++) { 
      // shift left
      in = _mm_slli_epi64(in, 8);

      // Insert the left most byte 'shift' times
      in = _mm_insert_epi8(in, byte0, 8);
      in = _mm_insert_epi8(in, byte1, 0);
    }
    
    // Fill the rest of the buffer with bytes from input
    for (xReadPos = 0; iByteCount < 8; iByteCount++) {
      // Get bytes
      byte0 = *(yLine0 + xReadPos);
      byte1 = *(yLine1 + xReadPos++);

      // shift left
      in = _mm_slli_epi64(in, 8);

      // insert bytes
      in = _mm_insert_epi8(in, byte0, 8);
      in = _mm_insert_epi8(in, byte1, 0);
    }
    int xRefBuf = 0; // The buffer is now ready to start at xRef of 0
#if PROFILE_FUNCTION
    loadTicks_hor += rdtsc() - lastTicks1;
#endif
    
    // The input buffer is set up. Here we go:
    for (int xP=0; xP < PicWidthInSamplesCurLayerY; xP++) {
      // Map from this layer x to refLayer x
      MAP_X_TO_REF_PHASE(xP, xRef, xPhase);

#if PROFILE_FUNCTION
      lastTicks1 = rdtsc();
#endif
      while (xRef != xRefBuf) {
        // We have to load bytes until xRefBuf and xRef match again
        if (xReadPos < PicWidthInSamplesRefLayerY) {  // Otherwise repeat last byte (padding)
          // Get bytes
          byte0 = *(yLine0 + xReadPos);
          byte1 = *(yLine1 + xReadPos++);
        }

        // shift left
        in = _mm_slli_epi64(in, 8);

        // insert bytes
        in = _mm_insert_epi8(in, byte0, 8);
        in = _mm_insert_epi8(in, byte1, 0);

        xRefBuf++;
      }
#if PROFILE_FUNCTION
      loadTicks_hor += rdtsc() - lastTicks1;
#endif

#if PROFILE_FUNCTION
      lastTicks1 = rdtsc();
#endif
      // Get the correct filter (according to xPhase)
      filter = _mm_load_si128((__m128i *)(fL_16bit + xPhase * 8));

      // multiply
      result = _mm_maddubs_epi16(in, filter);
      
      // Add up results
      result = _mm_add_epi16(result, _mm_srli_si128(result, 4));
      result = _mm_add_epi16(result, _mm_srli_si128(result, 2));

#if PROFILE_FUNCTION
      arithmeticTicks_hor += rdtsc() - lastTicks1;
      lastTicks1 = rdtsc();
#endif

      // Save results
      *(yLineOut0 + xP ) = _mm_extract_epi16(result, 4);
      *(yLineOut1 + xP ) = _mm_extract_epi16(result, 0);

#if PROFILE_FUNCTION
      saveTicks_hor += rdtsc() - lastTicks1;
#endif
    }

    // Go to next y lines
    yLine0 += srcstride * 2;
    yLine1 += srcstride * 2;
    yLineOut0 += tmpStride * 2;
    yLineOut1 += tmpStride * 2;
  }

  // -------- Vertical upsampling ------------ ////
#if PROFILE_FUNCTION
  horizontalTicks += rdtsc() - lastTicksPart;
  lastTicksPart = rdtsc();
#endif
  int yRef, yPhase;
  MAP_Y_TO_REF_PHASE(0, yRef, yPhase);
  yRef -= 3;  // For the filter (3 coefficients on the top)
  assert(yRef > -8 && yRef <= 0); // No shift for more than 8 bytes. (Is this even possible?)
  shift = abs(yRef);

  int16_t *yLineIn;
  uint8_t  *yLineOut;
  for (int xP=0; xP < PicWidthInSamplesCurLayerY; xP++) {
    // Get top y buffer for x position
    yLineIn = tmp + xP;
    yLineOut = dst + xP;
    
    int yReadPos = 0;

    // Fill the input buffer (in) and apply padding.
#if PROFILE_FUNCTION
    lastTicks1 = rdtsc();
#endif
    int iByteCount;
    int byte = *(yLineIn); // Read top val

    for (iByteCount = 0; iByteCount < shift; iByteCount++) { 
      // shift left
      in = _mm_slli_si128(in, 2);

      // Insert the top most byte 'shift' times
      in = _mm_insert_epi16(in, byte, 0);
    }
    
    // Fill the rest of the buffer with bytes from input
    for (; iByteCount < 8; iByteCount++) {
      // Get bytes
      byte = *(yLineIn); yLineIn += tmpStride;  // Read val and goto next y line
      yReadPos++;

      // shift left
      in = _mm_slli_si128(in, 2);

      // insert bytes
      in = _mm_insert_epi16(in, byte, 0);
    }
    int yRefBuf = 0; // The buffer is now ready to start at yRef of 0
#if PROFILE_FUNCTION
    loadTicks_ver += rdtsc() - lastTicks1;
#endif

    // The input buffer is set up. Here we go:
    for (int y = 0; y < PicHeightInSamplesCurLayerY; y++) {
      // Map from this layer y to refLayer y
      MAP_X_TO_REF_PHASE(y, yRef, yPhase);

#if PROFILE_FUNCTION
      lastTicks1 = rdtsc();
#endif
      while (yRef != yRefBuf) {
        // We have to load bytes until yRefBuf and yRef match again
        if (yReadPos < PicHeightInSamplesRefLayerY) {  // Otherwise repeat last byte (padding)
          // Get bytes
         byte = *(yLineIn); yLineIn += tmpStride;  // Read val and goto next y line
         yReadPos++;
        }

        // shift left
        in = _mm_slli_si128(in, 2);

        // insert bytes
        in = _mm_insert_epi16(in, byte, 0);
        
        yRefBuf++;
      }
#if PROFILE_FUNCTION
      loadTicks_ver += rdtsc() - lastTicks1;
      lastTicks1 = rdtsc();
#endif

      // Get the correct filter (according to yPhase)
      filter = _mm_load_si128((__m128i *)(fL_16bit + xPhase * 8));

      // multiply (16 bit input, out 32bit)
      result = _mm_madd_epi16 (in, filter);

      // Add up results
      result = _mm_add_epi32(result, _mm_srli_si128(result, 8));
      result = _mm_add_epi32(result, _mm_srli_si128(result, 4));
#if PROFILE_FUNCTION
      arithmeticTicks_ver += rdtsc() - lastTicks1;
      lastTicks1 = rdtsc();
#endif
            
      // Save results
      *(yLineOut) = Clip3(0, 255, (( _mm_extract_epi32(result, 0) + offset ) >> shift2) );
      yLineOut += dststride;

#if PROFILE_FUNCTION
      saveTicks_ver += rdtsc() - lastTicks1;
#endif

    }
  }

#if PROFILE_FUNCTION
  verticalTicks += rdtsc() - lastTicksPart;
  }

  uint64_t wholeTime = rdtsc() - beforeWholeFunction;
  double d_wholeTime = wholeTime;
  double p1 = allocBuffer / d_wholeTime * 100;
  double p2 = horizontalTicks / d_wholeTime * 100;
  double p3 = verticalTicks / d_wholeTime * 100;
  // Print values
  printf("Mem Alloc: %" PRIu64 " (%f)\n", allocBuffer, p1);
  printf("Hor Filt : %" PRIu64 " (%f)\n", horizontalTicks, p2);
  printf("Ver Filt : %" PRIu64 " (%f)\n", verticalTicks, p3);
  printf("\n");

  uint64_t verSum = loadTicks_ver + arithmeticTicks_ver + saveTicks_ver;
  uint64_t horSum = loadTicks_hor + arithmeticTicks_hor + saveTicks_hor;
  double d_verSum = verSum;
  double d_horSum = horSum;
  double p4 = loadTicks_ver / d_verSum * 100;
  double p5 = arithmeticTicks_ver / d_verSum * 100;
  double p6 = saveTicks_ver / d_verSum * 100;
  double p7 = loadTicks_hor / d_horSum * 100;
  double p8 = arithmeticTicks_hor / d_horSum * 100;
  double p9 = saveTicks_hor / d_horSum * 100;

  printf("Load  Ver: %" PRIu64 " (%f)\n", loadTicks_ver, p4);
  printf("Load  Hor: %" PRIu64 " (%f)\n", loadTicks_hor, p7);
  printf("Arith Ver: %" PRIu64 " (%f)\n", arithmeticTicks_ver, p5);
  printf("Arith Hor: %" PRIu64 " (%f)\n", arithmeticTicks_hor, p8);
  printf("Save  Ver: %" PRIu64 " (%f)\n", saveTicks_ver, p6);
  printf("Save  Hor: %" PRIu64 " (%f)\n", saveTicks_hor, p9);
#endif
#if MEASURE_EXECUTION_TIME
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
  printf("Upsampling Y from (%dx%d) to (%dx%d) took %d us\n", PicWidthInSamplesRefLayerY, PicHeightInSamplesRefLayerY, PicWidthInSamplesCurLayerY, PicHeightInSamplesCurLayerY, duration);
#endif
}

void resampling_process_of_luma_block_sse_8bit(const uint8_t *src, ptrdiff_t src_stride, int16_t src_height,
                                               int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                               int x_dst, int y_dst, const int *position_params, int16_t *buffer)
{
  int BitDepthRefLayerY = position_params[8];
  int BitDepthCurrY     = position_params[9];
  
  // As the function name indicates this function can only upsample from 8bit to 8bit
  assert(BitDepthRefLayerY == 8);
  assert(BitDepthCurrY == 8);

  int clipMax           = (1 << (BitDepthCurrY)) - 1;

  int xRef16, xRefBuf, xRef, xPhase, xP;
  int yRef16, yRefBuf, yRef, yPhase, yP;
  const uint8_t  *rlPicSampleL;
  int16_t *rsLumaSample;
  int16_t  *tmpSample;

  // Calculate the position of the top left point in the reference
  int x_src = ((((x_dst - position_params[0]) * position_params[4] + position_params[6] + (1 << 11)) >> 12) + position_params[2]) >> 4;  // (H 63)
  int y_src = ((((y_dst - position_params[1]) * position_params[5] + position_params[7] + (1 << 11)) >> 12) + position_params[3]) >> 4;  // (H 64)

  assert(src_height >= 1);
  int processing_height = src_height + 3 + 4; // padding (3 on left, 4 on right)
  // Process in y blocks of 8, then process the remainder

  __m128i r0, r1, r2, r3, r4, r5, r6, r7, r8;
  __m128i filter;

  // Horizontal filtering
  for (int x=0; x < dst_width; x++) {
    xP = x_dst + x;

    // 1.
    // H.8.1.4.1.3 Derivation process for reference layer sample location in units of 1/16-th sample
    // The position_params array contains the precomputed values needed for this.
    xRef16 = (((xP - position_params[0]) * position_params[4] + position_params[6] + (1 << 11)) >> 12) + position_params[2];  // (H 63)
        
    // 2. The variables xRef and xPhase are derived as follows:
    xRef   = xRef16 >> 4;  // (H 29)
    xPhase = xRef16 % 16;  // (H 30)
    
    // Get pointers to source and destination
    xRefBuf = xRef - x_src;
    rlPicSampleL = src - 3*src_stride - 3 + xRefBuf;
    tmpSample    = buffer + x * (MAX_CU_SIZE+8);
    
    // Load the correct filter accoding to xPhase
    filter = _mm_load_si128((__m128i *)(fL_8bit_2x + xPhase*16));

    int y, outAddrOffset=0;
    for (y=processing_height; y >= 8; y-=8, outAddrOffset+=8) {

      // Load 8 y lines
      
      // Load 2 8 bit values from two lines and pack 2 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_unpacklo_epi64(r2, r3);

      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r4 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r3 = _mm_unpacklo_epi64(r3, r4);

      // Multiply the 4 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      r2 = _mm_maddubs_epi16(r2, filter);
      r3 = _mm_maddubs_epi16(r3, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r2 = _mm_hadd_epi16(r2, r3);
      r0 = _mm_hadd_epi16(r0, r2);

      // Save 8 resulting values to tmp
      _mm_store_si128((__m128i*) (tmpSample+outAddrOffset), r0);
    }

    // Process the remaining y lines (y lines remaining)
    if (y == 7) {
      // Load 7 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_unpacklo_epi64(r2, r3);

      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      
      // Multiply the 4 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      r2 = _mm_maddubs_epi16(r2, filter);
      r3 = _mm_maddubs_epi16(r3, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r2 = _mm_hadd_epi16(r2, r3);
      r0 = _mm_hadd_epi16(r0, r2);

      // Save 7 resulting values to tmp (last value is 0 but it is written anyways)
      _mm_store_si128((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 6) {
      // Load 6 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_unpacklo_epi64(r2, r3);

      // Multiply the 4´3 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      r2 = _mm_maddubs_epi16(r2, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r2 = _mm_hadd_epi16(r2, _mm_setzero_si128());
      r0 = _mm_hadd_epi16(r0, r2);

      // Save 6 resulting values to tmp (last 2 values are 0 but it is written anyways)
      _mm_store_si128((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 5) {
      // Load 5 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;

      // Multiply the 4´3 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      r2 = _mm_maddubs_epi16(r2, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r2 = _mm_hadd_epi16(r2, _mm_setzero_si128());
      r0 = _mm_hadd_epi16(r0, r2);

      // Save 5 resulting values to tmp (last 3 values are 0 but it is written anyways)
      _mm_store_si128((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 4) {
      // Load 4 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_unpacklo_epi64(r1, r2);

      // Multiply the 2 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      
      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 4 resulting values to tmp
      _mm_storel_epi64((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 3) {
      // Load 3 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      
      // Multiply the 2 buffers by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      
      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 3 resulting values to tmp (last value is 0 but it is written anyways)
      _mm_storel_epi64((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 2) {
      // Load 2 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, r1);

      // Multiply the buffer by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      
      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 2 resulting values to tmp (last 2 values are 0 but it is written anyways)
      _mm_storel_epi64((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else if (y == 1) {
      // Load 1 y lines
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleL);    rlPicSampleL += src_stride;
      r0 = _mm_unpacklo_epi64(r0, _mm_setzero_si128());
      
      // Multiply the buffer by the double filter
      r0 = _mm_maddubs_epi16(r0, filter);
      
      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 1 resulting value to tmp (last 3 values are 0 but it is written anyways)
      _mm_storel_epi64((__m128i*) (tmpSample+outAddrOffset), r0);
    }
    else
      assert(y==0);
  }

  // Horizontal upsampling done. Perform vertical upsampling.
  // The temporary buffer (buffer) is flipped (x and y switched).
  // In this second filtering we will just perform the same flipping again.
  
  // 4. The variables shift1, shift2 and offset are derived as follows:
  //int shift2 = 12;     // (20 - BitDepthCurrY)  (H 34)
  //int offset = 2048;   // (1 << (shift2 - 1))   (H 35)
  __m128i offset = _mm_set_epi32(2048, 2048, 2048, 2048);
  // To clipt the output to 0...255
  __m128i clip8 = _mm_set_epi32(255, 255, 255, 255);

  int16_t *out;
  for (int y = 0; y < dst_height; y++) {
    yP = y_dst + y;

    // 1.
    // H.8.1.4.1.3 Derivation process for reference layer sample location in units of 1/16-th sample
    // The position_params array contains the precomputed values needed for this.
    yRef16 = (((yP - position_params[1]) * position_params[5] + position_params[7] + (1 << 11)) >> 12) + position_params[3];  // (H 64)

    // 3. The variables yRef and yPhase are derived as follows:
    yPhase = yRef16 % 16;  // (H 32)
    yRef   = yRef16 >> 4;  // (H 31)

    
    // Get pointer to source and destination
    yRefBuf = yRef - y_src;
    tmpSample = buffer + yRefBuf;
    out = dst + y * dst_stride;

    
    // Load the correct filter accoding to xPhase
    filter = _mm_load_si128((__m128i *)(fL_16bit + yPhase * 8));

    int x, outAddrOffset = 0;
    for (x = dst_width; x >= 8; x -= 8, outAddrOffset += 8) {

      // Load 8 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r3 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r4 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r5 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r6 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r7 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 8 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16 (r0, filter);
      r1 = _mm_madd_epi16 (r1, filter);
      r2 = _mm_madd_epi16 (r2, filter);
      r3 = _mm_madd_epi16 (r3, filter);
      r4 = _mm_madd_epi16 (r4, filter);
      r5 = _mm_madd_epi16 (r5, filter);
      r6 = _mm_madd_epi16 (r6, filter);
      r7 = _mm_madd_epi16 (r7, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);
      r4 = _mm_hadd_epi32(r4, r5);
      r6 = _mm_hadd_epi32(r6, r7);

      r0 = _mm_hadd_epi32(r0, r2);
      r4 = _mm_hadd_epi32(r4, r6);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r4 = _mm_add_epi32(r4, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r4 = _mm_srai_epi32(r4, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r4 = _mm_min_epi32(r4, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r4 = _mm_slli_epi32(r4, 6);

      r0 = _mm_hadd_epi16(r0, r4);  // The upper 16 bit of each 16 bit result should be 0

      // Save 8 resulting values to dst
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    // Process the remaining x lines (x lines remaining)
    if (x == 7) {
      // Load 7 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r3 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r4 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r5 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r6 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 7 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);
      r4 = _mm_madd_epi16(r4, filter);
      r5 = _mm_madd_epi16(r5, filter);
      r6 = _mm_madd_epi16(r6, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);
      r4 = _mm_hadd_epi32(r4, r5);
      r6 = _mm_hadd_epi32(r6, _mm_setzero_si128());

      r0 = _mm_hadd_epi32(r0, r2);
      r4 = _mm_hadd_epi32(r4, r6);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r4 = _mm_add_epi32(r4, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r4 = _mm_srai_epi32(r4, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r4 = _mm_min_epi32(r4, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r4 = _mm_slli_epi32(r4, 6);

      r0 = _mm_hadd_epi16(r0, r4);  // The upper 16 bit of each 16 bit result should be 0

      // Save 8 resulting values to dst. The last value is 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 6) {
      // Load 6 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r3 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r4 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r5 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 6 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);
      r4 = _mm_madd_epi16(r4, filter);
      r5 = _mm_madd_epi16(r5, filter);
      
      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);
      r4 = _mm_hadd_epi32(r4, r5);
      
      r0 = _mm_hadd_epi32(r0, r2);
      r4 = _mm_hadd_epi32(r4, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r4 = _mm_add_epi32(r4, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r4 = _mm_srai_epi32(r4, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r4 = _mm_min_epi32(r4, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r4 = _mm_slli_epi32(r4, 6);

      r0 = _mm_hadd_epi16(r0, r4);  // The upper 16 bit of each 16 bit result should be 0

      // Save 8 resulting values to dst. The last 2 values are 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 5) {
      // Load 5 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r3 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r4 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 5 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);
      r4 = _mm_madd_epi16(r4, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);
      r4 = _mm_hadd_epi32(r4, _mm_setzero_si128());

      r0 = _mm_hadd_epi32(r0, r2);
      r4 = _mm_hadd_epi32(r4, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r4 = _mm_add_epi32(r4, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r4 = _mm_srai_epi32(r4, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r4 = _mm_min_epi32(r4, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r4 = _mm_slli_epi32(r4, 6);

      r0 = _mm_hadd_epi16(r0, r4);  // The upper 16 bit of each 16 bit result should be 0

      // Save 8 resulting values to dst. The last 3 values are 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 4) {
      // Load 4 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r3 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 4 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);

      r0 = _mm_hadd_epi32(r0, r2);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());  // The upper 16 bit of each 16 bit result should be 0

      // Save 4 resulting values to dst.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);

    }
    else if (x == 3) {
      // Load 3 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r2 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 3 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, _mm_setzero_si128());

      r0 = _mm_hadd_epi32(r0, r2);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());  // The upper 16 bit of each 16 bit result should be 0

      // Save 4 resulting values to dst. The last value is 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 2) {
      // Load 2 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);
      r1 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 2 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r0 = _mm_hadd_epi32(r0, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());  // The upper 16 bit of each 16 bit result should be 0

      // Save 4 resulting values to dst. The last 2 values are 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 1) {
      // Load 1 x lines
      r0 = _mm_loadu_si128((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE + 8);

      // Multiply the 2 buffers by the filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, _mm_setzero_si128());
      r0 = _mm_hadd_epi32(r0, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());  // The upper 16 bit of each 16 bit result should be 0

     // Save 4 resulting values to dst. The last 3 values are 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else
      assert(x==0);
  }
}

// The filters in 8 bit fromat. 4 times the same filter to perform 4 multiplications with the filter at once.
ALIGNED_16(const int8_t) fC_8bit_4x[256] = {
     0, 64,  0,  0,  0, 64,  0,  0,  0, 64,  0,  0,  0, 64,  0,  0,
    -2, 62,  4,  0, -2, 62,  4,  0, -2, 62,  4,  0, -2, 62,  4,  0,
    -2, 58, 10, -2, -2, 58, 10, -2, -2, 58, 10, -2, -2, 58, 10, -2,
    -4, 56, 14, -2, -4, 56, 14, -2, -4, 56, 14, -2, -4, 56, 14, -2,
    -4, 54, 16, -2, -4, 54, 16, -2, -4, 54, 16, -2, -4, 54, 16, -2,
    -6, 52, 20, -2, -6, 52, 20, -2, -6, 52, 20, -2, -6, 52, 20, -2,
    -6, 46, 28, -4, -6, 46, 28, -4, -6, 46, 28, -4, -6, 46, 28, -4,
    -4, 42, 30, -4, -4, 42, 30, -4, -4, 42, 30, -4, -4, 42, 30, -4,
    -4, 36, 36, -4, -4, 36, 36, -4, -4, 36, 36, -4, -4, 36, 36, -4,
    -4, 30, 42, -4, -4, 30, 42, -4, -4, 30, 42, -4, -4, 30, 42, -4,
    -4, 28, 46, -6, -4, 28, 46, -6, -4, 28, 46, -6, -4, 28, 46, -6,
    -2, 20, 52, -6, -2, 20, 52, -6, -2, 20, 52, -6, -2, 20, 52, -6,
    -2, 16, 54, -4, -2, 16, 54, -4, -2, 16, 54, -4, -2, 16, 54, -4,
    -2, 14, 56, -4, -2, 14, 56, -4, -2, 14, 56, -4, -2, 14, 56, -4,
    -2, 10, 58, -2, -2, 10, 58, -2, -2, 10, 58, -2, -2, 10, 58, -2,
     0,  4, 62, -2,  0,  4, 62, -2,  0,  4, 62, -2,  0,  4, 62, -2 };

ALIGNED_16(const int16_t) fC_16bit_2x[128] = {
     0, 64,  0,  0,  0, 64,  0,  0,
    -2, 62,  4,  0, -2, 62,  4,  0,
    -2, 58, 10, -2, -2, 58, 10, -2,
    -4, 56, 14, -2, -4, 56, 14, -2,
    -4, 54, 16, -2, -4, 54, 16, -2,
    -6, 52, 20, -2, -6, 52, 20, -2,
    -6, 46, 28, -4, -6, 46, 28, -4,
    -4, 42, 30, -4, -4, 42, 30, -4,
    -4, 36, 36, -4, -4, 36, 36, -4,
    -4, 30, 42, -4, -4, 30, 42, -4,
    -4, 28, 46, -6, -4, 28, 46, -6,
    -2, 20, 52, -6, -2, 20, 52, -6,
    -2, 16, 54, -4, -2, 16, 54, -4,
    -2, 14, 56, -4, -2, 14, 56, -4,
    -2, 10, 58, -2, -2, 10, 58, -2,
     0,  4, 62, -2,  0,  4, 62, -2 };

void resampling_process_of_chroma_block_sse_8bit(const uint8_t *src, ptrdiff_t src_stride, int16_t src_height,
                                                 int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                                 int x_dst, int y_dst, const int *position_params, int16_t *buffer)
{
  int BitDepthRefLayerC = position_params[8];
  int BitDepthCurrC = position_params[9];

  // As the function name indicates this function can only upsample from 8bit to 8bit
  assert(BitDepthRefLayerC == 8);
  assert(BitDepthCurrC == 8);

  int clipMax = (1 << (BitDepthCurrC)) - 1;

  int xRef16, xRefBuf, xRef, xPhase, xP;
  int yRef16, yRefBuf, yRef, yPhase, yP;
  const uint8_t  *rlPicSampleC;
  int16_t *rsChromaSample;
  int16_t  *tmpSample;

  // Calculate the position of the top left point in the reference
  int x_src = ((((x_dst - position_params[0]) * position_params[4] + position_params[6] + (1 << 11)) >> 12) + position_params[2]) >> 4;  // (H 63)
  int y_src = ((((y_dst - position_params[1]) * position_params[5] + position_params[7] + (1 << 11)) >> 12) + position_params[3]) >> 4;  // (H 64)

  assert(src_height >= 1);
  int processing_height = src_height + 1 + 3; // padding (1 on left, 2 on right)
                                              // Process in y blocks of 8, then process the remainder

  __m128i r0, r1, r2, r3, r4, r5, r6, r7, r8;
  __m128i filter;

  // Horizontal filtering
  for (int x = 0; x < dst_width; x++) {
    xP = x_dst + x;

    // 1.
    // H.8.1.4.1.3 Derivation process for reference layer sample location in units of 1/16-th sample
    // The position_params array contains the precomputed values needed for this.
    xRef16 = (((xP - position_params[0]) * position_params[4] + position_params[6] + (1 << 11)) >> 12) + position_params[2];  // (H 63)

    // 2. The variables xRef and xPhase are derived as follows:
    xRef = xRef16 >> 4;  // (H 29)
    xPhase = xRef16 % 16;  // (H 30)

    // Get pointers to source and destination
    xRefBuf = xRef - x_src;
    rlPicSampleC = src - src_stride - 1 + xRefBuf;
    tmpSample = buffer + x * (MAX_CU_SIZE/2 + 8);

    // Load the correct filter accoding to xPhase
    filter = _mm_load_si128((__m128i *)(fC_8bit_4x + xPhase * 16));

    int y, outAddrOffset = 0;
    for (y = processing_height; y >= 8; y -= 8, outAddrOffset += 8) {

      // Load 8 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);
      
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r3 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_unpacklo_epi32(r2, r3);

      r1 = _mm_unpacklo_epi64(r1, r2);

      // Multiply the 2 buffers by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);
      
      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);
      
      // Save 8 resulting values to tmp
      _mm_store_si128((__m128i*) (tmpSample + outAddrOffset), r0);
    }

    // Process the remaining y lines (y lines remaining)
    if (y == 7) {
      // Load 7 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;

      r1 = _mm_unpacklo_epi64(r1, r2);
      
      // Multiply the 2 buffers by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);

      // Save 8 resulting values to tmp. Last value is ignored later.
      _mm_store_si128((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 6) {
      // Load 6 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      // Multiply the 2 buffers by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);

      // Save 8 resulting values to tmp. Last 2 are ignored later.
      _mm_store_si128((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 5) {
      // Load 5 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;

      // Multiply the 2 buffers by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);
      r1 = _mm_maddubs_epi16(r1, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, r1);

      // Save 8 resulting values to tmp. Last 3 are ignored later.
      _mm_store_si128((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 4) {
      // Load 4 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r2 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r0 = _mm_unpacklo_epi64(r0, r1);

      // Multiply the buffer by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 4 resulting values to tmp.
      _mm_storel_epi64((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 3) {
      // Load 3 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_unpacklo_epi32(r1, r2);

      r0 = _mm_unpacklo_epi64(r0, r1);

      // Multiply the buffer by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 4 resulting values to tmp. Last value will be ignored late.
      _mm_storel_epi64((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 2) {
      // Load 2 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r1 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;
      r0 = _mm_unpacklo_epi32(r0, r1);

      // Multiply the buffer by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 4 resulting values to tmp. Last 2 values will be ignored later.
      _mm_storel_epi64((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else if (y == 1) {
      // Load 2 y lines

      // Load 4 8 bit values from four lines and pack 4 into one 128 bits value
      r0 = _mm_loadl_epi64((__m128i *)rlPicSampleC);    rlPicSampleC += src_stride;

      // Multiply the buffer by the qudruple filter
      r0 = _mm_maddubs_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi16(r0, _mm_setzero_si128());

      // Save 4 resulting values to tmp. Last 3 values will be ignored later.
      _mm_storel_epi64((__m128i*) (tmpSample + outAddrOffset), r0);
    }
    else
      assert(y==0);
  }

  //for (int y = processing_height; y > 0; y -= 1) {
  //  printf("\nx %d ", y);
  //  for (int x = 0; x < dst_width; x++) { 
  //    int outOffset = processing_height - y;
  //    printf("%d ", *(buffer + x * (MAX_CU_SIZE/2 + 8) + outOffset));
  //  }
  //}

  // DEBUG. Dump the temp buffers to file
  //FILE *fp = fopen("temp_array.txt", "wb");
  //int nrBytesY = PicHeightInSamplesRefLayerY * PicWidthInSamplesCurLayerY;
  //int16_t *srcY = tmp;
  //fwrite(srcY, sizeof(int16_t), nrBytesY, fp);
  //fclose(fp);

  // Horizontal upsampling done. Perform vertical upsampling.
  // The temporary buffer (buffer) is flipped (x and y switched).
  // In this second filtering we will just perform the same flipping again.

  // 4. The variables shift1, shift2 and offset are derived as follows:
  //int shift2 = 12;     // (20 - BitDepthCurrC)  (H 34)
  //int offset = 2048;   // (1 << (shift2 - 1))   (H 35)
  __m128i offset = _mm_set_epi32(2048, 2048, 2048, 2048);
  // To clipt the output to 0...255
  __m128i clip8 = _mm_set_epi32(255, 255, 255, 255);

  int16_t *out;
  for (int y = 0; y < dst_height; y++) {
    yP = y_dst + y;

    // 1.
    // H.8.1.4.1.3 Derivation process for reference layer sample location in units of 1/16-th sample
    // The position_params array contains the precomputed values needed for this.
    yRef16 = (((yP - position_params[1]) * position_params[5] + position_params[7] + (1 << 11)) >> 12) + position_params[3];  // (H 64)

    // 3. The variables yRef and yPhase are derived as follows:
    yPhase = yRef16 & 15;;  // (H 32)
    yRef   = yRef16 >> 4;  // (H 31)

    // Get pointer to source and destination
    yRefBuf = yRef - y_src;
    tmpSample = buffer + yRefBuf;
    out = dst + y * dst_stride;

    // Load the correct filter accoding to xPhase
    filter = _mm_load_si128((__m128i *)(fC_16bit_2x + yPhase * 8));

    int x, outAddrOffset = 0;
    for (x = dst_width; x >= 8; x -= 8, outAddrOffset += 8) {
      // Load 8 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r3 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_unpacklo_epi64(r2, r3);

      r3 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r4 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r3 = _mm_unpacklo_epi64(r3, r4);

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r2 = _mm_add_epi32(r2, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r2 = _mm_srai_epi32(r2, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r2 = _mm_min_epi32(r2, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r2 = _mm_slli_epi32(r2, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 8 resulting 16 bit values to dst
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }

    // Process the remaining x lines (x lines remaining)
    if (x == 7) {
      // Load 7 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r3 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_unpacklo_epi64(r2, r3);

      r3 = _mm_loadl_epi64((__m128i *)tmpSample);

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);
      r3 = _mm_madd_epi16(r3, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, r3);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r2 = _mm_add_epi32(r2, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r2 = _mm_srai_epi32(r2, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r2 = _mm_min_epi32(r2, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r2 = _mm_slli_epi32(r2, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 7 resulting 16 bit values to dst. The last value should be 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 6) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r3 = _mm_loadl_epi64((__m128i *)tmpSample);
      r2 = _mm_unpacklo_epi64(r2, r3);

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r2 = _mm_add_epi32(r2, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r2 = _mm_srai_epi32(r2, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r2 = _mm_min_epi32(r2, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r2 = _mm_slli_epi32(r2, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 6 resulting 16 bit values to dst. The last 2 values should be 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 5) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_unpacklo_epi64(r1, r2);

      r2 = _mm_loadl_epi64((__m128i *)tmpSample); 

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      r2 = _mm_madd_epi16(r2, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      r2 = _mm_hadd_epi32(r2, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      r2 = _mm_add_epi32(r2, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);
      r2 = _mm_srai_epi32(r2, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      r2 = _mm_min_epi32(r2, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      r2 = _mm_slli_epi32(r2, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 5 resulting 16 bit values to dst. The last 3 values should be 0.
      _mm_store_si128((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 4) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r2 = _mm_loadl_epi64((__m128i *)tmpSample); 
      r1 = _mm_unpacklo_epi64(r1, r2);

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);
      
      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);
      
      // Add offset
      r0 = _mm_add_epi32(r0, offset);
      
      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);
      
      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);
      
      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 4 resulting 16 bit values to dst.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 3) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r0 = _mm_unpacklo_epi64(r0, r1);

      r1 = _mm_loadl_epi64((__m128i *)tmpSample); 

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);
      r1 = _mm_madd_epi16(r1, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, r1);

      // Add offset
      r0 = _mm_add_epi32(r0, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 3 resulting 16 bit values to dst. The last value should be 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 2) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);    tmpSample += (MAX_CU_SIZE / 2 + 8);
      r1 = _mm_loadl_epi64((__m128i *)tmpSample); 
      r0 = _mm_unpacklo_epi64(r0, r1);

      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 2 resulting 16 bit values to dst. The last 2 values should be 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else if (x == 1) {
      // Load 6 x lines
      r0 = _mm_loadl_epi64((__m128i *)tmpSample);
      
      // Multiply the 4 buffers by the double filter (output 4 32 bit values)
      r0 = _mm_madd_epi16(r0, filter);

      // Add up results
      r0 = _mm_hadd_epi32(r0, _mm_setzero_si128());

      // Add offset
      r0 = _mm_add_epi32(r0, offset);

      // Shift right by 12 bits
      r0 = _mm_srai_epi32(r0, 12);

      // Clip to 255
      r0 = _mm_min_epi32(r0, clip8);

      // And shift left again by 6 bits
      r0 = _mm_slli_epi32(r0, 6);

      r0 = _mm_hadd_epi16(r0, r2);  // The upper 16 bit of each result should be 0

      // Save 1 resulting 16 bit values to dst. The last 3 values should be 0.
      _mm_storel_epi64((__m128i*) (out + outAddrOffset), r0);
    }
    else {
      assert(x == 0);
    }
  }

  //// DEBUG. Dump the output to file
  //static FILE *fp = NULL;
  //if (fp == NULL)
  //  fp = fopen("upsampleC_out.txt", "wb");
  //fprintf(fp, "C Upsampling - w,h %d,%d x,y %d,%d\n", dst_width, dst_height, x_dst, y_dst);
  //for (int y = 0; y < dst_height; y++) {
  //  fprintf(fp, "%d | ", y);
  //  for (int x = 0; x < dst_width; x++) {
  //    fprintf(fp, "%d ", *(dst + x + y * dst_stride));
  //  }
  //  fprintf(fp, "\n");
  //}

}

#endif