/*
 * H.265 video codec.
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

#include "image.h"
#include "decctx.h"
#include "encoder/encoder-context.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <limits>


#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef HAVE_SSE4_1
#define MEMORY_PADDING  8
#else
#define MEMORY_PADDING  0
#endif

#define STANDARD_ALIGNMENT 16

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

static const int alignment = 16;

LIBDE265_API void* de265_alloc_image_plane(struct de265_image* img, int cIdx,
                                           void* inputdata, int inputstride, void *userdata)
{
  int alignment = STANDARD_ALIGNMENT;
  int stride = (img->get_width(cIdx) + alignment-1) / alignment * alignment;
  int height = img->get_height(cIdx);

  uint8_t* p = (uint8_t *)ALLOC_ALIGNED_16(stride * height + MEMORY_PADDING);

  if (p==NULL) { return NULL; }

  img->set_image_plane(cIdx, p, stride, userdata);

  // copy input data if provided

  if (inputdata != NULL) {
    if (inputstride == stride) {
      memcpy(p, inputdata, stride*height);
    }
    else {
      for (int y=0;y<height;y++) {
        memcpy(p+y*stride, ((char*)inputdata) + inputstride*y, inputstride);
      }
    }
  }

  return p;
}


LIBDE265_API void de265_free_image_plane(struct de265_image* img, int cIdx)
{
  uint8_t* p = (uint8_t*)img->get_image_plane(cIdx);
  assert(p);
  FREE_ALIGNED(p);
}


static int  de265_image_get_buffer(de265_decoder_context* ctx,
                                   de265_image_spec* spec, de265_image* img, void* userdata)
{
  const int rawChromaWidth  = spec->width  / img->SubWidthC;
  const int rawChromaHeight = spec->height / img->SubHeightC;

  int luma_stride   = (spec->width    + spec->alignment-1) / spec->alignment * spec->alignment;
  int chroma_stride = (rawChromaWidth + spec->alignment-1) / spec->alignment * spec->alignment;

  assert(img->BitDepth_Y >= 8 && img->BitDepth_Y <= 16);
  assert(img->BitDepth_C >= 8 && img->BitDepth_C <= 16);

  int luma_bpl   = luma_stride   * ((img->BitDepth_Y+7)/8);
  int chroma_bpl = chroma_stride * ((img->BitDepth_C+7)/8);

  int luma_height   = spec->height;
  int chroma_height = rawChromaHeight;

  bool alloc_failed = false;

  uint8_t* p[3] = { 0,0,0 };
  p[0] = (uint8_t *)ALLOC_ALIGNED_16(luma_height   * luma_bpl   + MEMORY_PADDING);
  if (p[0]==NULL) { alloc_failed=true; }

  if (img->get_chroma_format() != de265_chroma_mono) {
    p[1] = (uint8_t *)ALLOC_ALIGNED_16(chroma_height * chroma_bpl + MEMORY_PADDING);
    p[2] = (uint8_t *)ALLOC_ALIGNED_16(chroma_height * chroma_bpl + MEMORY_PADDING);

    if (p[1]==NULL || p[2]==NULL) { alloc_failed=true; }
  }
  else {
    p[1] = NULL;
    p[2] = NULL;
    chroma_stride = 0;
  }

  if (alloc_failed) {
    for (int i=0;i<3;i++)
      if (p[i]) {
        FREE_ALIGNED(p[i]);
      }

    return 0;
  }

  img->set_image_plane(0, p[0], luma_stride, NULL);
  img->set_image_plane(1, p[1], chroma_stride, NULL);
  img->set_image_plane(2, p[2], chroma_stride, NULL);
  
  // Allocate buffers for prediction/residual/trCoeff
  decoder_context* decCtx = (decoder_context*)ctx;
  for (int i = 0; i < 3; i++)
  {
    if ((i==0 && !decCtx->param_internals_save_prediction) ||
        (i==1 && !decCtx->param_internals_save_residual) ||
        (i==2 && !decCtx->param_internals_save_tr_coeff))
      continue;

    p[0] = (uint8_t *)ALLOC_ALIGNED_16(luma_height * luma_bpl + MEMORY_PADDING);
    if (p[0]==NULL) { alloc_failed=true; }

    if (img->get_chroma_format() != de265_chroma_mono) {
      p[1] = (uint8_t *)ALLOC_ALIGNED_16(chroma_height * chroma_bpl + MEMORY_PADDING);
      p[2] = (uint8_t *)ALLOC_ALIGNED_16(chroma_height * chroma_bpl + MEMORY_PADDING);

      if (p[1]==NULL || p[2]==NULL) { alloc_failed=true; }
    }
    else {
      p[1] = NULL;
      p[2] = NULL;
      chroma_stride = 0;
    }

    if (alloc_failed) {
      for (int i=0;i<3;i++)
        if (p[i]) {
          FREE_ALIGNED(p[i]);
        }

      return 0;
    }

    // Set the buffers
    for (int j = 0; j < 3; j++)
    {
      if (i == 0)
        img->set_image_plane_prediction(j, p[j]);
      else if (i == 1)
        img->set_image_plane_residual(j, p[j]);
      else if (i == 2)
        img->set_image_plane_tr_coeff(j, p[j]);
    }
  }

  return 1;
}

static void de265_image_release_buffer(de265_decoder_context* ctx,
                                       de265_image* img, void* userdata)
{
  for (int i=0;i<3;i++) {
    uint8_t* p = (uint8_t*)img->get_image_plane(i);
    if (p) {
      FREE_ALIGNED(p);
    }

    p = (uint8_t*)img->get_image_plane_prediction(i);
    if (p) {
      FREE_ALIGNED(p);
    }

    p = (uint8_t*)img->get_image_plane_residual(i);
    if (p) {
      FREE_ALIGNED(p);
    }

    p = (uint8_t*)img->get_image_plane_tr_coeff(i);
    if (p) {
      FREE_ALIGNED(p);
    }
  }
}


de265_image_allocation de265_image::default_image_allocation = {
  de265_image_get_buffer,
  de265_image_release_buffer
};


void de265_image::set_image_plane(int cIdx, uint8_t* mem, int stride, void *userdata)
{
  pixels[cIdx] = mem;
  plane_user_data[cIdx] = userdata;

  if (cIdx==0) { this->stride        = stride; }
  else         { this->chroma_stride = stride; }
}

void de265_image::set_image_plane_prediction(int cIdx, uint8_t* mem)
{
  pixels_prediction[cIdx] = mem;
}

void de265_image::set_image_plane_residual(int cIdx, uint8_t* mem)
{
  pixels_residual[cIdx] = mem;
}

void de265_image::set_image_plane_tr_coeff(int cIdx, uint8_t* mem)
{
  pixels_tr_coeff[cIdx] = mem;
}

uint32_t de265_image::s_next_image_ID = 0;

de265_image::de265_image()
{
  ID = -1;
  removed_at_picture_id = 0; // picture not used, so we can assume it has been removed

  decctx = NULL;
  encctx = NULL;

  encoder_image_release_func = NULL;

  //alloc_functions.get_buffer = NULL;
  //alloc_functions.release_buffer = NULL;

  for (int c=0;c<3;c++) {
    pixels[c] = NULL;
    pixels_confwin[c] = NULL;
    pixels_confwin_prediction[c] = NULL;
    pixels_confwin_residual[c] = NULL;
    pixels_confwin_tr_coeff[c] = NULL;
    plane_user_data[c] = NULL;

    pixels_prediction[c] = NULL;
    pixels_residual[c] = NULL;
    pixels_tr_coeff[c] = NULL;
  }

  width=height=0;

  pts = 0;
  user_data = NULL;

  ctb_progress = NULL;

  integrity = INTEGRITY_NOT_DECODED;

  picture_order_cnt_lsb = -1; // undefined
  PicOrderCntVal = -1; // undefined
  PicState = UnusedForReference;
  PicOutputFlag = false;

  nThreadsQueued   = 0;
  nThreadsRunning  = 0;
  nThreadsBlocked  = 0;
  nThreadsFinished = 0;
  nThreadsTotal    = 0;

  de265_mutex_init(&mutex);
  de265_cond_init(&finished_cond);

  bIlRefPic = false;
  equalPictureSizeAndOffsetFlag = false;
  interLayerMotionPredictionEnabled = false;
  for (int i = 0; i<10; i++) {
    il_scaling_parameters[i] = -1;
  }
}


de265_error de265_image::alloc_image(int w,int h, enum de265_chroma c,
                                     std::shared_ptr<const seq_parameter_set> sps, bool allocMetadata,
                                     decoder_context* dctx,
                                     encoder_context* ectx,
                                     de265_PTS pts, void* user_data,
                                     bool useCustomAllocFunc,
                                     bool interLayerReferencePicture)
{
  //if (allocMetadata) { assert(sps); }
  if (allocMetadata) { assert(sps); }

  if (sps) { this->sps = sps; }

  release(); /* TODO: review code for efficient allocation when arrays are already
                allocated to the requested size. Without the release, the old image-data
                will not be freed. */

  ID = s_next_image_ID++;
  removed_at_picture_id = std::numeric_limits<int32_t>::max();

  decctx = dctx;
  encctx = ectx;

  // --- allocate image buffer ---

  chroma_format= c;
  bIlRefPic = interLayerReferencePicture;

  width = w;
  height = h;
  chroma_width = w;
  chroma_height= h;

  this->user_data = user_data;
  this->pts = pts;

  de265_image_spec spec;

  int WinUnitX, WinUnitY;

  switch (chroma_format) {
  case de265_chroma_mono: WinUnitX=1; WinUnitY=1; break;
  case de265_chroma_420:  WinUnitX=2; WinUnitY=2; break;
  case de265_chroma_422:  WinUnitX=2; WinUnitY=1; break;
  case de265_chroma_444:  WinUnitX=1; WinUnitY=1; break;
  default:
    assert(0);
  }

  switch (chroma_format) {
  case de265_chroma_420:
    spec.format = de265_image_format_YUV420P8;
    chroma_width  = (chroma_width +1)/2;
    chroma_height = (chroma_height+1)/2;
    SubWidthC  = 2;
    SubHeightC = 2;
    break;

  case de265_chroma_422:
    spec.format = de265_image_format_YUV422P8;
    chroma_width = (chroma_width+1)/2;
    SubWidthC  = 2;
    SubHeightC = 1;
    break;

  case de265_chroma_444:
    spec.format = de265_image_format_YUV444P8;
    SubWidthC  = 1;
    SubHeightC = 1;
    break;

  case de265_chroma_mono:
    spec.format = de265_image_format_mono8;
    chroma_width = 0;
    chroma_height= 0;
    SubWidthC  = 1;
    SubHeightC = 1;
    break;

  default:
    assert(false);
    break;
  }

  if (sps) {
    assert(sps->SubWidthC  == SubWidthC);
    assert(sps->SubHeightC == SubHeightC);
  }

  spec.width  = w;
  spec.height = h;
  spec.alignment = STANDARD_ALIGNMENT;


  // conformance window cropping

  int left   = sps ? sps->conf_win_left_offset : 0;
  int right  = sps ? sps->conf_win_right_offset : 0;
  int top    = sps ? sps->conf_win_top_offset : 0;
  int bottom = sps ? sps->conf_win_bottom_offset : 0;

  width_confwin = width - (left+right)*WinUnitX;
  height_confwin= height- (top+bottom)*WinUnitY;
  chroma_width_confwin = chroma_width -left-right;
  chroma_height_confwin= chroma_height-top-bottom;

  spec.crop_left  = left *WinUnitX;
  spec.crop_right = right*WinUnitX;
  spec.crop_top   = top   *WinUnitY;
  spec.crop_bottom= bottom*WinUnitY;

  spec.visible_width = width_confwin;
  spec.visible_height= height_confwin;


  BitDepth_Y = (sps==NULL) ? 8 : sps->BitDepth_Y;
  BitDepth_C = (sps==NULL) ? 8 : sps->BitDepth_C;

  bpp_shift[0] = (BitDepth_Y <= 8) ? 0 : 1;
  bpp_shift[1] = (BitDepth_C <= 8) ? 0 : 1;
  bpp_shift[2] = bpp_shift[1];

  if (interLayerReferencePicture) {
    // Do not allocate an actual picture or metadata.
    return DE265_OK;
  }

  // allocate memory and set conformance window pointers

  void* alloc_userdata = NULL;
  if (decctx) alloc_userdata = decctx->param_image_allocation_userdata;
  if (encctx) alloc_userdata = encctx->param_image_allocation_userdata; // actually not needed

  if (encctx && useCustomAllocFunc) {
    encoder_image_release_func = encctx->release_func;

    // if we do not provide a release function, use our own

    if (encoder_image_release_func == NULL) {
      image_allocation_functions = de265_image::default_image_allocation;
    }
    else {
      image_allocation_functions.get_buffer     = NULL;
      image_allocation_functions.release_buffer = NULL;
    }
  }
  else if (decctx && useCustomAllocFunc) {
    image_allocation_functions = decctx->param_image_allocation_functions;
  }
  else {
    image_allocation_functions = de265_image::default_image_allocation;
  }

  bool mem_alloc_success = true;

  if (image_allocation_functions.get_buffer != NULL) {
    mem_alloc_success = image_allocation_functions.get_buffer(decctx, &spec, this,
                                                              alloc_userdata);

    pixels_confwin[0] = pixels[0] + left*WinUnitX + top*WinUnitY*stride;
    pixels_confwin[1] = pixels[1] + left + top*chroma_stride;
    pixels_confwin[2] = pixels[2] + left + top*chroma_stride;

    if (pixels_prediction[0])
    {
      pixels_confwin_prediction[0] = pixels_prediction[0] + left*WinUnitX + top*WinUnitY*stride;
      pixels_confwin_prediction[1] = pixels_prediction[1] + left + top*chroma_stride;
      pixels_confwin_prediction[2] = pixels_prediction[2] + left + top*chroma_stride;
    }
    if (pixels_residual[0])
    {
      pixels_confwin_residual[0] = pixels_residual[0] + left*WinUnitX + top*WinUnitY*stride;
      pixels_confwin_residual[1] = pixels_residual[1] + left + top*chroma_stride;
      pixels_confwin_residual[2] = pixels_residual[2] + left + top*chroma_stride;
    }
    if (pixels_tr_coeff[0])
    {
      pixels_confwin_tr_coeff[0] = pixels_tr_coeff[0] + left*WinUnitX + top*WinUnitY*stride;
      pixels_confwin_tr_coeff[1] = pixels_tr_coeff[1] + left + top*chroma_stride;
      pixels_confwin_tr_coeff[2] = pixels_tr_coeff[2] + left + top*chroma_stride;
    }

    // check for memory shortage

    if (!mem_alloc_success)
      {
        return DE265_ERROR_OUT_OF_MEMORY;
      }
  }

  //alloc_functions = *allocfunc;
  //alloc_userdata  = userdata;

  // --- allocate decoding info arrays ---

  if (allocMetadata) {
    // For inter layer reference pictures only allocate cb_info and pb_info
      
    // intra pred mode
    mem_alloc_success &= intraPredMode.alloc(sps->PicWidthInMinPUs, sps->PicHeightInMinPUs,
                                              sps->Log2MinPUSize);

    mem_alloc_success &= intraPredModeC.alloc(sps->PicWidthInMinPUs, sps->PicHeightInMinPUs,
                                              sps->Log2MinPUSize);

    // tu info
    mem_alloc_success &= tu_info.alloc(sps->PicWidthInTbsY, sps->PicHeightInTbsY,
                                        sps->Log2MinTrafoSize);

    // deblk info
    int deblk_w = (sps->pic_width_in_luma_samples +3)/4;
    int deblk_h = (sps->pic_height_in_luma_samples+3)/4;
    mem_alloc_success &= deblk_info.alloc(deblk_w, deblk_h, 2);
    
    // cb info
    mem_alloc_success &= cb_info.alloc(sps->PicWidthInMinCbsY, sps->PicHeightInMinCbsY,
                                       sps->Log2MinCbSizeY);

    // pb info
    int puWidth  = sps->PicWidthInMinCbsY  << (sps->Log2MinCbSizeY -2);
    int puHeight = sps->PicHeightInMinCbsY << (sps->Log2MinCbSizeY -2);
    mem_alloc_success &= pb_info.alloc(puWidth,puHeight, 2);

    // CTB info
    if (ctb_info.data_size != sps->PicSizeInCtbsY) {
      delete[] ctb_progress;

      mem_alloc_success &= ctb_info.alloc(sps->PicWidthInCtbsY, sps->PicHeightInCtbsY,
                                          sps->Log2CtbSizeY);

      ctb_progress = new de265_progress_lock[ ctb_info.data_size ];
    }

    // check for memory shortage
    if (!mem_alloc_success) {
      return DE265_ERROR_OUT_OF_MEMORY;
    }
  }

  return DE265_OK;
}


de265_image::~de265_image()
{
  release();

  // free progress locks

  if (ctb_progress) {
    delete[] ctb_progress;
  }

  de265_cond_destroy(&finished_cond);
  de265_mutex_destroy(&mutex);
}


void de265_image::release()
{
  // free image memory

  if (pixels[0])
    {
      if (encoder_image_release_func != NULL) {
        encoder_image_release_func(encctx, this,
                                   encctx->param_image_allocation_userdata);
      }
      else {
        image_allocation_functions.release_buffer(decctx, this,
                                                decctx ?
                                                  decctx->param_image_allocation_userdata :
                                                  NULL);
      }

      for (int i=0;i<3;i++)
        {
          pixels[i] = NULL;
          pixels_confwin[i] = NULL;
          pixels_confwin_prediction[i] = NULL;
          pixels_confwin_residual[i] = NULL;
          pixels_confwin_tr_coeff[i] = NULL;
        }
    }

  // free slices

  for (int i=0;i<slices.size();i++) {
    delete slices[i];
  }
  slices.clear();
}


void de265_image::fill_image(int y,int cb,int cr)
{
  if (y>=0) {
    memset(pixels[0], y, stride * height);
  }

  if (cb>=0) {
    memset(pixels[1], cb, chroma_stride * chroma_height);
  }

  if (cr>=0) {
    memset(pixels[2], cr, chroma_stride * chroma_height);
  }
}


de265_error de265_image::copy_image(const de265_image* src)
{
  /* TODO: actually, since we allocate the image only for internal purpose, we
     do not have to call the external allocation routines for this. However, then
     we have to track for each image how to release it again.
     Another option would be to safe the copied data not in an de265_image at all.
  */

  de265_error err = alloc_image(src->width, src->height, src->chroma_format, src->sps, false,
                                src->decctx, src->encctx, src->pts, src->user_data, false);
  if (err != DE265_OK) {
    return err;
  }

  copy_lines_from(src, 0, src->height);

  return err;
}

void de265_image::set_lower_layer_picture(const de265_image* src)
{
  // Save pointer to lower layer reference
  assert( bIlRefPic );
  ilRefPic = src;
  bIlRefPic = true;

  // Copy the pointers to the slice segment headers.
  // TODO: Is this a good idea?
  slices.clear();
  for (int i = 0; i < src->slices.size(); i++) {
    slices.push_back(src->slices.at(i));
  }
}

void de265_image::set_inter_layer_metadata_scaling_parameters(int scaling_parameters[10])
{
  for (int i = 0; i < 10; i++) {
    il_scaling_parameters[i] = scaling_parameters[i];
  }
}

void de265_image::set_inter_layer_upsampling_parameters(int upsampling_params[2][10])
{
  for (int j = 0; j < 2; j++) {
    for (int i = 0; i < 10; i++) {
      il_upsampling_parameters[j][i] = upsampling_params[j][i];
    }
  }
}

// end = last line + 1
void de265_image::copy_lines_from(const de265_image* src, int first, int end)
{
  if (end > src->height) end=src->height;

  assert(first % 2 == 0);
  assert(end   % 2 == 0);

  int luma_bpp   = (sps->BitDepth_Y+7)/8;
  int chroma_bpp = (sps->BitDepth_C+7)/8;

  if (src->stride == stride) {
    memcpy(pixels[0]      + first*stride * luma_bpp,
           src->pixels[0] + first*src->stride * luma_bpp,
           (end-first)*stride * luma_bpp);
  }
  else {
    for (int yp=first;yp<end;yp++) {
      memcpy(pixels[0]+yp*stride * luma_bpp,
             src->pixels[0]+yp*src->stride * luma_bpp,
             src->width * luma_bpp);
    }
  }

  int first_chroma = first / src->SubHeightC;
  int end_chroma   = end   / src->SubHeightC;

  if (src->chroma_format != de265_chroma_mono) {
    if (src->chroma_stride == chroma_stride) {
      memcpy(pixels[1]      + first_chroma*chroma_stride * chroma_bpp,
             src->pixels[1] + first_chroma*chroma_stride * chroma_bpp,
             (end_chroma-first_chroma) * chroma_stride * chroma_bpp);
      memcpy(pixels[2]      + first_chroma*chroma_stride * chroma_bpp,
             src->pixels[2] + first_chroma*chroma_stride * chroma_bpp,
             (end_chroma-first_chroma) * chroma_stride * chroma_bpp);
    }
    else {
      for (int y=first_chroma;y<end_chroma;y++) {
        memcpy(pixels[1]+y*chroma_stride * chroma_bpp,
               src->pixels[1]+y*src->chroma_stride * chroma_bpp,
               src->chroma_width * chroma_bpp);
        memcpy(pixels[2]+y*chroma_stride * chroma_bpp,
               src->pixels[2]+y*src->chroma_stride * chroma_bpp,
               src->chroma_width * chroma_bpp);
      }
    }
  }
}

void de265_image::get_pixel_pointers_from(de265_image *src)
{
  if (width != src->get_width() || height != src->get_height()) {
    assert( false );
  }
  assert( bIlRefPic );

  // Pixel data
  pixels[0] = src->pixels[0];
  pixels[1] = src->pixels[1];
  pixels[2] = src->pixels[2];
  // Conformance windows pixels data
  pixels_confwin[0] = src->pixels_confwin[0];
  pixels_confwin[1] = src->pixels_confwin[1];
  pixels_confwin[2] = src->pixels_confwin[2];
  // Strides
  stride        = src->stride;
  chroma_stride = src->chroma_stride;
}

//void de265_image::upsample_image_from(decoder_context* ctx, de265_image* rlPic, int upsampling_params[2][10])
//{
//  assert( bIlRefPic );
//
//  int src_size[2] = {rlPic->get_width(), rlPic->get_height()};
//  int dst_size[2] = { get_width(), get_height() };
//  ctx->acceleration.resampling_process_of_luma_sample_values(rlPic->get_image_plane(0), rlPic->get_luma_stride(), src_size,
//                                                              get_image_plane(0), get_luma_stride(), dst_size,
//                                                              upsampling_params[0] );
//
//  // Chroma
//  src_size[0] = rlPic->get_width(1);
//  src_size[1] = rlPic->get_height(1);
//  dst_size[0] = get_width(1);
//  dst_size[1] = get_height(1);
//  ctx->acceleration.resampling_process_of_chroma_sample_values(rlPic->get_image_plane(1), rlPic->get_chroma_stride(), src_size,
//                                                              get_image_plane(1), get_chroma_stride(), dst_size,
//                                                              upsampling_params[1] );
//  ctx->acceleration.resampling_process_of_chroma_sample_values(rlPic->get_image_plane(2), rlPic->get_chroma_stride(), src_size,
//                                                              get_image_plane(2), get_chroma_stride(), dst_size,
//                                                              upsampling_params[1] );
//}

//void de265_image::colour_mapping(decoder_context* ctx, de265_image* rlPic, colour_mapping_table *map, int colourMappingParams[2])
//{
//  assert( bIlRefPic );
//
//  int PicWidthInSamplesRefLayerY  = rlPic->get_width();
//  int PicHeightInSamplesRefLayerY = rlPic->get_height();
//  int PicHeightInSamplesRefLayerC = rlPic->get_height(1);
//  int PicWidthInSamplesRefLayerC  = rlPic->get_width(1);
//
//  int SubWidthRefLayerC = colourMappingParams[0];
//  int SubHeightRefLayerC = colourMappingParams[1];
//
//  int maxValOut = (1 << map->BitDepthCmOutputY) - 1;
//
//  uint8_t* src_Y  = rlPic->get_image_plane(0);
//  uint8_t* src_cb = rlPic->get_image_plane(1);
//  uint8_t* src_cr = rlPic->get_image_plane(2);
//  int strideY = rlPic->get_luma_stride();
//  int strideC = rlPic->get_chroma_stride();
//
//  uint8_t* dst_Y  = get_image_plane(0);
//  
//  for (int xP = 0; xP < PicWidthInSamplesRefLayerY - 1; xP++) {
//    for (int yP = 0; yP < PicHeightInSamplesRefLayerY - 1; yP++) {
//      
//      // The chroma sample location ( xPC, yPC ) is set equal to ( xP / SubWidthRefLayerC, yP / SubHeightRefLayerC ).
//      int xPC = xP / SubWidthRefLayerC;
//      int yPC = yP / SubHeightRefLayerC;
//
//      // 1. The value of cmLumaSample is derived by applying the following ordered steps:
//      int yShift2Idx = map->BitDepthCmInputY - map->cm_octant_depth - map->cm_y_part_num_log2;  // (H 80)
//      int cShift2Idx = map->BitDepthCmInputC - map->cm_octant_depth;                            // (H 81)
//
//      // 2. The variables nMappingShift and nMappingOffset are derived as follows:
//      int nMappingShift = 10 + map->BitDepthCmInputY - map->BitDepthCmOutputY;  // (H 82)
//      int nMappingOffset = 1 << ( nMappingShift - 1 );                          // (H 83)
//
//      // 3. The variables tempCb and tempCr are derived as follows:
//      int tempCb, tempCr;
//
//      // Get pointers to y line for yPC
//      uint8_t* rlPicSampleCb_yPC = src_cb + yPC * strideC;
//      uint8_t* rlPicSampleCr_yPC = src_cr + yPC * strideC;
//      
//      if (SubWidthRefLayerC == 2 && SubHeightRefLayerC == 2) {
//        if ((xP % 2) == 0 && (yP % 2) == 0) {
//          int yP2C = libde265_max( 0, yPC - 1 );                                      // (H 84)
//          
//          // Get poitner to y line for yPC2
//          uint8_t* rlPicSampleCb_yP2C = src_cb + yP2C * strideC;
//          uint8_t* rlPicSampleCr_yP2C = src_cr + yP2C * strideC;
//
//          tempCb = (rlPicSampleCb_yPC[xPC] * 3 + rlPicSampleCb_yP2C[xPC] + 2) >> 2;   // (H 85)
//          tempCr = (rlPicSampleCr_yP2C[xPC] * 3 + rlPicSampleCr_yP2C[xPC] + 2) >> 2;  // (H 86)
//        }
//        else if ((xP % 2) == 0 && (yP % 2) == 1) {
//          int yP2C = libde265_min( yPC + 1, PicHeightInSamplesRefLayerC - 1 );        // (H 87)
//
//          // Get poitner to y line for yPC2
//          uint8_t* rlPicSampleCb_yP2C = src_cb + yP2C * strideC;
//          uint8_t* rlPicSampleCr_yP2C = src_cr + yP2C * strideC;
//
//          tempCb = (rlPicSampleCb_yPC[xPC] * 3 + rlPicSampleCb_yP2C[xPC] + 2 ) >> 2;   // (H 88)
//          tempCr = (rlPicSampleCr_yPC[xPC] * 3 + rlPicSampleCr_yP2C[xPC] + 2 ) >> 2;   // (H 89)
//        }
//        else if ((xP % 2) == 1 && (yP % 2) == 0) {
//          int xP2C = libde265_min( xPC + 1, PicWidthInSamplesRefLayerC - 1 );          // (H 90)
//          int yP2C = libde265_max( 0, yPC - 1 );                                       // (H 91)
//
//          // Get poitner to y line for yPC2
//          uint8_t* rlPicSampleCb_yP2C = src_cb + yP2C * strideC;
//          uint8_t* rlPicSampleCr_yP2C = src_cr + yP2C * strideC;
//
//          tempCb = (rlPicSampleCb_yP2C[xPC] + rlPicSampleCb_yP2C[xP2C] + (rlPicSampleCb_yPC[xPC] + rlPicSampleCb_yPC[xP2C]) * 3 + 4) >> 3;  // (H 92)
//          tempCr = (rlPicSampleCr_yP2C[xPC] + rlPicSampleCr_yP2C[xP2C] + (rlPicSampleCr_yPC[yPC] + rlPicSampleCr_yPC[xP2C]) * 3 + 4) >> 3;  // (H 93)
//        }
//        else if ((xP % 2) == 1 && (yP % 2) == 1) {
//          int xP2C = libde265_min( xPC + 1, PicWidthInSamplesRefLayerC - 1 );          // (H 94)
//          int yP2C = libde265_min( yPC + 1, PicHeightInSamplesRefLayerC - 1 );         // (H 95)
//
//          // Get poitner to y line for yPC2
//          uint8_t* rlPicSampleCb_yP2C = src_cb + yP2C * strideC;
//          uint8_t* rlPicSampleCr_yP2C = src_cr + yP2C * strideC;
//
//          tempCb = ((rlPicSampleCb_yPC[xPC] + rlPicSampleCb_yPC[xP2C]) * 3 + rlPicSampleCb_yP2C[xPC] + rlPicSampleCb_yP2C[xP2C] + 4) >> 3;  // (H 96)
//          tempCr = ((rlPicSampleCr_yPC[xPC] + rlPicSampleCr_yPC[xP2C]) * 3 + rlPicSampleCr_yP2C[xPC] + rlPicSampleCr_yP2C[xP2C] + 4) >> 3;  // (H 97)
//        }
//      }
//      else if (SubWidthRefLayerC == 2) {
//        if ((xP % 2) == 1) {
//          int xP2C = libde265_min( xPC + 1, PicWidthInSamplesRefLayerC - 1 );          // (H 98)
//          
//          tempCb = (rlPicSampleCb_yPC[xPC] + rlPicSampleCb_yPC[xP2C] + 1 ) >> 1;       // (H 99)
//          tempCr = (rlPicSampleCr_yPC[xPC] + rlPicSampleCr_yPC[xP2C] + 1 ) >> 1;       // (H 100)
//        }
//        else {
//          tempCb = rlPicSampleCb_yPC[xPC];                                             // (H 101)
//          tempCr = rlPicSampleCr_yPC[xPC];                                             // (H 102)
//        }
//      }
//      else {
//        tempCb = rlPicSampleCb_yPC[xPC];                                               // (H 103)
//        tempCr = rlPicSampleCr_yPC[xPC];                                               // (H 104)
//      }
//
//      // Get pointer to input/output yP line
//      uint8_t* cmLumaSample_yP = dst_Y + yP * stride;
//      uint8_t* rlPicSampleY_yP = src_Y + yP * strideY;
//
//      // 4. The value of cmLumaSample is derived as follows:
//      int idxY = rlPicSampleY_yP[xP] >> yShift2Idx;
//      int idxCb = (map->cm_octant_depth == 1) ? (tempCb >= map->CMThreshU) : (tempCb >> cShift2Idx);
//      int idxCr = (map->cm_octant_depth == 1) ? (tempCr >= map->CMThreshV) : (tempCr >> cShift2Idx);
//
//      cmLumaSample_yP[xP] = Clip3( 0, maxValOut,
//                           ((map->LutY[idxY][idxCb][idxCr][0] * rlPicSampleY_yP[xP] + map->LutY[idxY][idxCb][idxCr][1] * tempCb + 
//                             map->LutY[idxY][idxCb][idxCr][2] * tempCr + nMappingOffset ) >> nMappingShift ) + 
//                             map->LutY[idxY][idxCb][idxCr][3]); // (H 108) (H-109)
//    }
//  }
//}

void de265_image::exchange_pixel_data_with(de265_image& b)
{
  for (int i=0;i<3;i++) {
    std::swap(pixels[i], b.pixels[i]);
    std::swap(pixels_confwin[i], b.pixels_confwin[i]);
    std::swap(plane_user_data[i], b.plane_user_data[i]);
  }

  std::swap(stride, b.stride);
  std::swap(chroma_stride, b.chroma_stride);
  std::swap(image_allocation_functions, b.image_allocation_functions);
}


void de265_image::thread_start(int nThreads)
{
  de265_mutex_lock(&mutex);

  //printf("nThreads before: %d %d\n",nThreadsQueued, nThreadsTotal);

  nThreadsQueued += nThreads;
  nThreadsTotal += nThreads;

  //printf("nThreads after: %d %d\n",nThreadsQueued, nThreadsTotal);

  de265_mutex_unlock(&mutex);
}

void de265_image::thread_run(const thread_task* task)
{
  //printf("run thread %s\n", task->name().c_str());

  de265_mutex_lock(&mutex);
  nThreadsQueued--;
  nThreadsRunning++;
  de265_mutex_unlock(&mutex);
}

void de265_image::thread_blocks()
{
  de265_mutex_lock(&mutex);
  nThreadsRunning--;
  nThreadsBlocked++;
  de265_mutex_unlock(&mutex);
}

void de265_image::thread_unblocks()
{
  de265_mutex_lock(&mutex);
  nThreadsBlocked--;
  nThreadsRunning++;
  de265_mutex_unlock(&mutex);
}

void de265_image::thread_finishes(const thread_task* task)
{
  //printf("finish thread %s\n", task->name().c_str());

  de265_mutex_lock(&mutex);

  nThreadsRunning--;
  nThreadsFinished++;
  assert(nThreadsRunning >= 0);

  if (nThreadsFinished==nThreadsTotal) {
    de265_cond_broadcast(&finished_cond, &mutex);
  }

  de265_mutex_unlock(&mutex);
}

void de265_image::wait_for_progress(thread_task* task, int ctbx,int ctby, int progress)
{
  const int ctbW = sps->PicWidthInCtbsY;

  wait_for_progress(task, ctbx + ctbW*ctby, progress);
}

void de265_image::wait_for_progress(thread_task* task, int ctbAddrRS, int progress)
{
  if (task==NULL) { return; }

  de265_progress_lock* progresslock = &ctb_progress[ctbAddrRS];
  if (progresslock->get_progress() < progress) {
    thread_blocks();

    assert(task!=NULL);
    task->state = thread_task::Blocked;

    /* TODO: check whether we are the first blocked task in the list.
       If we are, we have to conceal input errors.
       Simplest concealment: do not block.
    */

    progresslock->wait_for_progress(progress);
    task->state = thread_task::Running;
    thread_unblocks();
  }
}


void de265_image::wait_for_completion()
{
  de265_mutex_lock(&mutex);
  while (nThreadsFinished!=nThreadsTotal) {
    de265_cond_wait(&finished_cond, &mutex);
  }
  de265_mutex_unlock(&mutex);
}

bool de265_image::debug_is_completed() const
{
  return nThreadsFinished==nThreadsTotal;
}



void de265_image::clear_metadata()
{
  // TODO: maybe we could avoid the memset by ensuring that all data is written to
  // during decoding (especially log2CbSize), but it is unlikely to be faster than the memset.

  cb_info.clear();
  //tu_info.clear();  // done on the fly
  ctb_info.clear();
  deblk_info.clear();

  // --- reset CTB progresses ---

  for (int i=0;i<ctb_info.data_size;i++) {
    ctb_progress[i].reset(CTB_PROGRESS_NONE);
  }
}


void de265_image::set_mv_info(int x,int y, int nPbW,int nPbH, const PBMotion& mv)
{
  assert(!bIlRefPic);
  int log2PuSize = 2;

  int xPu = x >> log2PuSize;
  int yPu = y >> log2PuSize;
  int wPu = nPbW >> log2PuSize;
  int hPu = nPbH >> log2PuSize;

  int stride = pb_info.width_in_units;

  for (int pby=0;pby<hPu;pby++)
    for (int pbx=0;pbx<wPu;pbx++)
      {
        pb_info[ xPu+pbx + (yPu+pby)*stride ] = mv;
      }
}

void de265_image::set_pred_mode(int x, int y, int nPbW, int nPbH, enum PredMode mode)
{
  int log2PuSize = 2;

  int cbX = x >> cb_info.log2unitSize;
  int cbY = y >> cb_info.log2unitSize;

  int wPu = nPbW >> cb_info.log2unitSize;
  int hPu = nPbH >> cb_info.log2unitSize;

  for (int cby = cbY; cby<cbY + hPu; cby++) {
    for (int cbx=cbX;cbx<cbX+wPu;cbx++) {
      assert( cbx + cby*cb_info.width_in_units < cb_info.size() );
      cb_info[ cbx + cby*cb_info.width_in_units ].PredMode = mode;
    }
  }
}


bool de265_image::available_zscan(int xCurr,int yCurr, int xN,int yN) const
{
  if (xN<0 || yN<0) return false;
  if (xN>=sps->pic_width_in_luma_samples ||
      yN>=sps->pic_height_in_luma_samples) return false;

  int minBlockAddrN = pps->MinTbAddrZS[ (xN>>sps->Log2MinTrafoSize) +
                                        (yN>>sps->Log2MinTrafoSize) * sps->PicWidthInTbsY ];
  int minBlockAddrCurr = pps->MinTbAddrZS[ (xCurr>>sps->Log2MinTrafoSize) +
                                           (yCurr>>sps->Log2MinTrafoSize) * sps->PicWidthInTbsY ];

  if (minBlockAddrN > minBlockAddrCurr) return false;

  int xCurrCtb = xCurr >> sps->Log2CtbSizeY;
  int yCurrCtb = yCurr >> sps->Log2CtbSizeY;
  int xNCtb = xN >> sps->Log2CtbSizeY;
  int yNCtb = yN >> sps->Log2CtbSizeY;

  if (get_SliceAddrRS(xCurrCtb,yCurrCtb) !=
      get_SliceAddrRS(xNCtb,   yNCtb)) {
    return false;
  }

  if (pps->TileIdRS[xCurrCtb + yCurrCtb*sps->PicWidthInCtbsY] !=
      pps->TileIdRS[xNCtb    + yNCtb   *sps->PicWidthInCtbsY]) {
    return false;
  }

  return true;
}


bool de265_image::available_pred_blk(int xC,int yC, int nCbS, int xP, int yP,
                                     int nPbW, int nPbH, int partIdx, int xN,int yN) const
{
  logtrace(LogMotion,"C:%d;%d P:%d;%d N:%d;%d size=%d;%d\n",xC,yC,xP,yP,xN,yN,nPbW,nPbH);

  int sameCb = (xC <= xN && xN < xC+nCbS &&
                yC <= yN && yN < yC+nCbS);

  bool availableN;

  if (!sameCb) {
    availableN = available_zscan(xP,yP,xN,yN);
  }
  else {
    availableN = !(nPbW<<1 == nCbS && nPbH<<1 == nCbS &&  // NxN
                   partIdx==1 &&
                   yN >= yC+nPbH && xN < xC+nPbW);  // xN/yN inside partIdx 2
  }

  if (availableN && get_pred_mode(xN,yN) == MODE_INTRA) {
    availableN = false;
  }

  return availableN;
}

const PBMotion& de265_image::get_mv_info_lower_layer(int x, int y) const
{
  assert(ilRefPic != NULL);

  if (equalPictureSizeAndOffsetFlag) {
    // No scaling necessary
    return ilRefPic->get_mv_info(x, y);
  }

  // Invoke the upsampling process for motion vectors 

  int PicWidthInSamplesRefLayerY  = ilRefPic->get_width();
  int PicHeightInSamplesRefLayerY = ilRefPic->get_height();

  int ScaledRefRegionWidthInSamplesY  = il_scaling_parameters[6];
  int RefLayerRegionWidthInSamplesY   = il_scaling_parameters[7];
  int ScaledRefRegionHeightInSamplesY = il_scaling_parameters[8];
  int RefLayerRegionHeightInSamplesY  = il_scaling_parameters[9];

  // 1. The center location ( xPCtr, yPCtr ) of the luma prediction block is derived as follows:
  int xPCtr = x + 8;  // (H 65)
  int yPCtr = y + 8;  // (H 66)

  // 2. The variables xRef and yRef are derived as follows:
  int xRef = (((xPCtr - il_scaling_parameters[0]) * il_scaling_parameters[2] + (1 << 15)) >> 16 ) + il_scaling_parameters[4];  // (H 67)
  int yRef = (((yPCtr - il_scaling_parameters[1]) * il_scaling_parameters[3] + (1 << 15)) >> 16 ) + il_scaling_parameters[5];  // (H 68)

  // 3. The rounded reference layer luma sample location ( xRL, yRL ) is derived as follows:
  int xRL = ((xRef + 4) >> 4) << 4;  // (H 69)
  int yRL = ((yRef + 4) >> 4) << 4;  // (H 70)

  // 4. Upsample the prediction mode. (H 71)
  PredMode rsPredMode;
  if( xRL < 0 || xRL >= PicWidthInSamplesRefLayerY || yRL < 0 || yRL >= PicHeightInSamplesRefLayerY ) {
    rsPredMode = MODE_INTRA;
  }
  else {
    rsPredMode = ilRefPic->get_pred_mode(xRL, yRL);
  }

  // 5. Upsample the motion vectors and prediction flags
  PBMotion mv_dst;
  if (rsPredMode != MODE_INTRA) {
    const PBMotion mv_src = ilRefPic->get_mv_info(xRL, yRL);
    // For X being each of 0 and 1...
    for (int l=0; l<2; l++) {
      // RefIdx, predFlag
      mv_dst.refIdx[l] = mv_src.refIdx[l];     // (H 72)
      mv_dst.predFlag[l] = mv_src.predFlag[l]; // (H 73)

      // Motion vector. X-component.
      if (ScaledRefRegionWidthInSamplesY != RefLayerRegionWidthInSamplesY) {
        int rlMvLX = mv_src.mv[l].x;
        int scaleMVX = Clip3( -4096, 4095, ((ScaledRefRegionWidthInSamplesY << 8) + (RefLayerRegionWidthInSamplesY >> 1)) / RefLayerRegionWidthInSamplesY); // (H 74)
        mv_dst.mv[l].x = Clip3( -32768, 32767, Sign( scaleMVX *	rlMvLX) * ((abs_value(scaleMVX * rlMvLX) + 127) >> 8)); // (H 75)
      }
      else {
        mv_dst.mv[l].x = mv_src.mv[l].x; // (H 76)
      }

      // Motion vector. Y-component.
      if (ScaledRefRegionHeightInSamplesY != RefLayerRegionHeightInSamplesY) {
        int rlMvLX = mv_src.mv[l].y;
        int scaleMVX = Clip3( -4096, 4095, ((ScaledRefRegionHeightInSamplesY << 8) + (RefLayerRegionHeightInSamplesY >> 1)) / RefLayerRegionHeightInSamplesY); // (H 77)
        mv_dst.mv[l].y = Clip3( -32768, 32767, Sign( scaleMVX *	rlMvLX) * ((abs_value(scaleMVX * rlMvLX) + 127) >> 8)); // (H 78)
      }
      else {
        mv_dst.mv[l].y = mv_src.mv[l].y;  // (H 79)
      }
    }
  }
  else {
    // Otherwise (rsPredMode is equal to MODE_INTRA), the following applies:
    mv_dst.mv[0].x = 0;
    mv_dst.mv[0].y = 0;
    mv_dst.mv[1].x = 0;
    mv_dst.mv[1].y = 0;
    mv_dst.refIdx[0] = -1;
    mv_dst.refIdx[1] = -1;
    mv_dst.predFlag[0] = 0;
    mv_dst.predFlag[1] = 0;
  }

  return mv_dst;
}

PredMode de265_image::get_pred_mode_lower_layer(int x, int y) const
{
  assert(ilRefPic != NULL);

  if (equalPictureSizeAndOffsetFlag) {
    // No scaling necessary
    return ilRefPic->get_pred_mode(x, y);
  }

  int PicWidthInSamplesRefLayerY  = ilRefPic->get_width();
  int PicHeightInSamplesRefLayerY = ilRefPic->get_height();

  // 1. The center location ( xPCtr, yPCtr ) of the luma prediction block is derived as follows:
  int xPCtr = x + 8;  // (H 65)
  int yPCtr = y + 8;  // (H 66)

  // 2. The variables xRef and yRef are derived as follows:
  int xRef = (((xPCtr - il_scaling_parameters[0]) * il_scaling_parameters[2] + (1 << 15)) >> 16 ) + il_scaling_parameters[4];  // (H 67)
  int yRef = (((yPCtr - il_scaling_parameters[1]) * il_scaling_parameters[3] + (1 << 15)) >> 16 ) + il_scaling_parameters[5];  // (H 68)

  // 3. The rounded reference layer luma sample location ( xRL, yRL ) is derived as follows:
  int xRL = ((xRef + 4) >> 4) << 4;  // (H 69)
  int yRL = ((yRef + 4) >> 4) << 4;  // (H 70)

  // 4. Upsample the prediction mode. (H 71)
  PredMode rsPredMode;
  if( xRL < 0 || xRL >= PicWidthInSamplesRefLayerY || yRL < 0 || yRL >= PicHeightInSamplesRefLayerY ) {
    rsPredMode = MODE_INTRA;
  }
  else {
    rsPredMode = ilRefPic->get_pred_mode(xRL, yRL);
  }
  return rsPredMode;
}