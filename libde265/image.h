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

#ifndef DE265_IMAGE_H
#define DE265_IMAGE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <memory>
#ifdef HAVE_STDBOOL_H
#include <stdbool.h>
#endif

#include "libde265/de265.h"
#include "libde265/en265.h"
#include "libde265/sps.h"
#include "libde265/pps.h"
#include "libde265/motion.h"
#include "libde265/threads.h"
#include "libde265/slice.h"
#include "libde265/nal.h"

enum PictureState {
  UnusedForReference,
  UsedForShortTermReference,
  UsedForLongTermReference
};


/* TODO:
   At INTEGRITY_DERIVED_FROM_FAULTY_REFERENCE images, we can check the SEI hash, whether
   the output image is correct despite the faulty reference, and set the state back to correct.
*/
#define INTEGRITY_CORRECT 0
#define INTEGRITY_UNAVAILABLE_REFERENCE 1
#define INTEGRITY_NOT_DECODED 2
#define INTEGRITY_DECODING_ERRORS 3
#define INTEGRITY_DERIVED_FROM_FAULTY_REFERENCE 4

#define SEI_HASH_UNCHECKED 0
#define SEI_HASH_CORRECT   1
#define SEI_HASH_INCORRECT 2

#define TU_FLAG_NONZERO_COEFF  (1<<7)
#define TU_FLAG_SPLIT_TRANSFORM_MASK  0x1F

#define DEBLOCK_FLAG_VERTI (1<<4)
#define DEBLOCK_FLAG_HORIZ (1<<5)
#define DEBLOCK_PB_EDGE_VERTI (1<<6)
#define DEBLOCK_PB_EDGE_HORIZ (1<<7)
#define DEBLOCK_BS_MASK     0x03


#define CTB_PROGRESS_NONE      0
#define CTB_PROGRESS_PREFILTER 1
#define CTB_PROGRESS_DEBLK_V   2
#define CTB_PROGRESS_DEBLK_H   3
#define CTB_PROGRESS_SAO       4

class decoder_context;

template <class DataUnit> class MetaDataArray
{
 public:
  MetaDataArray() { data=NULL; data_size=0; log2unitSize=0; width_in_units=0; height_in_units=0; }
  ~MetaDataArray() { free(data); }

  LIBDE265_CHECK_RESULT bool alloc(int w,int h, int _log2unitSize) {
    int size = w*h;

    if (size != data_size) {
      free(data);
      data = (DataUnit*)malloc(size * sizeof(DataUnit));
      if (data == NULL) {
        data_size = 0;
        return false;
      }
      data_size = size;
    }

    width_in_units = w;
    height_in_units = h;

    log2unitSize = _log2unitSize;

    return data != NULL;
  }

  void clear() {
    if (data) memset(data, 0, sizeof(DataUnit) * data_size);
  }

  void copy(const MetaDataArray *src) {
    if (data && src->data_size == data_size) memcpy(data, src->data, data_size * sizeof(DataUnit));
  }

  const DataUnit& get(int x,int y) const {
    int unitX = x>>log2unitSize;
    int unitY = y>>log2unitSize;

    assert(unitX >= 0 && unitX < width_in_units);
    assert(unitY >= 0 && unitY < height_in_units);

    return data[ unitX + unitY*width_in_units ];
  }

  DataUnit& get(int x,int y) {
    int unitX = x>>log2unitSize;
    int unitY = y>>log2unitSize;

    assert(unitX >= 0 && unitX < width_in_units);
    assert(unitY >= 0 && unitY < height_in_units);

    return data[ unitX + unitY*width_in_units ];
  }

  void set(int x,int y, const DataUnit& d) {
    int unitX = x>>log2unitSize;
    int unitY = y>>log2unitSize;

    assert(unitX >= 0 && unitX < width_in_units);
    assert(unitY >= 0 && unitY < height_in_units);

    data[ unitX + unitY*width_in_units ] = d;
  }

  DataUnit& operator[](int idx) { return data[idx]; }
  const DataUnit& operator[](int idx) const { return data[idx]; }

  int size() const { return data_size; }

  // private:
  DataUnit* data;
  int data_size;
  int log2unitSize;
  int width_in_units;
  int height_in_units;
};

#define SET_CB_BLK(x,y,log2BlkWidth,  Field,value)              \
  int cbX = x >> cb_info.log2unitSize; \
  int cbY = y >> cb_info.log2unitSize; \
  int width = 1 << (log2BlkWidth - cb_info.log2unitSize);           \
  for (int cby=cbY;cby<cbY+width;cby++)                             \
    for (int cbx=cbX;cbx<cbX+width;cbx++)                           \
      {                                                             \
        cb_info[ cbx + cby*cb_info.width_in_units ].Field = value;  \
      }

#define CLEAR_TB_BLK(x,y,log2BlkWidth)              \
  int tuX = x >> tu_info.log2unitSize; \
  int tuY = y >> tu_info.log2unitSize; \
  int width = 1 << (log2BlkWidth - tu_info.log2unitSize);           \
  for (int tuy=tuY;tuy<tuY+width;tuy++)                             \
    for (int tux=tuX;tux<tuX+width;tux++)                           \
      {                                                             \
        tu_info[ tux + tuy*tu_info.width_in_units ] = 0;  \
      }


typedef struct {
  uint16_t SliceAddrRS;
  uint16_t SliceHeaderIndex; // index into array to slice header for this CTB

  sao_info saoInfo;
  bool     deblock;         // this CTB has to be deblocked

  // The following flag helps to quickly check whether we have to
  // check all conditions in the SAO filter or whether we can skip them.
  bool     has_pcm_or_cu_transquant_bypass; // pcm or transquant_bypass is used in this CTB
} CTB_info;


typedef struct {
  uint8_t log2CbSize : 3;   /* [0;6] (1<<log2CbSize) = 64
                               This is only set in the top-left corner of the CB.
                               The other values should be zero.
                               TODO: in the encoder, we have to clear to zero.
                               Used in deblocking and QP-scale decoding */
  uint8_t PartMode : 3;     // (enum PartMode)  [0;7] set only in top-left of CB
                            // Used for spatial merging candidates in current frame
                            // and for deriving interSplitFlag in decoding.

  uint8_t ctDepth : 2;      // [0:3]? (for CTB size 64: 0:64, 1:32, 2:16, 3:8)
                            // Used for decoding/encoding split_cu flag.

  // --- byte boundary ---
  uint8_t PredMode : 2;     // (enum PredMode)  [0;2] must be saved for past images
                            // Used in motion decoding.
  uint8_t pcm_flag : 1;     // Stored for intra-prediction / SAO
  uint8_t cu_transquant_bypass : 1; // Stored for SAO
  // note: 4 bits left

  // --- byte boundary ---
  int8_t  QP_Y;  // Stored for QP prediction

} CB_ref_info;




struct de265_image {
  de265_image();
  ~de265_image();


  de265_error alloc_image(int w,int h, enum de265_chroma c,
                          std::shared_ptr<const seq_parameter_set> sps,
                          bool allocMetadata,
                          decoder_context* dctx,
                          class encoder_context* ectx,
                          de265_PTS pts, void* user_data,
                          bool useCustomAllocFunctions,
                          bool interLayerReferencePicture=false);

  //de265_error alloc_encoder_data(const seq_parameter_set* sps);

  bool is_allocated() const { return pixels[0] != NULL; }
  
  void release();

  void set_headers(std::shared_ptr<video_parameter_set> _vps,
                   std::shared_ptr<seq_parameter_set>   _sps,
                   std::shared_ptr<pic_parameter_set>   _pps) {
    vps = _vps;
    sps = _sps;
    pps = _pps;
  }

  void fill_image(int y,int u,int v);
  de265_error copy_image(const de265_image* src);
  void copy_lines_from(const de265_image* src, int first, int end);
  void exchange_pixel_data_with(de265_image&);

  uint32_t get_ID() const { return ID; }

  /* */ uint8_t* get_image_plane(int cIdx)       { return pixels[cIdx]; }
  const uint8_t* get_image_plane(int cIdx) const { return pixels[cIdx]; }

  // Get the image plane for the prediction and residual signal
  /* */ uint8_t* get_image_plane_prediction(int cIdx)       { return pixels_prediction[cIdx]; }
  const uint8_t* get_image_plane_prediction(int cIdx) const { return pixels_prediction[cIdx]; }
  /* */ uint8_t* get_image_plane_residual(int cIdx)       { return pixels_residual[cIdx]; }
  const uint8_t* get_image_plane_residual(int cIdx) const { return pixels_residual[cIdx]; }
  /* */ uint8_t* get_image_plane_tr_coeff(int cIdx)       { return pixels_tr_coeff[cIdx]; }
  const uint8_t* get_image_plane_tr_coeff(int cIdx) const { return pixels_tr_coeff[cIdx]; }

  void set_image_plane(int cIdx, uint8_t* mem, int stride, void *userdata);

  // Set the image planes for the prediction/residual/trCoeff planes
  void set_image_plane_prediction(int cIdx, uint8_t* mem);
  void set_image_plane_residual(int cIdx, uint8_t* mem);
  void set_image_plane_tr_coeff(int cIdx, uint8_t* mem);

  uint8_t* get_image_plane_at_pos(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return pixels[cIdx] + xpos + ypos*stride;
  }


  /// xpos;ypos in actual plane resolution
  template <class pixel_t>
  pixel_t* get_image_plane_at_pos_NEW(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return (pixel_t*)(pixels[cIdx] + (xpos + ypos*stride)*sizeof(pixel_t));
  }

  template <class pixel_t>
  pixel_t* get_image_plane_prediction_at_pos_NEW(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return (pixel_t*)(pixels_prediction[cIdx] + (xpos + ypos*stride)*sizeof(pixel_t));
  }

  template <class pixel_t>
  pixel_t* get_image_plane_residual_at_pos_NEW(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return (pixel_t*)(pixels_residual[cIdx] + (xpos + ypos*stride)*sizeof(pixel_t));
  }

  template <class pixel_t>
  pixel_t* get_image_plane_tr_coeff_at_pos_NEW(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return (pixel_t*)(pixels_tr_coeff[cIdx] + (xpos + ypos*stride)*sizeof(pixel_t));
  }

  const uint8_t* get_image_plane_at_pos(int cIdx, int xpos,int ypos) const
  {
    int stride = get_image_stride(cIdx);
    return pixels[cIdx] + xpos + ypos*stride;
  }

  void* get_image_plane_at_pos_any_depth(int cIdx, int xpos,int ypos)
  {
    int stride = get_image_stride(cIdx);
    return pixels[cIdx] + ((xpos + ypos*stride) << bpp_shift[cIdx]);
  }

  const void* get_image_plane_at_pos_any_depth(int cIdx, int xpos,int ypos) const
  {
    int stride = get_image_stride(cIdx);
    return pixels[cIdx] + ((xpos + ypos*stride) << bpp_shift[cIdx]);
  }

  /* Number of pixels in one row (not number of bytes).
   */
  int get_image_stride(int cIdx) const
  {
    if (cIdx==0) return stride;
    else         return chroma_stride;
  }

  int get_luma_stride() const { return stride; }
  int get_chroma_stride() const { return chroma_stride; }

  int get_width (int cIdx=0) const { return cIdx==0 ? width  : chroma_width;  }
  int get_height(int cIdx=0) const { return cIdx==0 ? height : chroma_height; }

  enum de265_chroma get_chroma_format() const { return chroma_format; }

  int get_bit_depth(int cIdx) const {
    if (cIdx==0) return sps->BitDepth_Y;
    else         return sps->BitDepth_C;
  }

  int get_bytes_per_pixel(int cIdx) const {
    return (get_bit_depth(cIdx)+7)/8;
  }

  bool high_bit_depth(int cIdx) const {
    return get_bit_depth(cIdx)>8;
  }

  bool can_be_released() const { return PicOutputFlag==false && PicState==UnusedForReference; }


  void add_slice_segment_header(slice_segment_header* shdr) {
    shdr->slice_index = slices.size();
    slices.push_back(shdr);
  }


  bool available_zscan(int xCurr,int yCurr, int xN,int yN) const;

  bool available_pred_blk(int xC,int yC, int nCbS,
                          int xP, int yP, int nPbW, int nPbH, int partIdx,
                          int xN,int yN) const;


  static de265_image_allocation default_image_allocation;

  void printBlk(const char* title, int x0,int y0,int blkSize,int cIdx) const {
    ::printBlk(title, get_image_plane_at_pos(cIdx,x0,y0),
               blkSize, get_image_stride(cIdx));
  }

private:
  uint32_t ID;
  static uint32_t s_next_image_ID;

  uint8_t* pixels[3];
  uint8_t  bpp_shift[3];  // 0 for 8 bit, 1 for 16 bit

  // The pixel data for the prediction and residual components.
  // These are only allocated and filled if explicitly enabled before deocoding is started.
  uint8_t* pixels_prediction[3];
  uint8_t* pixels_residual[3];
  uint8_t* pixels_tr_coeff[3];

  enum de265_chroma chroma_format;

  int width, height;  // size in luma pixels

  int chroma_width, chroma_height;
  int stride, chroma_stride;

  /// Multilayer extension
public:
        bool is_inter_layer_reference_picture()       { return bIlRefPic; }
  const bool is_inter_layer_reference_picture() const { return bIlRefPic; }

  // Set the pointer to the lower layer reference
  void set_lower_layer_picture(const de265_image* src);
  // Get pointers to the reconstruction pixel data from the source
  void get_pixel_pointers_from(de265_image *src);
  
  // Set the parameters needed for metadata upsampling
  void set_inter_layer_metadata_scaling_parameters( int scaling_parameters[10] );
  // Set the given upsampling parameters
  void set_inter_layer_upsampling_parameters(int upsampling_params[2][10]);
  // If set to true, the functions get_pred_mode and get_mv_info will return the upsampled lower layer metadata.
  // Remember to set_inter_layer_metadata_scaling_parameters() and set_lower_layer_picture() first.
  void setInterLayerMotionPredictionEnabled(bool ilPred) { interLayerMotionPredictionEnabled = ilPred; }

  // Set/get if upsampling has to be performed
  void setEqualPictureSizeAndOffsetFlag(bool f) { equalPictureSizeAndOffsetFlag = f; }
  bool getEqualPictureSizeAndOffsetFlag()       { return equalPictureSizeAndOffsetFlag; }
  const bool getEqualPictureSizeAndOffsetFlag() const { return equalPictureSizeAndOffsetFlag; }

  const de265_image* get_il_refPic() const { return ilRefPic; };
  const int *get_upsampling_parameters(bool bChroma) const { return bChroma ? il_upsampling_parameters[1] : il_upsampling_parameters[0]; }
  
  //// Upsample the image from the source using the given upsampling parameters
  //void upsample_image_from(decoder_context* ctx, de265_image* rlPic, int upsampling_params[2][10]);
  //
    
  //// The colour mapping process as specified in clause H.8.1.4.3 is invoked
  //void colour_mapping(decoder_context* ctx, de265_image* rlPic, colour_mapping_table *map, int colourMappingParams[2]);

private:
  bool bIlRefPic;
  bool equalPictureSizeAndOffsetFlag;  // Is upsampling required? (True for SNR scalability)
  bool interLayerMotionPredictionEnabled;
  int  il_scaling_parameters[10];
  int  il_upsampling_parameters[2][10];
  // Pointer to the lower layer reference picture.
  // This is used by get_SliceHeaderIndex to retrive the header of the lower layer reference.
  const de265_image* ilRefPic;

  const PBMotion& get_mv_info_lower_layer(int x, int y) const;  // Get MV from lower layer
  PredMode get_pred_mode_lower_layer(int x, int y) const;       // Get PredMode from lower layer

public:
  uint8_t BitDepth_Y, BitDepth_C;
  uint8_t SubWidthC, SubHeightC;
  std::vector<slice_segment_header*> slices;

public:

  // --- conformance cropping window ---

  uint8_t* pixels_confwin[3];   // pointer to pixels in the conformance window
  uint8_t* pixels_confwin_prediction[3]; 
  uint8_t* pixels_confwin_residual[3];
  uint8_t* pixels_confwin_tr_coeff[3];

  int width_confwin, height_confwin;
  int chroma_width_confwin, chroma_height_confwin;

  // --- decoding info ---

  // If PicOutputFlag==false && PicState==UnusedForReference, image buffer is free.

  int  picture_order_cnt_lsb;
  int  PicOrderCntVal;
  enum PictureState PicState;
  bool PicOutputFlag;

  int32_t removed_at_picture_id;

  const video_parameter_set& get_vps() const { return *vps; }
  const seq_parameter_set& get_sps() const { return *sps; }
  const pic_parameter_set& get_pps() const { return *pps; }

  bool has_vps() const { return (bool)vps; }
  bool has_sps() const { return (bool)sps; }
  bool has_pps() const { return (bool)pps; }

  std::shared_ptr<const seq_parameter_set> get_shared_sps() { return sps; }

  //std::shared_ptr<const seq_parameter_set> get_shared_sps() const { return sps; }
  //std::shared_ptr<const pic_parameter_set> get_shared_pps() const { return pps; }

  decoder_context*    decctx;
  class encoder_context*    encctx;

  int number_of_ctbs() const { return ctb_info.size(); }

private:
  // The image also keeps a reference to VPS/SPS/PPS, because when decoding is delayed,
  // the currently active parameter sets in the decctx might already have been replaced
  // with new parameters.
  std::shared_ptr<const video_parameter_set> vps;
  std::shared_ptr<const seq_parameter_set>   sps;  // the SPS used for decoding this image
  std::shared_ptr<const pic_parameter_set>   pps;  // the PPS used for decoding this image

  MetaDataArray<CTB_info>    ctb_info;
  MetaDataArray<CB_ref_info> cb_info;
  MetaDataArray<PBMotion>    pb_info;
  MetaDataArray<uint8_t>     intraPredMode;
  MetaDataArray<uint8_t>     intraPredModeC;
  MetaDataArray<uint8_t>     tu_info;
  MetaDataArray<uint8_t>     deblk_info;

public:
  // --- meta information ---

  de265_PTS pts;
  void*     user_data;
  void*     plane_user_data[3];  // this is logically attached to the pixel data pointers
  de265_image_allocation image_allocation_functions; // the functions used for memory allocation
  void (*encoder_image_release_func)(en265_encoder_context*,
                                     de265_image*,
                                     void* userdata);

  uint8_t integrity; /* Whether an error occured while the image was decoded.
                        When generated, this is initialized to INTEGRITY_CORRECT,
                        and changed on decoding errors.
                      */
  bool sei_hash_check_result;

  nal_header nal_hdr;

  // --- multi core ---

  de265_progress_lock* ctb_progress; // ctb_info_size

  void mark_all_CTB_progress(int progress) {
    for (int i=0;i<ctb_info.data_size;i++) {
      ctb_progress[i].set_progress(progress);
    }
  }


  void thread_start(int nThreads);
  void thread_run(const thread_task*);
  void thread_blocks();
  void thread_unblocks();
  /* NOTE: you should not access any data in the thread_task after
     calling this, as this function may unlock other threads that
     will push this image to the output queue and free all decoder data. */
  void thread_finishes(const thread_task*);

  void wait_for_progress(thread_task* task, int ctbx,int ctby, int progress);
  void wait_for_progress(thread_task* task, int ctbAddrRS, int progress);

  void wait_for_completion();  // block until image is decoded by background threads
  bool debug_is_completed() const;
  int  num_threads_active() const { return nThreadsRunning + nThreadsBlocked; } // for debug only

  //private:
  int   nThreadsQueued;
  int   nThreadsRunning;
  int   nThreadsBlocked;
  int   nThreadsFinished;
  int   nThreadsTotal;

  // ALIGNED_8(de265_sync_int tasks_pending); // number of tasks pending to complete decoding
  de265_mutex mutex;
  de265_cond  finished_cond;

public:

  /* Clear all CTB/CB/PB decoding data of this image.
     All CTB's processing states are set to 'unprocessed'.
  */
  void clear_metadata();


  // --- CB metadata access ---

  void set_pred_mode(int x,int y, int log2BlkWidth, enum PredMode mode)
  {
    SET_CB_BLK(x,y,log2BlkWidth, PredMode, mode);
  }

  void set_pred_mode(int x,int y, int nPbW,int nPbH, enum PredMode mode);
  
  void fill_pred_mode(enum PredMode mode)
  {
    for (int i=0;i<cb_info.data_size;i++)
      { cb_info[i].PredMode = MODE_INTRA; }
  }

  enum PredMode get_pred_mode(int x,int y) const
  {
    if (bIlRefPic)
      return get_pred_mode_lower_layer(x,y);
    else
      return (enum PredMode)cb_info.get(x,y).PredMode;
  }

  uint8_t get_cu_skip_flag(int x,int y) const
  {
    return get_pred_mode(x,y)==MODE_SKIP;
  }

  void set_pcm_flag(int x,int y, int log2BlkWidth, uint8_t value=1)
  {
    SET_CB_BLK(x,y,log2BlkWidth, pcm_flag, value);

    // TODO: in the encoder, we somewhere have to clear this
    ctb_info.get(x,y).has_pcm_or_cu_transquant_bypass = true;
  }

  int  get_pcm_flag(int x,int y) const
  {
    return cb_info.get(x,y).pcm_flag;
  }

  void set_cu_transquant_bypass(int x,int y, int log2BlkWidth, uint8_t value=1)
  {
    SET_CB_BLK(x,y,log2BlkWidth, cu_transquant_bypass, value);

    // TODO: in the encoder, we somewhere have to clear this
    ctb_info.get(x,y).has_pcm_or_cu_transquant_bypass = true;
  }

  int  get_cu_transquant_bypass(int x,int y) const
  {
    return cb_info.get(x,y).cu_transquant_bypass;
  }

  void set_log2CbSize(int x0, int y0, int log2CbSize, bool fill)
  {
    // In theory, we could assume that remaining cb_info blocks are initialized to zero.
    // But in corrupted streams, slices may overlap and set contradicting log2CbSizes.
    // We also need this for encoding.
    if (fill) {
      SET_CB_BLK(x0,y0,log2CbSize, log2CbSize, 0);
    }

    cb_info.get(x0,y0).log2CbSize = log2CbSize;
  }

  int  get_log2CbSize(int x0, int y0) const
  {
    return (enum PredMode)cb_info.get(x0,y0).log2CbSize;
  }

  // coordinates in CB units
  int  get_log2CbSize_cbUnits(int xCb, int yCb) const
  {
    return (enum PredMode)cb_info[ xCb + yCb*cb_info.width_in_units ].log2CbSize;
  }

  void set_PartMode(int x,int y, enum PartMode mode)
  {
    cb_info.get(x,y).PartMode = mode;
  }

  enum PartMode get_PartMode(int x,int y) const
  {
    return (enum PartMode)cb_info.get(x,y).PartMode;
  }

  void set_ctDepth(int x,int y, int log2BlkWidth, int depth)
  {
    SET_CB_BLK(x,y,log2BlkWidth, ctDepth, depth);
  }

  int get_ctDepth(int x,int y) const
  {
    return cb_info.get(x,y).ctDepth;
  }

  void set_QPY(int x,int y, int log2BlkWidth, int QP_Y)
  {
    SET_CB_BLK (x, y, log2BlkWidth, QP_Y, QP_Y);
  }

  int  get_QPY(int x0,int y0) const
  {
    return cb_info.get(x0,y0).QP_Y;
  }

  // --- TU metadata access ---

  void set_split_transform_flag(int x0,int y0,int trafoDepth)
  {
    tu_info.get(x0,y0) |= (1<<trafoDepth);
  }

  void clear_split_transform_flags(int x0,int y0,int log2CbSize)
  {
    CLEAR_TB_BLK (x0,y0, log2CbSize);
  }

  int  get_split_transform_flag(int x0,int y0,int trafoDepth) const
  {
    return (tu_info.get(x0,y0) & (1<<trafoDepth));
  }

  void set_nonzero_coefficient(int x,int y, int log2TrafoSize)
  {
    const int tuX = x >> tu_info.log2unitSize;
    const int tuY = y >> tu_info.log2unitSize;
    const int width = 1 << (log2TrafoSize - tu_info.log2unitSize);

    for (int tuy=tuY;tuy<tuY+width;tuy++)
      for (int tux=tuX;tux<tuX+width;tux++)
        {
          tu_info[ tux + tuy*tu_info.width_in_units ] |= TU_FLAG_NONZERO_COEFF;
        }
  }

  int  get_nonzero_coefficient(int x,int y) const
  {
    return tu_info.get(x,y) & TU_FLAG_NONZERO_COEFF;
  }


  // --- intraPredMode metadata access ---

  enum IntraPredMode get_IntraPredMode(int x,int y) const
  {
    return (enum IntraPredMode)intraPredMode.get(x,y);
  }

  enum IntraPredMode get_IntraPredMode_atIndex(int idx) const
  {
    return (enum IntraPredMode)intraPredMode[idx];
  }

  void set_IntraPredMode(int PUidx,int log2blkSize, enum IntraPredMode mode)
  {
    int pbSize = 1<<(log2blkSize - intraPredMode.log2unitSize);

    for (int y=0;y<pbSize;y++)
      for (int x=0;x<pbSize;x++)
        intraPredMode[PUidx + x + y*intraPredMode.width_in_units] = mode;
  }

  void set_IntraPredMode(int x0,int y0,int log2blkSize,
                         enum IntraPredMode mode)
  {
    int pbSize = 1<<(log2blkSize - intraPredMode.log2unitSize);
    int PUidx  = (x0>>sps->Log2MinPUSize) + (y0>>sps->Log2MinPUSize)*sps->PicWidthInMinPUs;

    for (int y=0;y<pbSize;y++)
      for (int x=0;x<pbSize;x++) {
        assert(x < sps->PicWidthInMinPUs);
        assert(y < sps->PicHeightInMinPUs);

        int idx = PUidx + x + y*intraPredMode.width_in_units;
        assert(idx<intraPredMode.data_size);
        intraPredMode[idx] = mode;
      }
  }


  enum IntraPredMode get_IntraPredModeC(int x,int y) const
  {
    return (enum IntraPredMode)(intraPredModeC.get(x,y) & 0x3f);
  }

  bool is_IntraPredModeC_Mode4(int x,int y) const
  {
    return intraPredModeC.get(x,y) & 0x80;
  }

  void set_IntraPredModeC(int x0,int y0,int log2blkSize, enum IntraPredMode mode,
                          bool is_mode4)
  {
    uint8_t combinedValue = mode;
    if (is_mode4) combinedValue |= 0x80;

    int pbSize = 1<<(log2blkSize - intraPredMode.log2unitSize);
    int PUidx  = (x0>>sps->Log2MinPUSize) + (y0>>sps->Log2MinPUSize)*sps->PicWidthInMinPUs;

    for (int y=0;y<pbSize;y++)
      for (int x=0;x<pbSize;x++) {
        assert(x<sps->PicWidthInMinPUs);
        assert(y<sps->PicHeightInMinPUs);

        int idx = PUidx + x + y*intraPredModeC.width_in_units;
        assert(idx<intraPredModeC.data_size);
        intraPredModeC[idx] = combinedValue;
      }
  }


  /*
  // NOTE: encoder only
  void set_ChromaIntraPredMode(int x,int y,int log2BlkWidth, enum IntraChromaPredMode mode)
  {
    SET_CB_BLK (x, y, log2BlkWidth, intra_chroma_pred_mode, mode);
  }

  // NOTE: encoder only
  enum IntraChromaPredMode get_ChromaIntraPredMode(int x,int y) const
  {
    return (enum IntraChromaPredMode)(cb_info.get(x,y).intra_chroma_pred_mode);
  }
  */

  // --- CTB metadata access ---

  // address of first CTB in slice
  void set_SliceAddrRS(int ctbX, int ctbY, int SliceAddrRS)
  {
    int idx = ctbX + ctbY*ctb_info.width_in_units;
    ctb_info[idx].SliceAddrRS = SliceAddrRS;
  }

  int  get_SliceAddrRS(int ctbX, int ctbY) const
  {
    return ctb_info[ctbX + ctbY*ctb_info.width_in_units].SliceAddrRS;
  }

  int  get_SliceAddrRS_atCtbRS(int ctbRS) const
  {
    return ctb_info[ctbRS].SliceAddrRS;
  }


  void set_SliceHeaderIndex(int x, int y, int SliceHeaderIndex)
  {
    ctb_info.get(x,y).SliceHeaderIndex = SliceHeaderIndex;
  }

  int  get_SliceHeaderIndex(int x, int y) const
  {
    if (bIlRefPic) {
      // Get the slice header index from the lower layer reference picture

      if (equalPictureSizeAndOffsetFlag) {
        // SNR. No sclaing required.
        return ilRefPic->get_SliceHeaderIndex(x, y);
      }
      else {
        // Get the corresponding lower layer position.
      
        // 1. The center location ( xPCtr, yPCtr ) of the luma prediction block is derived as follows:
        int xPCtr = x + 8;  // (H 65)
        int yPCtr = y + 8;  // (H 66)

        // 2. The variables xRef and yRef are derived as follows:
        int xRef = (((xPCtr - il_scaling_parameters[0]) * il_scaling_parameters[2] + (1 << 15)) >> 16 ) + il_scaling_parameters[4];  // (H 67)
        int yRef = (((yPCtr - il_scaling_parameters[1]) * il_scaling_parameters[3] + (1 << 15)) >> 16 ) + il_scaling_parameters[5];  // (H 68)

        // 3. The rounded reference layer luma sample location ( xRL, yRL ) is derived as follows:
        int xRL = ((xRef + 4) >> 4) << 4;  // (H 69)
        int yRL = ((yRef + 4) >> 4) << 4;  // (H 70)

        return ilRefPic->get_SliceHeaderIndex(xRL, yRL);
      }
    }
    return ctb_info.get(x,y).SliceHeaderIndex;
  }

  int  get_SliceHeaderIndexCtb(int ctbX, int ctbY) const
  {
    return ctb_info[ctbX + ctbY*ctb_info.width_in_units].SliceHeaderIndex;
  }

  int  get_SliceHeaderIndex_atIndex(int ctb) const
  {
    return ctb_info[ctb].SliceHeaderIndex;
  }

  bool is_SliceHeader_available(int x,int y) const
  {
    int idx = ctb_info.get(x,y).SliceHeaderIndex;
    return idx >= 0 && idx < slices.size();
  }

  slice_segment_header* get_SliceHeader(int x, int y)
  {
    int idx = get_SliceHeaderIndex(x,y);
    if (idx >= slices.size()) { return NULL; }
    return slices[idx];
  }

  slice_segment_header* get_SliceHeaderCtb(int ctbX, int ctbY)
  {
    int idx = get_SliceHeaderIndexCtb(ctbX,ctbY);
    if (idx >= slices.size()) { return NULL; }
    return slices[idx];
  }

  const slice_segment_header* get_SliceHeaderCtb(int ctbX, int ctbY) const
  {
    int idx = get_SliceHeaderIndexCtb(ctbX,ctbY);
    if (idx >= slices.size()) { return NULL; }
    return slices[idx];
  }

  void set_sao_info(int ctbX,int ctbY,const sao_info* saoinfo)
  {
    sao_info* sao = &ctb_info[ctbX + ctbY*ctb_info.width_in_units].saoInfo;

    memcpy(sao,
           saoinfo,
           sizeof(sao_info));
  }

  const sao_info* get_sao_info(int ctbX,int ctbY) const
  {
    return &ctb_info[ctbX + ctbY*ctb_info.width_in_units].saoInfo;
  }


  void set_CtbDeblockFlag(int ctbX, int ctbY, bool flag)
  {
    int idx = ctbX + ctbY*ctb_info.width_in_units;
    ctb_info[idx].deblock = flag;
  }

  bool get_CtbDeblockFlag(int ctbX, int ctbY) const
  {
    return ctb_info[ctbX + ctbY*ctb_info.width_in_units].deblock;
  }


  bool get_CTB_has_pcm_or_cu_transquant_bypass(int ctbX,int ctbY) const
  {
    int idx = ctbX + ctbY*ctb_info.width_in_units;
    return ctb_info[idx].has_pcm_or_cu_transquant_bypass;
  }



  // --- DEBLK metadata access ---

  int  get_deblk_width()  const { return deblk_info.width_in_units; }
  int  get_deblk_height() const { return deblk_info.height_in_units; }

  void    set_deblk_flags(int x0,int y0, uint8_t flags)
  {
    const int xd = x0/4;
    const int yd = y0/4;

    if (xd<deblk_info.width_in_units &&
        yd<deblk_info.height_in_units) {
      deblk_info[xd + yd*deblk_info.width_in_units] |= flags;
    }
  }

  uint8_t get_deblk_flags(int x0,int y0) const
  {
    const int xd = x0/4;
    const int yd = y0/4;

    return deblk_info[xd + yd*deblk_info.width_in_units];
  }

  void    set_deblk_bS(int x0,int y0, uint8_t bS)
  {
    uint8_t* data = &deblk_info[x0/4 + y0/4*deblk_info.width_in_units];
    *data &= ~DEBLOCK_BS_MASK;
    *data |= bS;
  }

  uint8_t get_deblk_bS(int x0,int y0) const
  {
    return deblk_info[x0/4 + y0/4*deblk_info.width_in_units] & DEBLOCK_BS_MASK;
  }


  // --- PB metadata access ---

  const PBMotion& get_mv_info(int x,int y) const
  {
    if (bIlRefPic) {
      if (interLayerMotionPredictionEnabled)
        // Get mv info from lower layer image
        return get_mv_info_lower_layer(x,y);
      else {
        // Return invalid vector
        PBMotion mv_dst;
        mv_dst.mv[0].x = 0;
        mv_dst.mv[0].y = 0;
        mv_dst.mv[1].x = 0;
        mv_dst.mv[1].y = 0;
        mv_dst.refIdx[0] = -1;
        mv_dst.refIdx[1] = -1;
        mv_dst.predFlag[0] = 0;
        mv_dst.predFlag[1] = 0;
        return mv_dst;
      }
    }
    else
      return pb_info.get(x,y);
  }

  void set_mv_info(int x,int y, int nPbW,int nPbH, const PBMotion& mv);

  // --- value logging ---

  void printBlk(int x0,int y0, int cIdx, int log2BlkSize);

  // --- functions for retrieving internals ---

  void internals_get_CTB_Info_Layout(int *widthInUnits, int *heightInUnits, int *log2UnitSize) const
  {
	  *widthInUnits = ctb_info.width_in_units;
	  *heightInUnits = ctb_info.height_in_units;
	  *log2UnitSize = ctb_info.log2unitSize;
  }

  void internals_get_sliceIdx(uint16_t *idxArray) const
  {
	  for (int i = 0; i < ctb_info.size(); i++)
	  {
		  idxArray[i] = ctb_info[i].SliceAddrRS;
	  }
  }

  void internals_get_CB_Info_Layout(int *widthInUnits, int *heightInUnits, int *log2UnitSize) const
  {
	*widthInUnits = cb_info.width_in_units;
	*heightInUnits = cb_info.height_in_units;
	*log2UnitSize = cb_info.log2unitSize;
  }

  /* Get the coding block info. The values are stored in the following way:
   * This function returns an array of values, sorted in a grid.
   * You can get the layout of the grid by calling internals_get_CB_Info_Layout.
   * The values are then arranged in a raster scan so the conversion from 
   * a unit position (x,y) to an index is: y*widthInUnits+x

   * Each value in this array contains the following infomation:
   * Bits 0:2 - The cb size in log2 pixels. Every CB can span multiple values in this array.
                Only the top left most unit contains a value. All others are set to 0. (3 Bits)
   * Bits 3:5 - The part size (3 Bits)
   * Bits 6:7 - The prediction mode (0: INTRA, 1: INTER, 2: SKIP) (2 Bits)
   * Bit  8   - PCM flag (1 Bit)
   * Bit  9   - CU Transquant bypass flag (1 Bit)
  */
  void internals_get_CB_info(uint16_t *idxArray) const
  {
	  for (int i = 0; i < cb_info.size(); i++)
	  {
		  uint16_t ret_val = 0;
		  CB_ref_info cb_inf = cb_info[i];

		  int size = cb_inf.log2CbSize;
		  assert( size >= 0 && size < 8 );  // 3 bits
		  ret_val = size;

		  int partMode = cb_inf.PartMode;
		  assert( partMode >= 0 && partMode < 8 );  // 3 bits
		  ret_val |= (partMode << 3);

		  int predMode = cb_inf.PredMode;
		  assert( predMode == 0 || predMode == 1 || predMode == 2 );
		  ret_val |= (predMode << 6);

		  bool pcmFlag = cb_inf.pcm_flag;
		  if (pcmFlag)
			  ret_val |= 256;	// Set bit 8

		  bool bCUTransQuantBypass = cb_inf.cu_transquant_bypass;
		  if (bCUTransQuantBypass )
			  ret_val |= 512;	// Set bit 9

		  idxArray[i] = ret_val;
	  }
  }

  void internals_get_PB_Info_layout(int *widthInUnits, int *heightInUnits, int *log2UnitSize) const
  {
	*widthInUnits = pb_info.width_in_units;
	*heightInUnits = pb_info.height_in_units;
	*log2UnitSize = pb_info.log2unitSize;
  }

  // Get pb info (prediction flags, reference indices and motion vectors)
  // if the refIdx is -1, the reference is unused
  void internals_get_PB_info(int16_t *refPOC0, int16_t *refPOC1, int16_t *x0, int16_t *y0, int16_t *x1, int16_t *y1) const
  {
	  for (int i = 0; i < pb_info.size(); i++)
	  {
		  PBMotion mv = pb_info[i];

		  // Get X/Y
		  int pbSizePix = (1 << pb_info.log2unitSize);
		  int iXPix = (i % pb_info.width_in_units) * pbSizePix;
		  int iYPix = (i / pb_info.width_in_units) * pbSizePix;
		  
		  // Get slice header
		  int shdr_idx = get_SliceHeaderIndex(iXPix, iYPix);
		  int refPOC[2] = {0,0};
		  if (shdr_idx < slices.size()) {
			  // Slice found
			  slice_segment_header *shrd = slices[shdr_idx];
			  
			  // Get reference POCs
			  if (mv.refIdx[0] >= 0 && mv.refIdx[0] < 16)
				refPOC[0] = shrd->RefPicList_POC[0][mv.refIdx[0]];
			  if (mv.refIdx[1] >= 0 && mv.refIdx[1] < 16)
				refPOC[1] = shrd->RefPicList_POC[1][mv.refIdx[1]];
		  }

		  if (mv.predFlag[0] != 1) {
			  // No predictin info for this block
			  refPOC0[i] = -1;
			  x0[i] = 0;
			  y0[i] = 0;
		  }
		  else {
			  refPOC0[i] = refPOC[0];
			  x0[i] = mv.mv[0].x;
			  y0[i] = mv.mv[0].y;
		  }

		  if (mv.predFlag[1] != 1) {
			  // No predictin info for this block
			  refPOC1[i] = -1;
			  x1[i] = 0;
			  y1[i] = 0;
		  }
		  else {
			  refPOC1[i] = refPOC[1];
			  x1[i] = mv.mv[1].x;
			  y1[i] = mv.mv[1].y;
		  }
	  }
  }

  // Get intra direction metadata array layout
  void internals_get_IntraDir_Info_Layer(int *widthInUnits, int *heightInUnits, int *log2UnitSize) const
  {
	*widthInUnits = intraPredMode.width_in_units;
	*heightInUnits = intraPredMode.height_in_units;
	*log2UnitSize = intraPredMode.log2unitSize;
  }

  // Get the luma and chroma intra dir array
  void internals_get_intraDir_info(uint8_t *intraDir, uint8_t *intraDirChroma) const
  {
	  for (int i = 0; i < intraPredMode.size(); i++)
	  {
		  intraDir[i] = intraPredMode[i];
		  intraDirChroma[i] = intraPredModeC[i];
	  }
  }

  // Get intra direction metadata array layout
  void internals_get_TUInfo_Info_layout(int *widthInUnits, int *heightInUnits, int *log2UnitSize) const
  {
	*widthInUnits = tu_info.width_in_units;
	*heightInUnits = tu_info.height_in_units;
	*log2UnitSize = tu_info.log2unitSize;
  }

  // Get the luma and chroma intra dir array
  void internals_get_TUInfo_info(uint8_t *tuInfo) const
  {
	  for (int i = 0; i < tu_info.size(); i++)
	  {
		  tuInfo[i] = tu_info[i];
	  }
  }

  void internals_set_residual_CU_signal_to_zero(int x, int y, int w)
  {
    internals_set_CU_signal_to_zero(x, y, w, 1);
  }

  void internals_set_tr_coeff_CU_signal_to_zero(int x, int y, int w)
  {
    internals_set_CU_signal_to_zero(x, y, w, 2);
  }

private:
  // Set the residual data for the CU (Luma and Chroma) to the 'middle' value depending on the bit depth (128 for 8 bit)
  // The parameters are: The x and y pixel position of the CU in the frame. w: the width/height of the CU.
  void internals_set_CU_signal_to_zero(int x, int y, int w, int plane)
  {
    for (int cIdx = 0; cIdx < 3; cIdx++)
    {
      if ((plane == 0 && get_image_plane_prediction(cIdx) == NULL) ||
          (plane == 1 && get_image_plane_residual(cIdx) == NULL) ||
          (plane == 2 && get_image_plane_tr_coeff(cIdx) == NULL))
        continue;

      int bd = get_bit_depth(cIdx);
      int stride = get_image_stride(cIdx);
      int xT = (cIdx == 0) ? x : x / SubWidthC;
      int yT = (cIdx == 0) ? y : y / SubHeightC;
      int width = (cIdx == 0) ? w : w / SubHeightC;
      int height = (cIdx == 0) ? w : w / SubHeightC;
      int imageHeight = get_height(cIdx);
      if (imageHeight < yT + height)
        height = imageHeight - yT;
      if (bd <= 8)
      {
        // One byte per value
        uint8_t *pix;
        if (plane == 0)
          pix = get_image_plane_prediction_at_pos_NEW<uint8_t>(cIdx, xT,yT);
        else if (plane == 1)
          pix = get_image_plane_residual_at_pos_NEW<uint8_t>(cIdx, xT,yT);
        else if (plane == 2)
          pix = get_image_plane_tr_coeff_at_pos_NEW<uint8_t>(cIdx, xT,yT);
        else
          return;

        uint8_t middleValue = (1 << (bd-1));
        for (int y = 0; y < height; y++)
          memset(pix + y*stride, middleValue, width);
      }
      else if (bd <= 16)
      {
        // Two bytes per value
        uint16_t *pix;
        if (plane == 0)
          pix = get_image_plane_prediction_at_pos_NEW<uint16_t>(cIdx, xT,yT);
        else if (plane == 1)
          pix = get_image_plane_residual_at_pos_NEW<uint16_t>(cIdx, xT,yT);
        else if (plane == 2)
          pix = get_image_plane_tr_coeff_at_pos_NEW<uint16_t>(cIdx, xT,yT);
        else
          return;

        // There is no memset 16 function so set the middle value manually
        uint16_t middleValue = (1 << (bd-1));
        for (int y = 0; y < height; y++)
          for (int x = 0; x < width; x++)
            pix[x + y * stride] = middleValue;
      }
      else
        // > 16 bit not supported
        assert(false);
    }
  }
};


#endif
