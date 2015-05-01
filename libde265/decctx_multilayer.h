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

#ifndef DE265_DECCTX_MULTILAYER_H
#define DE265_DECCTX_MULTILAYER_H

#include "libde265/vps.h"
#include "libde265/sps.h"
#include "libde265/pps.h"
#include "libde265/nal.h"
#include "libde265/slice.h"
#include "libde265/image.h"
#include "libde265/motion.h"
#include "libde265/de265.h"
#include "libde265/dpb.h"
#include "libde265/sei.h"
#include "libde265/threads.h"
#include "libde265/acceleration.h"
#include "libde265/nal-parser.h"
#include "libde265/decctx.h"

class decoder_context_multilayer {
 public:
  decoder_context_multilayer();
  ~decoder_context_multilayer();

  de265_error get_warning();
  void reset();
  de265_error decode(int* more);

  // Get the layer decoder context with the given layer_id.
  // Create it if does not exist yet.
  decoder_context* get_layer_dec(int layer_id);

  // Get the number of picture in the output queue (sum over all layers)
  int num_pictures_in_output_queue();
  // Get next output picture. Lowest layers will be returned first. Return the layer that the image is from.
  de265_image* get_next_picture_in_output_queue(int* layerID);
  // Pop next output picture. Lowest layers will be poped first.
  void pop_next_picture_in_output_queue();

  // Flush data
  void flush_data();

  // --- parameters ---

  bool param_sei_check_hash;
  bool param_conceal_stream_errors;
  bool param_suppress_faulty_pictures;

  int  param_sps_headers_fd;
  int  param_vps_headers_fd;
  int  param_pps_headers_fd;
  int  param_slice_headers_fd;

  bool param_disable_deblocking;
  bool param_disable_sao;

  de265_acceleration accelerationFunction;

  NAL_Parser nal_parser;

  int limit_HighestTid;    // never switch to a layer above this one

  multilayer_decoder_parameters ml_dec_params;

protected:
  decoder_context* layer_decoders[MAX_LAYER_ID];
  int num_layer_decoders;

  std::vector<image_unit*> image_units[MAX_LAYER_ID]; // One image list per layer for the decoded images
};


#endif
