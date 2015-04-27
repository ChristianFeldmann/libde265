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

#include "decctx_multilayer.h"
#include "util.h"
#include "sao.h"
#include "sei.h"
#include "deblock.h"

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "fallback.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_SSE4_1
#include "x86/sse.h"
#endif

#ifdef HAVE_ARM
#include "arm/arm.h"
#endif

decoder_context_multilayer::decoder_context_multilayer()
{
  //memset(ctx, 0, sizeof(decoder_context_multilayer));

  // --- parameters ---

  param_sei_check_hash = false;
  param_conceal_stream_errors = true;
  param_suppress_faulty_pictures = false;

  param_disable_deblocking = false;
  param_disable_sao = false;
  //param_disable_mc_residual_idct = false;
  //param_disable_intra_residual_idct = false;

  // --- processing ---

  param_sps_headers_fd = -1;
  param_vps_headers_fd = -1;
  param_pps_headers_fd = -1;
  param_slice_headers_fd = -1;

  nrLayersToDecode = 1;

  // Init by creating one decoder (there has to be at least one layer in the bitstream)
  num_layer_decoders = 1;
  layer_decoders[0] = new decoder_context;
  layer_decoders[0]->setLayerID(0);
  for (int i=1; i<MAX_NR_LAYERS; i++)
    layer_decoders[i] = NULL;
}

decoder_context_multilayer::~decoder_context_multilayer()
{
}

de265_error decoder_context_multilayer::get_warning()
{
  // Loop through all layer decoders and return the first error that is not DE265_OK
  for (int i = 0; i < num_layer_decoders; i++)
  {
    de265_error err = layer_decoders[i]->get_warning();
    if (err != DE265_OK)
      return err;
  }
  return DE265_OK;
}

void decoder_context_multilayer::reset()
{
  // Reset all layer decoders
  // Todo: Not shure if this is the correct way to handle this
  for (int i = 0; i < num_layer_decoders; i++)
  {
    layer_decoders[i]->reset();
  }
}

int decoder_context_multilayer::num_pictures_in_output_queue()
{
  int nrPics = 0;
  for (int i = 0; i < num_layer_decoders; i++)
  {
    nrPics += layer_decoders[i]->num_pictures_in_output_queue();
  }
  return nrPics;
}

de265_image* decoder_context_multilayer::get_next_picture_in_output_queue(int* layerID)
{
  de265_image* img = NULL;
  for (int i = 0; i < num_layer_decoders; i++)
  {
    img = layer_decoders[i]->get_next_picture_in_output_queue();
    if (img != NULL) {
      *layerID = i;
      return img;
    }
  }
  *layerID = -1;
  return NULL;
}

void decoder_context_multilayer::pop_next_picture_in_output_queue()
{
  de265_image* img = NULL;
  for (int i = 0; i < num_layer_decoders; i++) {
    img = layer_decoders[i]->get_next_picture_in_output_queue();
    if (img != NULL) {
      layer_decoders[i]->pop_next_picture_in_output_queue();
      return;
    }
  }
}

de265_error decoder_context_multilayer::decode(int* more)
{
  // First push all NAL units that we have in the buffer to the different decoders
  if (nal_parser.get_NAL_queue_length()) {
    NAL_unit* nal = nal_parser.pop_from_NAL_queue();
    while (nal) {
      // Parse the header
      bitreader reader;
      bitreader_init(&reader, nal->data(), nal->size());
      nal_header nal_hdr;
      nal_read_header(&reader, &nal_hdr);

      if (nal_hdr.nuh_layer_id >= nrLayersToDecode) {
        // Discard all NAL units with nuh_layer_id > (nrLayersToDecode-1)
        nal_parser.free_NAL_unit(nal);
      }
      else {
        decoder_context* layerCtx = layer_decoders[nal_hdr.nuh_layer_id];
        if (layerCtx == NULL) {
          // The decoder for the layer nuh_layer_id does not exits yet.
          // Create it
          if (nal_hdr.nuh_layer_id != num_layer_decoders) {
            // NAL unit layer ids should be continuous in the bitstream.
            return DE265_ERROR_INVALID_LAYER_ID;
          }
          layer_decoders[nal_hdr.nuh_layer_id] = new decoder_context;
          layer_decoders[nal_hdr.nuh_layer_id]->setLayerID(nal_hdr.nuh_layer_id);
          num_layer_decoders++;

          // Get the pointer to the new decoder
          layerCtx = layer_decoders[nal_hdr.nuh_layer_id];
        }

        // Push the NAL unit to the correct layer decoder
        layerCtx->nal_parser.push_to_NAL_queue(nal); // The layer Ctx now owns this NAL unit and will take care of deleting it
      }

      // Process the next NAL unit
      nal = nal_parser.pop_from_NAL_queue();
    }
  }

  // Call the decode function of all layers
  *more = 0;
  int         layer_more = 1;
  de265_error layer_error;
  de265_error ret_error = DE265_OK;
  for (int i = 0; i < num_layer_decoders; i++) {
    layer_error = layer_decoders[i]->decode(&layer_more);

    if (layer_more != 0)
      *more = 1;

    if (layer_error != DE265_OK)
      ret_error = layer_error;
  }

  return ret_error;
}

void decoder_context_multilayer::flush_data()
{
  // Flush data and mark as end of stream
  nal_parser.flush_data();
  nal_parser.mark_end_of_stream();
  
  // Also mark end of stream for all decoders
  for (int i = 0; i < num_layer_decoders; i++) {
    layer_decoders[i]->nal_parser.mark_end_of_frame();
  }
}