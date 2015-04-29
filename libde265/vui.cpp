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

#include "vui.h"
#include "vps.h"

de265_error vps_vui_video_signal_info::read_video_signal_info(bitreader* reader)
{
  video_vps_format = get_bits(reader,3);
  video_full_range_vps_flag = get_bits(reader,1);
  colour_primaries_vps = get_bits(reader,8);
  transfer_characteristics_vps = get_bits(reader,8);
  matrix_coeffs_vps = get_bits(reader,8);

  return DE265_OK;
}

de265_error vps_vui_bsp_hrd_params::read_vps_vui_bsp_hrd_params(bitreader* reader,
                                                                video_parameter_set* vps)
{
  video_parameter_set_extension *vps_ext = &vps->vps_extension;

  vps_num_add_hrd_params = get_uvlc(reader);
  for( int i = vps->vps_num_hrd_parameters;
        i < vps->vps_num_hrd_parameters + vps_num_add_hrd_params; i++ ) {
    if (i > 0) {
      cprms_add_present_flag[ i ] = get_bits(reader,1);
    }
    num_sub_layer_hrd_minus1[i] = get_uvlc(reader);
        
    hrd_params[i].read_hrd_parameters(reader, vps, cprms_add_present_flag[i], num_sub_layer_hrd_minus1[i]);
  }

  if( vps->vps_num_hrd_parameters + vps_num_add_hrd_params > 0 ) {
    for( int h = 1; h < vps_ext->NumOutputLayerSets; h++ ) {
      num_signalled_partitioning_schemes[ h ] = get_uvlc(reader);
      for( int j = 1; j < num_signalled_partitioning_schemes[ h ] + 1; j++ ) {
        num_partitions_in_scheme_minus1[ h ][ j ] = get_uvlc(reader);
        for( int k = 0; k <= num_partitions_in_scheme_minus1[ h ][ j ]; k++ ) {
          for( int r = 0; r < vps_ext->NumLayersInIdList[ vps_ext->OlsIdxToLsIdx[ h ] ]; r++ ) {
            layer_included_in_partition_flag[ h ][ j ][ k ][ r ] = get_bits(reader,1);
          }
        }
      }

      for( int i = 0; i < num_signalled_partitioning_schemes[ h ] + 1; i++ ) {
        for( int t = 0; t <= vps_ext->MaxSubLayersInLayerSetMinus1[ vps_ext->OlsIdxToLsIdx[ h ] ]; t++ ) {
          num_bsp_schedules_minus1[ h ][ i ][ t ] = get_uvlc(reader);
          for (int j = 0; j <= num_bsp_schedules_minus1[h][i][t]; j++) {
            for( int k = 0; k <= num_partitions_in_scheme_minus1[ h ][ i ]; k++ ) {
              if( vps->vps_num_hrd_parameters + vps_num_add_hrd_params > 1 ) {
                int nr_bits = ceil_log2( vps->vps_num_hrd_parameters + vps_num_add_hrd_params );
                bsp_hrd_idx[ h ][ i ][ t ][ j ][ k ] = get_bits(reader,nr_bits);
              }
              bsp_sched_idx[ h ][ i ][ t ][ j ][ k ] = get_uvlc(reader);
            }
          }
        }
      }
    }
  }

  return DE265_OK;
}

de265_error vps_vui::read_vps_vui( bitreader* reader, 
                                   video_parameter_set* vps)
{
  video_parameter_set_extension * vps_ext = &vps->vps_extension;

  cross_layer_pic_type_aligned_flag = get_bits(reader,1);
  cross_layer_irap_aligned_flag = vps_ext->vps_vui_present_flag;
  if (!cross_layer_pic_type_aligned_flag) {
    cross_layer_irap_aligned_flag = get_bits(reader,1);
  }
  if (cross_layer_irap_aligned_flag) {
all_layers_idr_aligned_flag = get_bits(reader, 1);
  }
  bit_rate_present_vps_flag = get_bits(reader, 1);
  pic_rate_present_vps_flag = get_bits(reader, 1);
  if (bit_rate_present_vps_flag || pic_rate_present_vps_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i < vps_ext->NumLayerSets; i++) {
      for (int j = 0; j <= vps_ext->MaxSubLayersInLayerSetMinus1[i]; j++) {
        if (bit_rate_present_vps_flag) {
          bit_rate_present_flag[i][j] = get_bits(reader, 1);
        }
        if (pic_rate_present_vps_flag) {
          pic_rate_present_flag[i][j] = get_bits(reader, 1);
        }
        if (bit_rate_present_flag[i][j]) {
          avg_bit_rate[i][j] = get_bits(reader, 16);
          max_bit_rate[i][j] = get_bits(reader, 16);
        }
        if (pic_rate_present_flag[i][j]) {
          constant_pic_rate_idc[i][j] = get_bits(reader, 2);
          avg_pic_rate[i][j] = get_bits(reader, 16);
        }
      }
    }
  }

  video_signal_info_idx_present_flag = get_bits(reader, 1);
  if (video_signal_info_idx_present_flag) {
    vps_num_video_signal_info_minus1 = get_bits(reader, 4);
  }
  for (int i = 0; i <= vps_num_video_signal_info_minus1; i++) {
    video_signal_info[i].read_video_signal_info(reader);
  }

  if (video_signal_info_idx_present_flag && vps_num_video_signal_info_minus1 > 0) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      vps_video_signal_info_idx[i] = get_bits(reader, 4);
    }
  }

  tiles_not_in_use_flag = get_bits(reader, 1);

  if (!tiles_not_in_use_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      tiles_in_use_flag[i] = get_bits(reader, 1);
      if (tiles_in_use_flag[i]) {
        loop_filter_not_across_tiles_flag[i] = get_bits(reader, 1);
      }
    }
    for (int i = vps->vps_base_layer_internal_flag ? 1 : 2; i <= vps->MaxLayersMinus1; i++) {
      for (int j = 0; j < vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]]; j++) {
        int layerIdx = vps_ext->LayerIdxInVps[vps_ext->IdDirectRefLayer[vps_ext->layer_id_in_nuh[i]][j]];
        if (tiles_in_use_flag[i] && tiles_in_use_flag[layerIdx]) {
          tile_boundaries_aligned_flag[i][j] = get_bits(reader, 1);
        }
      }
    }
  }

  wpp_not_in_use_flag = get_bits(reader, 1);
  if (!wpp_not_in_use_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      wpp_in_use_flag[i] = get_bits(reader, 1);
    }
  }
  single_layer_for_non_irap_flag = get_bits(reader, 1);
  higher_layer_irap_skip_flag = get_bits(reader, 1);
  ilp_restricted_ref_layers_flag = get_bits(reader, 1);

  if (ilp_restricted_ref_layers_flag) {
    for (int i = 1; i <= vps->MaxLayersMinus1; i++) {
      for (int j = 0; j < vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]]; j++) {
        if (vps->vps_base_layer_internal_flag ||
          vps_ext->IdDirectRefLayer[vps_ext->layer_id_in_nuh[i]][j] > 0) {
          min_spatial_segment_offset_plus1[i][j] = get_uvlc(reader);
          if (min_spatial_segment_offset_plus1[i][j] > 0) {
            ctu_based_offset_enabled_flag[i][j] = get_bits(reader, 1);
            if (ctu_based_offset_enabled_flag[i][j]) {
              min_horizontal_ctu_offset_plus1[i][j] = get_uvlc(reader);
            }
          }
        }
      }
    }
  }

  vps_vui_bsp_hrd_present_flag = get_bits(reader, 1);
  if (vps_vui_bsp_hrd_present_flag) {
    bsp_hrd_params.read_vps_vui_bsp_hrd_params(reader, vps);
  }

  for (int i = 1; i <= vps->MaxLayersMinus1; i++) {
    if (vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]] == 0) {
      base_layer_parameter_set_compatibility_flag[i] = get_bits(reader, 1);
    }
  }

  return DE265_OK;
}