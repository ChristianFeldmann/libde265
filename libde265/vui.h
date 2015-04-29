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

#ifndef DE265_VUI_H
#define DE265_VUI_H

#include "libde265/bitstream.h"
#include "libde265/de265.h"
#include "libde265/util.h"

// Structs defined in vps.h
struct video_parameter_set;
struct video_parameter_set_extension;
struct hrd_parameters;

struct vps_vui_bsp_hrd_params {
  de265_error read_vps_vui_bsp_hrd_params(bitreader* reader, video_parameter_set* vps);

  int  vps_num_add_hrd_params;
  
  bool_1d cprms_add_present_flag;
  int_1d  num_sub_layer_hrd_minus1;

  std::map<int, hrd_parameters> hrd_params;

  int_1d  num_signalled_partitioning_schemes;
  int_2d  num_partitions_in_scheme_minus1;
  int_4d  layer_included_in_partition_flag;
  int_3d  num_bsp_schedules_minus1;
  int_5d  bsp_hrd_idx;
  int_5d  bsp_sched_idx;
};

struct vps_vui_video_signal_info {
  de265_error read_video_signal_info(bitreader* reader);

  int  video_vps_format;
  bool video_full_range_vps_flag;
  int  colour_primaries_vps;
  int  transfer_characteristics_vps;
  int  matrix_coeffs_vps;
};

struct vps_vui {
  de265_error read_vps_vui( bitreader* reader, video_parameter_set* vps);

  bool    cross_layer_pic_type_aligned_flag;
  bool    cross_layer_irap_aligned_flag;
  bool    all_layers_idr_aligned_flag;
  bool    bit_rate_present_vps_flag;
  bool    pic_rate_present_vps_flag;
  bool_2d bit_rate_present_flag;
  bool_2d pic_rate_present_flag;
  int_2d  avg_bit_rate;
  int_2d  max_bit_rate;
  int_2d  constant_pic_rate_idc;
  int_2d  avg_pic_rate;
  bool    video_signal_info_idx_present_flag;
  int     vps_num_video_signal_info_minus1;

  std::map<int, vps_vui_video_signal_info> video_signal_info;

  int  vps_video_signal_info_idx[8];
  bool tiles_not_in_use_flag;
  bool tiles_in_use_flag[8];
  bool loop_filter_not_across_tiles_flag[8];
  bool tile_boundaries_aligned_flag[8][8];
  bool wpp_not_in_use_flag;
  bool wpp_in_use_flag[8];
  bool single_layer_for_non_irap_flag;
  bool higher_layer_irap_skip_flag;
  bool ilp_restricted_ref_layers_flag;
  int  min_spatial_segment_offset_plus1[8][8];
  bool ctu_based_offset_enabled_flag   [8][8];
  int  min_horizontal_ctu_offset_plus1 [8][8];
  bool vps_vui_bsp_hrd_present_flag;

  vps_vui_bsp_hrd_params bsp_hrd_params;

  bool base_layer_parameter_set_compatibility_flag[8];
};

#endif