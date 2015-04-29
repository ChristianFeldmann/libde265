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

#ifndef DE265_VPS_H
#define DE265_VPS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_STDBOOL_H
#include <stdbool.h>
#endif

#include "libde265/bitstream.h"
#include "libde265/de265.h"
#include <vector>
#include <map>

#define MAX_TEMPORAL_SUBLAYERS 8
typedef std::map<int, bool> bool_1d;
typedef std::map<int, bool_1d> bool_2d;
typedef std::map<int, int>    int_1d;
typedef std::map<int, int_1d> int_2d;
typedef std::map<int, int_2d> int_3d;
typedef std::map<int, int_3d> int_4d;
typedef std::map<int, int_4d> int_5d;
//typedef std::map<int, std::map<int, std::map<int, std::map<int, int>>>> int_4d;
//typedef std::map<int, std::map<int, std::map<int, std::map<int, std::map<int, int>>>>> int_5d;
typedef std::map<int, std::map<int, char>> char_2d;

/*
struct bit_rate_pic_rate_info {
  char bit_rate_info_present_flag[8];
  char pic_rate_info_present_flag[8];

  int avg_bit_rate[8];
  int max_bit_rate[8];

  char constant_pic_rate_idc[8];
  int  avg_pic_rate[8];

};

void read_bit_rate_pic_rate_info(bitreader* reader,
                                 struct bit_rate_pic_rate_info* hdr,
                                 int TempLevelLow,
                                 int TempLevelHigh);

void dump_bit_rate_pic_rate_info(struct bit_rate_pic_rate_info* hdr,
                                 int TempLevelLow,
                                 int TempLevelHigh);
*/

struct profile_data {
  // --- profile ---

  char sub_layer_profile_present_flag;

  char sub_layer_profile_space;
  char sub_layer_tier_flag;
  char sub_layer_profile_idc;

  char sub_layer_profile_compatibility_flag[32];

  char sub_layer_progressive_source_flag;
  char sub_layer_interlaced_source_flag;
  char sub_layer_non_packed_constraint_flag;
  char sub_layer_frame_only_constraint_flag;


  // --- level ---

  char sub_layer_level_present_flag;
  int  sub_layer_level_idc;
};


struct profile_tier_level {
  int general_profile_space;
  int general_tier_flag;
  int general_profile_idc;

  char general_profile_compatibility_flag[32];

  char general_progressive_source_flag;
  char general_interlaced_source_flag;
  char general_non_packed_constraint_flag;
  char general_frame_only_constraint_flag;

  int general_level_idc;

  struct profile_data profile[MAX_TEMPORAL_SUBLAYERS];
};

typedef struct {
  int_1d  bit_rate_value_minus1;
  int_1d  cpb_size_value_minus1;
  int_1d  cpb_size_du_value_minus1;
  int_1d  bit_rate_du_value_minus1;
  bool_1d cbr_flag;
} sub_layer_hrd_parameters;

typedef struct {
  bool commonInfPresentFlag;

  // Common info
  bool nal_hrd_parameters_present_flag;
  bool vcl_hrd_parameters_present_flag;
  bool sub_pic_hrd_params_present_flag;
  int  tick_divisor_minus2;
  int  du_cpb_removal_delay_increment_length_minus1;
  bool sub_pic_cpb_params_in_pic_timing_sei_flag;
  int  dpb_output_delay_du_length_minus1;
  int  bit_rate_scale;
  int  cpb_size_scale;
  int  cpb_size_du_scale;
  int  initial_cpb_removal_delay_length_minus1;
  int  au_cpb_removal_delay_length_minus1;
  int  dpb_output_delay_length_minus1;
  // end common info

  bool_1d  fixed_pic_rate_general_flag;
  bool_1d  fixed_pic_rate_within_cvs_flag;
  int_1d   elemental_duration_in_tc_minus1;
  bool_1d  low_delay_hrd_flag;
  int_1d   cpb_cnt_minus1;

  std::map<int, sub_layer_hrd_parameters> sub_layer_hrd;

} hrd_parameters;

typedef struct {
  int conf_win_vps_left_offset;
  int conf_win_vps_right_offset;
  int conf_win_vps_top_offset;
  int conf_win_vps_bottom_offset;
} conformance_window;

typedef struct {
  int  pic_width_vps_in_luma_samples;
  int  pic_height_vps_in_luma_samples;
  bool chroma_and_bit_depth_vps_present_flag;

  de265_chroma chroma_format_vps_idc;

  bool separate_colour_plane_vps_flag;
  int  bit_depth_vps_luma_minus8;
  int  bit_depth_vps_chroma_minus8;
  bool conformance_window_vps_flag;

  conformance_window m_conformanceWindowVps;
} rep_format;

typedef struct {
  int  video_vps_format;
  bool video_full_range_vps_flag;
  int  colour_primaries_vps;
  int  transfer_characteristics_vps;
  int  matrix_coeffs_vps;
} vps_vui_video_signal_info;

typedef struct {
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
} vps_vui_bsp_hrd_params;

typedef struct {
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
} vps_vui;

typedef struct {
  bool_1d  sub_layer_flag_info_present_flag;
  bool_2d  sub_layer_dpb_info_present_flag;
  int_3d   max_vps_dec_pic_buffering_minus1;
  int_2d   max_vps_num_reorder_pics;
  int_2d   max_vps_latency_increase_plus1;
} dpb_size_table;

typedef struct {
  std::map<int, profile_tier_level> vps_ext_PTL;

  bool     splitting_flag;
  bool     scalability_mask_flag[16];
  int      dimension_id_len_minus1[16];
  bool     vps_nuh_layer_id_present_flag;
  int      layer_id_in_nuh[8];
  int      dimension_id[8][16];
  int      view_id_len;
  int      view_id_val[8];
  bool     direct_dependency_flag[8][8];
  int      num_add_layer_sets;
  int_2d   highest_layer_idx_plus1;
  bool     vps_sub_layers_max_minus1_present_flag;
  int      sub_layers_vps_max_minus1[8];
  bool     max_tid_ref_present_flag;
  int      max_tid_il_ref_pics_plus1[7][8];
  bool     default_ref_layers_active_flag;
  int      vps_num_profile_tier_level_minus1;
  bool_1d  vps_profile_present_flag;
  
  int      num_add_olss;
  int      default_output_layer_idc;
  int_1d   layer_set_idx_for_ols_minus1;
  bool_2d  output_layer_flag;

  int_2d   profile_tier_level_idx;

  bool_1d  alt_output_layer_flag;
  int      vps_num_rep_formats_minus1;

  std::map<int, rep_format> vps_ext_rep_format;

  bool   rep_format_idx_present_flag;
  int    vps_rep_format_idx[16];
  bool   max_one_active_ref_layer_flag;
  bool   vps_poc_lsb_aligned_flag;
  bool   poc_lsb_not_present_flag[8];

  dpb_size_table dpb_size_table;
  
  int   direct_dep_type_len_minus2;
  bool  direct_dependency_all_layers_flag;
  int   direct_dependency_all_layers_type;
  int   direct_dependency_type[8][8];
  int   vps_non_vui_extension_length;
  bool  vps_vui_present_flag;

  vps_vui vui;
} video_parameter_set_extension;

typedef struct {
  int vps_max_dec_pic_buffering;
  int vps_max_num_reorder_pics;
  int vps_max_latency_increase;
} layer_data;

typedef struct {
  int  video_parameter_set_id;
  bool vps_base_layer_internal_flag;
  bool vps_base_layer_available_flag;
  int  vps_max_layers;
  int  vps_max_sub_layers;
  int  vps_temporal_id_nesting_flag;
  struct profile_tier_level profile_tier_level;
  //struct bit_rate_pic_rate_info bit_rate_pic_rate_info;
  int  vps_sub_layer_ordering_info_present_flag;

  layer_data layer[MAX_TEMPORAL_SUBLAYERS];

  uint8_t vps_max_layer_id;
  int     vps_num_layer_sets;

  char_2d  layer_id_included_flag;

  char     vps_timing_info_present_flag;
  uint32_t vps_num_units_in_tick;
  uint32_t vps_time_scale;
  char     vps_poc_proportional_to_timing_flag;

  int vps_num_ticks_poc_diff_one;
  int vps_num_hrd_parameters;

  uint16_t hrd_layer_set_idx[1024];
  char     cprms_present_flag[1024];

  std::map<int, hrd_parameters> hrd_params;

  bool vps_extension_flag;
  video_parameter_set_extension vps_extension;

  bool vps_extension2_flag;

} video_parameter_set;


de265_error read_profile_tier_level(bitreader* reader,
                             bool profile_present_flag,
                             struct profile_tier_level* hdr,
                             int max_sub_layers);

void dump_profile_tier_level(const struct profile_tier_level* hdr,
                             int max_sub_layers, FILE* fh);

de265_error read_vps(struct decoder_context* ctx, bitreader* reader, video_parameter_set* vps);
de265_error parse_rep_format(decoder_context* ctx, bitreader* reader, rep_format *rep);
de265_error read_video_signal_info(bitreader* reader, vps_vui_video_signal_info* info);

de265_error read_hrd_parameters(bitreader* reader,
                                hrd_parameters *hrd,
                                video_parameter_set *vps, 
                                bool commonInfPresentFlag, 
                                int maxNumSubLayersMinus1);

de265_error read_sub_layer_hrd_parameters(bitreader* reader,
                                   sub_layer_hrd_parameters *sub_hrd,
                                   hrd_parameters *hrd,
                                   int subLayerId);

de265_error read_vps_extension(decoder_context* ctx, bitreader* reader, video_parameter_set *vps);

void dump_vps(video_parameter_set*, int fd);
void dump_profile_tier_level(const struct profile_tier_level* hdr,
                             int max_sub_layers, FILE* fh);

#endif
