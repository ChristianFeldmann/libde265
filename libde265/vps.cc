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

#include "vps.h"
#include "util.h"
#include "decctx.h"

#include <assert.h>


de265_error read_vps(decoder_context* ctx, bitreader* reader, video_parameter_set* vps)
{
  int vlc;

  vps->video_parameter_set_id = vlc = get_bits(reader, 4);
  if (vlc >= DE265_MAX_VPS_SETS) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  vps->vps_base_layer_internal_flag = get_bits(reader, 1);
  vps->vps_base_layer_available_flag = get_bits(reader, 1);

  vps->vps_max_layers = vlc = get_bits(reader,6) +1;
  if (vlc > 63) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE; // vps_max_layers_minus1 (range 0...63)

  vps->vps_max_sub_layers = vlc = get_bits(reader,3) +1;
  if (vlc >= MAX_TEMPORAL_SUBLAYERS) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  vps->vps_temporal_id_nesting_flag = get_bits(reader,1);
  
  vlc = get_bits(reader,16);  // vps_reserved_0xffff_16bits
  if (vlc != 0xffff) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  read_profile_tier_level(reader, true, &vps->profile_tier_level,
                          vps->vps_max_sub_layers);

  /*
  read_bit_rate_pic_rate_info(reader, &vps->bit_rate_pic_rate_info,
                              0, vps->vps_max_sub_layers-1);
  */

  vps->vps_sub_layer_ordering_info_present_flag = get_bits(reader,1);
  //assert(vps->vps_max_sub_layers-1 < MAX_TEMPORAL_SUBLAYERS);

  int firstLayerRead = vps->vps_sub_layer_ordering_info_present_flag ? 0 : (vps->vps_max_sub_layers-1);

  for (int i=firstLayerRead;i<vps->vps_max_sub_layers;i++) {
    vps->layer[i].vps_max_dec_pic_buffering = get_uvlc(reader);
    vps->layer[i].vps_max_num_reorder_pics  = get_uvlc(reader);
    vps->layer[i].vps_max_latency_increase  = get_uvlc(reader);

    if (vps->layer[i].vps_max_dec_pic_buffering == UVLC_ERROR ||
        vps->layer[i].vps_max_num_reorder_pics  == UVLC_ERROR ||
        vps->layer[i].vps_max_latency_increase  == UVLC_ERROR) {
      return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
  }

  if (!vps->vps_sub_layer_ordering_info_present_flag) {
    assert(firstLayerRead < MAX_TEMPORAL_SUBLAYERS);

    for (int i=0;i<firstLayerRead;i++) {
      vps->layer[i].vps_max_dec_pic_buffering = vps->layer[firstLayerRead].vps_max_dec_pic_buffering;
      vps->layer[i].vps_max_num_reorder_pics  = vps->layer[firstLayerRead].vps_max_num_reorder_pics;
      vps->layer[i].vps_max_latency_increase  = vps->layer[firstLayerRead].vps_max_latency_increase;
    }
  }


  vps->vps_max_layer_id = get_bits(reader,6);
  vps->vps_num_layer_sets = get_uvlc(reader);
  if (vps->vps_num_layer_sets==UVLC_ERROR) {
    return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
  }
  vps->vps_num_layer_sets++;

  if (vps->vps_num_layer_sets<0 ||
      vps->vps_num_layer_sets>=1024) {
    ctx->add_warning(DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE, false);
    return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
  }

  for (int i=1; i <= vps->vps_num_layer_sets-1; i++)
    for (int j=0; j <= vps->vps_max_layer_id; j++)
      {
        vps->layer_id_included_flag[i][j] = get_bits(reader,1);
      }

  vps->vps_timing_info_present_flag = get_bits(reader,1);

  if (vps->vps_timing_info_present_flag) {
    vps->vps_num_units_in_tick = get_bits(reader,32);
    vps->vps_time_scale        = get_bits(reader,32);
    vps->vps_poc_proportional_to_timing_flag = get_bits(reader,1);

    if (vps->vps_poc_proportional_to_timing_flag) {
      vps->vps_num_ticks_poc_diff_one = get_uvlc(reader);
      if (vps->vps_num_ticks_poc_diff_one == UVLC_ERROR) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
      vps->vps_num_ticks_poc_diff_one++;

      vps->vps_num_hrd_parameters     = get_uvlc(reader);
      if (vps->vps_num_hrd_parameters == UVLC_ERROR) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
      if (vps->vps_num_hrd_parameters >= 1024) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }

      for (int i=0; i<vps->vps_num_hrd_parameters; i++) {
        vps->hrd_layer_set_idx[i] = get_uvlc(reader);
        if (vps->hrd_layer_set_idx[i] == UVLC_ERROR) {
          return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
        }

        if (i > 0) {
          vps->cprms_present_flag[i] = get_bits(reader,1);
        }

        read_hrd_parameters(reader, &vps->hrd_params[i], vps, vps->cprms_present_flag[i], vps->vps_max_sub_layers-1);
      }
    }
  }

  vps->vps_extension_flag = get_bits(reader,1);

  if (vps->vps_extension_flag) {
    // Parser the VPS extension
    read_vps_extension(ctx, reader, vps);

    vps->vps_extension2_flag = get_bits(reader,1);
    if (vps->vps_extension2_flag) {
      // All remaining bits are vps_extension_data_flag
      return DE265_OK;
    }
  }
  else {
    //set_default_vps_extension(vps); ??
  }

  return DE265_OK;
}

de265_error parse_rep_format(decoder_context* ctx, 
                             bitreader* reader, 
                             rep_format *rep)
{
  rep->pic_width_vps_in_luma_samples = get_bits(reader,16);
  rep->pic_height_vps_in_luma_samples = get_bits(reader,16);
  rep->chroma_and_bit_depth_vps_present_flag = get_bits(reader,1);
  if( rep->chroma_and_bit_depth_vps_present_flag ) {
    rep->chroma_format_vps_idc = (de265_chroma)get_bits(reader,2);
    if (rep->chroma_format_vps_idc == 3) {
      rep->separate_colour_plane_vps_flag = get_bits(reader,1);
    }
    rep->bit_depth_vps_luma_minus8 = get_bits(reader,4);
    rep->bit_depth_vps_chroma_minus8 = get_bits(reader,4);
  }
  rep->conformance_window_vps_flag = get_bits(reader,1);
  if( rep->conformance_window_vps_flag ) {
    rep->m_conformanceWindowVps.conf_win_vps_left_offset = get_uvlc(reader);
    rep->m_conformanceWindowVps.conf_win_vps_right_offset = get_uvlc(reader);
    rep->m_conformanceWindowVps.conf_win_vps_top_offset = get_uvlc(reader);
    rep->m_conformanceWindowVps.conf_win_vps_bottom_offset = get_uvlc(reader);
  }

  return DE265_OK;
}

de265_error read_video_signal_info(bitreader* reader,
                            vps_vui_video_signal_info* info)
{
  info->video_vps_format = get_bits(reader,3);
  info->video_full_range_vps_flag = get_bits(reader,1);
  info->colour_primaries_vps = get_bits(reader,8);
  info->transfer_characteristics_vps = get_bits(reader,8);
  info->matrix_coeffs_vps = get_bits(reader,8);

  return DE265_OK;
}

de265_error read_hrd_parameters(bitreader* reader,
                                hrd_parameters *hrd, 
                                video_parameter_set *vps, 
                                bool commonInfPresentFlag, 
                                int maxNumSubLayersMinus1)
{
  if (commonInfPresentFlag) {
    hrd->nal_hrd_parameters_present_flag = get_bits(reader,1);
    hrd->vcl_hrd_parameters_present_flag = get_bits(reader,1);
    if (hrd->nal_hrd_parameters_present_flag || hrd->vcl_hrd_parameters_present_flag) {
      hrd->sub_pic_hrd_params_present_flag = get_bits(reader,1);
      if (hrd->sub_pic_hrd_params_present_flag) {
          hrd->tick_divisor_minus2 = get_bits(reader,8);
          hrd->du_cpb_removal_delay_increment_length_minus1 = get_bits(reader,5);
          hrd->sub_pic_cpb_params_in_pic_timing_sei_flag = get_bits(reader,1);
          hrd->dpb_output_delay_du_length_minus1 = get_bits(reader,5);
      }
      hrd->bit_rate_scale = get_bits(reader,4);
      hrd->cpb_size_scale = get_bits(reader,4);
      if (hrd->sub_pic_hrd_params_present_flag) {
        hrd->cpb_size_du_scale = get_bits(reader,4);
      }
      hrd->initial_cpb_removal_delay_length_minus1 = get_bits(reader,5);
      hrd->au_cpb_removal_delay_length_minus1 = get_bits(reader,5);
      hrd->dpb_output_delay_length_minus1 = get_bits(reader,5);
    }
  }

  for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
    hrd->fixed_pic_rate_general_flag[ i ] = get_bits(reader,1);
    if (!hrd->fixed_pic_rate_general_flag[i]) {
      hrd->fixed_pic_rate_within_cvs_flag[ i ] = get_bits(reader,1);
    }
    if (hrd->fixed_pic_rate_within_cvs_flag[i]) {
      hrd->elemental_duration_in_tc_minus1[i] = get_uvlc(reader);
    }
    else {
      hrd->low_delay_hrd_flag[ i ] = get_bits(reader,1);
    }
    if (!hrd->low_delay_hrd_flag[i]) {
      hrd->cpb_cnt_minus1[i] = get_uvlc(reader);
    }
    if (hrd->nal_hrd_parameters_present_flag) {
      read_sub_layer_hrd_parameters(reader, &hrd->sub_layer_hrd[i], hrd, i);
    }
    if (hrd->vcl_hrd_parameters_present_flag) {
      read_sub_layer_hrd_parameters(reader, &hrd->sub_layer_hrd[i], hrd, i);
    }
  }

  return DE265_OK;
}

de265_error read_sub_layer_hrd_parameters(bitreader* reader,
                                   sub_layer_hrd_parameters *sub_hrd,
                                   hrd_parameters *hrd,
                                   int subLayerId)
{
  int CpbCnt = hrd->cpb_cnt_minus1[subLayerId];
  for( int i = 0; i <= CpbCnt; i++ ) {
    sub_hrd->bit_rate_value_minus1[ i ] = get_uvlc(reader);
    sub_hrd->cpb_size_value_minus1[ i ] = get_uvlc(reader);
    if( hrd->sub_pic_hrd_params_present_flag ) {
      sub_hrd->cpb_size_du_value_minus1[ i ] = get_uvlc(reader);
      sub_hrd->bit_rate_du_value_minus1[ i ] = get_uvlc(reader);
    }
    sub_hrd->cbr_flag[ i ] = get_bits(reader,1);
  }

  return DE265_OK;
}

de265_error read_vps_extension(decoder_context* ctx, bitreader* reader, video_parameter_set *vps)
{
  // Byte alignment (vps_extension_alignment_bit_equal_to_one)
  for (int nrBits = bits_to_byte_boundary(reader); nrBits > 0; nrBits--) {
    if (get_bits(reader, 1) != 1) {
      return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
  }

  video_parameter_set_extension *vps_ext = &vps->vps_extension;

  if (vps->vps_max_layers - 1 > 0 && vps->vps_base_layer_internal_flag) {
    read_profile_tier_level(reader, false, &vps_ext->vps_ext_PTL[0], vps->vps_max_sub_layers);
  }

  int vlc;
  vps_ext->splitting_flag = vlc = get_bits(reader,1);
  int NumScalabilityTypes = 0;
  for (int i = 0; i < 16; i++)
  {
    vps_ext->scalability_mask_flag[i] = vlc = get_bits(reader,1);
    NumScalabilityTypes += vlc;
  }

  for (int j = 0; j < NumScalabilityTypes - vps_ext->splitting_flag; j++)
  {
    vps_ext->dimension_id_len_minus1[j] = get_bits(reader,3);
  }

  if (vps_ext->splitting_flag) {
    int numBits = 0;
    for(int i = 0; i < NumScalabilityTypes - i; i++)
    {
      numBits += vps_ext->dimension_id_len_minus1[i] + 1;
    }
    if (numBits < 6) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    vps_ext->dimension_id_len_minus1[NumScalabilityTypes-1] = 6 - numBits;
    numBits = 6;
  }

  vps_ext->vps_nuh_layer_id_present_flag = get_bits(reader,1);

  int MaxLayersMinus1 = libde265_min(62, vps->vps_max_layers-1);

  vps_ext->layer_id_in_nuh[0] = 0;
  for (int i = 1; i <= MaxLayersMinus1; i++)
  {
    if (vps_ext->vps_nuh_layer_id_present_flag) {
      vps_ext->layer_id_in_nuh[i] = vlc = get_bits(reader,6);
      if (vlc > vps_ext->layer_id_in_nuh[i-1]) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
    else {
      vps_ext->layer_id_in_nuh[i] = i;
    }
    if (!vps_ext->splitting_flag) {
      for (int j = 0; j < NumScalabilityTypes; j++)
      {
        int nrBits = vps_ext->dimension_id_len_minus1[j] + 1;
        vps_ext->dimension_id[i][j] = get_bits(reader,nrBits);
      }
    }
  }

  // Standard:
  // For i from 0 to MaxLayersMinus1, inclusive, the variable LayerIdxInVps[ layer_id_in_nuh[ i ] ] is set equal to i.
  int_1d LayerIdxInVps;
  for (int i = 1; i <= MaxLayersMinus1; i++) {
    LayerIdxInVps[ vps_ext->layer_id_in_nuh[ i ] ] = 1;
  }

  vps_ext->view_id_len = vlc = get_bits(reader,4);
  if (vlc > 0) {
    
    // Standard F.7.4.3.1.1 (F-2) (JCTVC-R1008_v7)
    // The variable ScalabilityId[ i ][ smIdx ] specifying the identifier of the smIdx-th scalability dimension type of the i-th layer, and the variables ViewOrderIdx[ lId ], DependencyId[ lId ], and AuxId[ lId ] specifying the view order index, the spatial/quality scalability identifier, and the auxiliary identifier, respectively, of the layer with nuh_layer_id equal to lId  are derived as follows:
    int NumViews = 1;
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      int lId = vps_ext->layer_id_in_nuh[i];
      int_2d ScalabilityId;
      for (int smIdx = 0, j = 0; smIdx < 16; smIdx++)
      {
        if (vps_ext->scalability_mask_flag[smIdx])
          ScalabilityId[i][smIdx] = vps_ext->dimension_id[i][j++];
        else
          ScalabilityId[i][smIdx] = 0;
      }
      int_1d ViewOrderIdx;
      int_1d DependencyId;
      ViewOrderIdx[ lId ] = ScalabilityId[ i ][ 1 ];
      DependencyId[ lId ] = ScalabilityId[ i ][ 2 ];
      if (i > 0) {
        int newViewFlag  = 1;
        for (int j = 0; j < i; j++) {
          if (ViewOrderIdx[ lId ] == ViewOrderIdx[ vps_ext->layer_id_in_nuh[j] ] )
            newViewFlag  = 0;
        }
        NumViews += newViewFlag ;
      }
    }
    
    for (int i = 0; i < NumViews; i++)
    {
      vps_ext->view_id_val[i] = get_bits(reader,vlc);
    }
  }

  for (int i = 1; i <= MaxLayersMinus1; i++) {
    for (int j = 0; j < i; j++)
    {
      vps_ext->direct_dependency_flag[i][j] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3.1.1 (F-3) (JCTVC-R1008_v7)
  // The variable DependencyFlag[ i ][ j ] is derived as follows:
  bool_2d DependencyFlag;
  for( int i = 0; i <= MaxLayersMinus1; i++ ) {
    for( int j = 0; j <= MaxLayersMinus1; j++ ) {
      DependencyFlag[i][j] = vps_ext->direct_dependency_flag[i][j];
      for (int k = 0; k < i; k++) {
        if (vps_ext->direct_dependency_flag[i][k] && DependencyFlag[k][j])
          DependencyFlag[i][j] = true;
      }
    }
  }

  // Standard F.7.4.3.1.1 (F-4) (JCTVC-R1008_v7)
  // The variables NumDirectRefLayers[ iNuhLId ], IdDirectRefLayer[ iNuhLId ][ d ], NumRefLayers[ iNuhLId ], IdRefLayer[ iNuhLId ][ r ], NumPredictedLayers[ iNuhLId ], and IdPredictedLayer[ iNuhLId ][ p ] are derived as follows:
  int_1d NumDirectRefLayers;
  int_1d NumRefLayers;
  int_1d NumPredictedLayers;
  int_2d IdDirectRefLayer;
  int_2d IdRefLayer;
  int_2d IdPredictedLayer;
  for( int i = 0; i <= MaxLayersMinus1; i++ ) {
    int iNuhLId = vps_ext->layer_id_in_nuh[i];
    int d = 0; int r = 0; int p = 0;
    for( int j = 0; j <= MaxLayersMinus1; j++ ) {
      int jNuhLid = vps_ext->layer_id_in_nuh[ j ];
      if( vps_ext->direct_dependency_flag[i][j] )
        IdDirectRefLayer[ iNuhLId ][ d++ ] = jNuhLid;
      if( DependencyFlag[ i ][ j ] )
        IdRefLayer[ iNuhLId ][ r++ ] = jNuhLid;
      if( DependencyFlag[ j ][ i ] )
        IdPredictedLayer[ iNuhLId ][ p++ ] = jNuhLid;
    }
    NumDirectRefLayers[ iNuhLId ] = d;
    NumRefLayers[ iNuhLId ] = r;
    NumPredictedLayers[ iNuhLId ] = p;
  }

  // Standard F.7.4.3.1.1 (F-5) (JCTVC-R1008_v7)
  // The variables NumIndependentLayers, NumLayersInTreePartition[ i ], and TreePartitionLayerIdList[ i ][ j ] for i in the range of 0 to NumIndependentLayers − 1, inclusive, and j in the range of 0 to NumLayersInTreePartition[ i ] − 1, inclusive, are derived as follows:
  int_2d   TreePartitionLayerIdList;
  int_1d   NumLayersInTreePartition;
  bool     layerIdInListFlag[64];
  for (int i = 0; i <= 63; i++) {
    layerIdInListFlag[i] = false;
  }
  int k = 0;
  for ( int i = 0; i <= MaxLayersMinus1; i++ ) {
    int iNuhLId = vps_ext->layer_id_in_nuh[i];
    if( NumDirectRefLayers[ iNuhLId ] == 0 ) {
      TreePartitionLayerIdList[ k ][ 0 ] = iNuhLId;
      int h = 1;
      for( int j = 0; j < NumPredictedLayers[ iNuhLId ]; j++ ) {
        int predLId = IdPredictedLayer[ iNuhLId ][ j ];
        if( !layerIdInListFlag[ predLId ] ) {
          TreePartitionLayerIdList[ k ][ h++ ] = predLId;
          layerIdInListFlag[ predLId ] = 1;
        }
      }
      NumLayersInTreePartition[ k++ ] = h;
    }
  }
  int NumIndependentLayers = k;

  if (NumIndependentLayers > 1) {
    vps_ext->num_add_layer_sets = get_uvlc(reader);
  }

  for( int i = 0; i < vps_ext->num_add_layer_sets; i++ ) {
    for( int j = 1; j < NumIndependentLayers; j++ ) {
      int nr_bits = ceil_log2(NumLayersInTreePartition[j] + 1);
      vps_ext->highest_layer_idx_plus1[ i ][ j ] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3 - layer_id_included_flag - (7-3) (JCTVC-R1013_v6)
  int_1d NumLayersInIdList;
  int_2d LayerSetLayerIdList;
  // The value of NumLayersInIdList[ 0 ] is set equal to 1 and the value of LayerSetLayerIdList[ 0 ][ 0 ] is set equal to 0
  NumLayersInIdList[0] = 1;
  LayerSetLayerIdList[0][0] = 0;
  //For each value of i in the range of 1 to vps_num_layer_sets_minus1, inclusive, the variable NumLayersInIdList[ i ] and the layer identifier list LayerSetLayerIdList[ i ] are derived as follows:
  for (int i = 1; i <= vps->vps_num_layer_sets-1; i++) {
    int n = 0;
    for (int m = 0; m <= vps->vps_max_layer_id; m++) {
      if( vps->layer_id_included_flag[ i ][ m ] ) {
        LayerSetLayerIdList[ i ][ n++ ] = m;
      }
    }
    NumLayersInIdList[ i ] = n;

    // For each value of i in the range of 1 to vps_num_layer_sets_minus1, inclusive, NumLayersInIdList[ i ] shall be in the range of 1 to vps_max_layers_minus1 + 1, inclusive.
    if (n > vps->vps_max_layers) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
  }

  // Standard F.7.4.3.1.1 (F-8) (JCTVC-R1008_v7)
  // specifies the values of NumLayersInIdList[ vps_num_layer_sets_minus1 + 1 + i ] and LayerSetLayerIdList[ vps_num_layer_sets_minus1 + 1 + i ][ layerNum ] as follows
  for( int i = 0; i < vps_ext->num_add_layer_sets; i++ ) {
    int layerNum = 0;
    int lsIdx = vps->vps_num_layer_sets + i;
    for (int treeIdx = 1; treeIdx < NumIndependentLayers; treeIdx++) {
      for( int layerCnt = 0; layerCnt < vps_ext->highest_layer_idx_plus1[ i ][ treeIdx ]; layerCnt++ ) {
        LayerSetLayerIdList[ lsIdx ][ layerNum++ ] = TreePartitionLayerIdList[ treeIdx ][ layerCnt ];
      }
    }
    NumLayersInIdList[ lsIdx ] = layerNum;
  }

  vps_ext->vps_sub_layers_max_minus1_present_flag = vlc = get_bits(reader,1);
  if (vlc) {
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      vps_ext->sub_layers_vps_max_minus1[i] = get_bits(reader,3);
    }
  }

  vps_ext->max_tid_ref_present_flag = vlc = get_bits(reader,1);
  if (vlc) {
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      for( int j = i + 1; j <= MaxLayersMinus1; j++ ) {
        if (vps_ext->direct_dependency_flag[j][i])
          vps_ext->max_tid_il_ref_pics_plus1[i][j] = get_bits(reader,3);
      }
    }
  }

  vps_ext->default_ref_layers_active_flag = get_bits(reader,1);
  vps_ext->vps_num_profile_tier_level_minus1 = get_uvlc(reader);
  for( int i = vps->vps_base_layer_internal_flag ? 2 : 1;
             i <= vps_ext->vps_num_profile_tier_level_minus1; i++ ) {
    vps_ext->vps_profile_present_flag[i] = vlc = get_bits(reader,1);
    read_profile_tier_level(reader, vlc, &vps_ext->vps_ext_PTL[i], vps->vps_max_sub_layers);
  }

  // Standard F.7.4.3.1.1 (F-6) (JCTVC-R1008_v7)
  // The variable NumLayerSets is derived as follows:
  int NumLayerSets = vps->vps_num_layer_sets + vps_ext->num_add_layer_sets;
  if (NumLayerSets > 1) {
    vps_ext->num_add_olss = get_uvlc(reader);
    vps_ext->default_output_layer_idc = get_bits(reader,2);
  }

  // Local variables
  int_1d OlsIdxToLsIdx;
  bool_2d NecessaryLayerFlag;
  bool_2d OutputLayerFlag;
  int_1d  NumNecessaryLayers;
  int_1d NumOutputLayersInOutputLayerSet;
  int_1d OlsHighestOutputLayerId;

  vps_ext->output_layer_flag[0][0] = 1;
  OutputLayerFlag[0][0] = 1;

  int NumOutputLayerSets = vps_ext->num_add_olss + NumLayerSets;
  for( int i = 1; i < NumOutputLayerSets; i++ ) {
    if (NumLayerSets > 2 && i >= NumLayerSets) {
      int nr_bits = ceil_log2(NumLayerSets - 1);
      vps_ext->layer_set_idx_for_ols_minus1[i] = get_bits(reader,nr_bits);
    }
    int defaultOutputLayerIdc = libde265_min(vps_ext->default_output_layer_idc, 2);

    // Standard F.7.4.3.1.1 (F-10) (JCTVC-R1008_v7)
    // For i in the range of 0 to NumOutputLayerSets − 1, inclusive, the variable OlsIdxToLsIdx[ i ] is derived as specified in the following:
    OlsIdxToLsIdx[ i ] =  ( i < NumLayerSets ) ? i : ( vps_ext->layer_set_idx_for_ols_minus1[ i ] + 1 );

    if (i > vps->vps_num_layer_sets - 1 || defaultOutputLayerIdc == 2) {
      // i in the range of ( defaultOutputLayerIdc  = =  2 ) ? 0 : ( vps_num_layer_sets_minus1 + 1 ) to NumOutputLayerSets − 1, inclusive, 
      for (int j = 0; j < NumLayersInIdList[OlsIdxToLsIdx[i]]; j++) {
        vps_ext->output_layer_flag[ i ][ j ] = get_bits(reader,1);
        OutputLayerFlag[i][j] = vps_ext->output_layer_flag[ i ][ j ];
      }
    }
    else {
      for (int j = 0; j < NumLayersInIdList[OlsIdxToLsIdx[i]]; j++) {
        if (defaultOutputLayerIdc == 0 || defaultOutputLayerIdc == 1) {

          // with nuhLayerIdA being the highest value in LayerSetLayerIdList[ OlsIdxToLsIdx[ i ] ], 
          int nuhLayerIdA = 0;
          for (int l = 0; l < NumLayersInIdList[OlsIdxToLsIdx[i]]; l++) {
            if (LayerSetLayerIdList[OlsIdxToLsIdx[i]][l] > nuhLayerIdA) {
              nuhLayerIdA = LayerSetLayerIdList[OlsIdxToLsIdx[i]][l];
            }
          }
          if (defaultOutputLayerIdc == 0 || LayerSetLayerIdList[OlsIdxToLsIdx[i]][j] == nuhLayerIdA) {
            OutputLayerFlag[i][j] = true;
          }
          else {
            OutputLayerFlag[i][j] = false;
          }
        }
      }
    }

    // Standard F.7.4.3.1.1 (F-11) (JCTVC-R1008_v7)
    // The variable NumOutputLayersInOutputLayerSet[ i ] is derived as follows:
    NumOutputLayersInOutputLayerSet[ i ] = 0;
    for( int j = 0; j < NumLayersInIdList[ OlsIdxToLsIdx[ i ] ]; j++ ) {
      NumOutputLayersInOutputLayerSet[ i ] += OutputLayerFlag[ i ][ j ];
      if( OutputLayerFlag[ i ][ j ] ) {
        OlsHighestOutputLayerId[ i ] = LayerSetLayerIdList[ OlsIdxToLsIdx[ i ] ][ j ];
      }
    }
    
    // Standard F.7.4.3.1.1 (F-12) (JCTVC-R1008_v7)
    // The variables NumNecessaryLayers[ olsIdx ] and NecessaryLayerFlag[ olsIdx ][ lIdx ] are derived as follows:
    for( int olsIdx = 0; olsIdx < NumOutputLayerSets; olsIdx++ ) {
      int lsIdx = OlsIdxToLsIdx[ olsIdx ];
      for( int lsLayerIdx = 0; lsLayerIdx < NumLayersInIdList[ lsIdx ]; lsLayerIdx++ ) {
        NecessaryLayerFlag[ olsIdx ][ lsLayerIdx ] = 0;
      }
      for( int lsLayerIdx = 0; lsLayerIdx < NumLayersInIdList[ lsIdx ]; lsLayerIdx++ ) {
        if( OutputLayerFlag[ olsIdx ][ lsLayerIdx ] ) {
          NecessaryLayerFlag[ olsIdx ][ lsLayerIdx ] = 1;
          int currLayerId = LayerSetLayerIdList[ lsIdx ][ lsLayerIdx ];
          for( int rLsLayerIdx = 0; rLsLayerIdx < lsLayerIdx; rLsLayerIdx++ ) {
            int refLayerId = LayerSetLayerIdList[ lsIdx ][ rLsLayerIdx ];
            if( DependencyFlag[ LayerIdxInVps[ currLayerId ] ][ LayerIdxInVps[ refLayerId ] ] ) {
              NecessaryLayerFlag[ olsIdx ][ rLsLayerIdx ] = 1;
            }
          }
        }
      }
      NumNecessaryLayers[ olsIdx ] = 0;
      for (int lsLayerIdx = 0; lsLayerIdx < NumLayersInIdList[lsIdx]; lsLayerIdx++) {
        NumNecessaryLayers[ olsIdx ] += NecessaryLayerFlag[ olsIdx ][ lsLayerIdx ];
      }
    }

    for (int j = 0; j < NumLayersInIdList[OlsIdxToLsIdx[i]]; j++) {
      if( NecessaryLayerFlag[ i ][ j ] && vps_ext->vps_num_profile_tier_level_minus1 > 0 ) {
        int nr_bits = ceil_log2(vps_ext->vps_num_profile_tier_level_minus1 + 1 );
        vps_ext->profile_tier_level_idx[ i ][ j ] = get_bits(reader,nr_bits);
      }
    }

    if (NumOutputLayersInOutputLayerSet[i] == 1 && NumDirectRefLayers[OlsHighestOutputLayerId[i]] > 0) {
      vps_ext->alt_output_layer_flag[ i ] = get_bits(reader,1);
    }
  }

  vps_ext->vps_num_rep_formats_minus1 = get_uvlc(reader);
  for (int i = 0; i <= vps_ext->vps_num_rep_formats_minus1; i++) {
    de265_error err = parse_rep_format(ctx, reader, &vps_ext->vps_ext_rep_format[i]);
    if ( err != DE265_OK) {
      return err;
    }
  }

  if (vps_ext->vps_num_rep_formats_minus1 > 0) {
    vps_ext->rep_format_idx_present_flag = get_bits(reader,1);
  }
  if (vps_ext->rep_format_idx_present_flag) {
    for( int i = vps->vps_base_layer_internal_flag ? 1 : 0; i <= MaxLayersMinus1; i++ ) {
      int nr_bits = ceil_log2( vps_ext->vps_num_rep_formats_minus1 + 1);
      vps_ext->vps_rep_format_idx[ i ] = get_bits(reader,nr_bits);
    }
  }
  vps_ext->max_one_active_ref_layer_flag = get_bits(reader,1);
  vps_ext->vps_poc_lsb_aligned_flag = get_bits(reader,1);

  for (int i = 1; i <= MaxLayersMinus1; i++) {
    if( NumDirectRefLayers[ vps_ext->layer_id_in_nuh[ i ] ] == 0 ) {
      vps_ext->poc_lsb_not_present_flag[ i ] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3.1.1 (F-9) (JCTVC-R1008_v7)
  // The variable MaxSubLayersInLayerSetMinus1[ i ] is derived as follows:
  int_1d MaxSubLayersInLayerSetMinus1;
  for( int i = 0; i < NumLayerSets; i++ ) {
    int maxSlMinus1 = 0;
    for( int k = 0; k < NumLayersInIdList[ i ]; k++ ) {
      int lId = LayerSetLayerIdList[ i ][ k ];
      maxSlMinus1 = libde265_max( maxSlMinus1, vps_ext->sub_layers_vps_max_minus1[ LayerIdxInVps[ lId ] ] );
    }
    MaxSubLayersInLayerSetMinus1[ i ] = maxSlMinus1;
  }

  // dpb_size() (F.7.3.2.1.3) (JCTVC-R1008_v7)
  dpb_size_table* dpb = &vps_ext->dpb_size_table;
	for( int i = 1; i < NumOutputLayerSets; i++ ) {
    int currLsIdx = OlsIdxToLsIdx[ i ];
    dpb->sub_layer_flag_info_present_flag[ i ] = get_bits(reader,1);
    for( int j = 0; j <= MaxSubLayersInLayerSetMinus1[ currLsIdx ]; j++ ) {
      if (j > 0 && dpb->sub_layer_flag_info_present_flag[i]) {
        dpb->sub_layer_dpb_info_present_flag[ i ][ j ] = get_bits(reader,1);
      }
      dpb->sub_layer_dpb_info_present_flag[i][0] = true;
      if( dpb->sub_layer_dpb_info_present_flag[ i ][ j ] ) {
        for (int k = 0; k < NumLayersInIdList[currLsIdx]; k++) {
          if( NecessaryLayerFlag[ i ][ k ]  &&  ( vps->vps_base_layer_internal_flag || 
            ( LayerSetLayerIdList[ currLsIdx ][ k ] != 0 ) ) ) {
            dpb->max_vps_dec_pic_buffering_minus1[i][k][j] = get_uvlc(reader);
          }
        }
        dpb->max_vps_num_reorder_pics[ i ][ j ] = get_uvlc(reader);
        dpb->max_vps_latency_increase_plus1[ i ][ j ] = get_uvlc(reader);
      }
    }
  }
  // end dpb_size()

  vps_ext->direct_dep_type_len_minus2 = get_uvlc(reader);
  vps_ext->direct_dependency_all_layers_flag = get_bits(reader,1);
  if (vps_ext->direct_dependency_all_layers_flag) {
    int nr_bits = vps_ext->direct_dep_type_len_minus2 + 2;
    vps_ext->direct_dependency_all_layers_type = get_bits(reader,nr_bits);
  }
  else {
    for (int i = vps->vps_base_layer_internal_flag ? 1 : 2; i <= MaxLayersMinus1; i++) {
      for( int j = vps->vps_base_layer_internal_flag ? 0 : 1; j < i; j++ ) {
        if( vps_ext->direct_dependency_flag[ i ][ j ] ) {
          int nr_bits = vps_ext->direct_dep_type_len_minus2 + 2;
          vps_ext->direct_dependency_type[ i ][ j ] = get_bits(reader,nr_bits);
        }
      }
    }
  }

  vps_ext->vps_non_vui_extension_length = get_uvlc(reader);
  for (int i = 1; i <= vps_ext->vps_non_vui_extension_length; i++) {
    int vps_non_vui_extension_data_byte = get_bits(reader, 8);
  }
	vps_ext->vps_vui_present_flag = get_bits(reader,1);

  if (vps_ext->vps_vui_present_flag) {
    // Byte alignment (vps_vui_alignment_bit_equal_to_one)
    for (int nrBits = bits_to_byte_boundary(reader); nrBits > 0; nrBits--) {
      if (get_bits(reader, 1) != 1) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
    }

    // vps_vui()
    vps_vui* vui = &vps_ext->vui;

    vui->cross_layer_pic_type_aligned_flag = get_bits(reader,1);
    vui->cross_layer_irap_aligned_flag = vps_ext->vps_vui_present_flag;
    if (!vui->cross_layer_pic_type_aligned_flag) {
      vui->cross_layer_irap_aligned_flag = get_bits(reader,1);
    }
    if (vui->cross_layer_irap_aligned_flag) {
      vui->all_layers_idr_aligned_flag = get_bits(reader,1);
    }
    vui->bit_rate_present_vps_flag = get_bits(reader,1);
    vui->pic_rate_present_vps_flag = get_bits(reader,1);
    if (vui->bit_rate_present_vps_flag || vui->pic_rate_present_vps_flag) {
      for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i < NumLayerSets; i++) {
        for( int j = 0; j <= MaxSubLayersInLayerSetMinus1[ i ]; j++ ) {
          if (vui->bit_rate_present_vps_flag) {
            vui->bit_rate_present_flag[ i ][ j ] = get_bits(reader,1);
          }
          if (vui->pic_rate_present_vps_flag) {
            vui->pic_rate_present_flag[ i ][ j ] = get_bits(reader,1);
          }
          if( vui->bit_rate_present_flag[ i ][ j ] ) {
            vui->avg_bit_rate[ i ][ j ] = get_bits(reader,16);
            vui->max_bit_rate[ i ][ j ] = get_bits(reader,16);
          }
          if( vui->pic_rate_present_flag[ i ][ j ] ) {
            vui->constant_pic_rate_idc[ i ][ j ] = get_bits(reader,2);
            vui->avg_pic_rate[ i ][ j ] = get_bits(reader,16);
          }
        }
      }
    }

    vui->video_signal_info_idx_present_flag = get_bits(reader,1);
    if (vui->video_signal_info_idx_present_flag) {
      vui->vps_num_video_signal_info_minus1 = get_bits(reader,4);
    }
    for (int i = 0; i <= vui->vps_num_video_signal_info_minus1; i++) {
      read_video_signal_info(reader, &vui->video_signal_info[i] );
    }

    if (vui->video_signal_info_idx_present_flag && vui->vps_num_video_signal_info_minus1 > 0) {
      for( int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= MaxLayersMinus1; i++ ) {
        vui->vps_video_signal_info_idx[ i ] = get_bits(reader,4);
      }
    }

    vui->tiles_not_in_use_flag = get_bits(reader,1);

    if( !vui->tiles_not_in_use_flag ) {
      for( int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= MaxLayersMinus1; i++ ) {
        vui->tiles_in_use_flag[ i ] = get_bits(reader,1);
        if (vui->tiles_in_use_flag[i]) {
          vui->loop_filter_not_across_tiles_flag[ i ] = get_bits(reader,1);
        }
      }
      for (int i = vps->vps_base_layer_internal_flag ? 1 : 2; i <= MaxLayersMinus1; i++) {
        for( int j = 0; j < NumDirectRefLayers[ vps_ext->layer_id_in_nuh[ i ] ]; j++ ) {
          int layerIdx = LayerIdxInVps[ IdDirectRefLayer[ vps_ext->layer_id_in_nuh[ i ] ][ j ] ];
          if( vui->tiles_in_use_flag[ i ] && vui->tiles_in_use_flag[ layerIdx ] ) {
            vui->tile_boundaries_aligned_flag[ i ][ j ] = get_bits(reader,1);
          }
        }
      }
    }

    vui->wpp_not_in_use_flag = get_bits(reader,1);
    if (!vui->wpp_not_in_use_flag) {
      for( int i = vps->vps_base_layer_internal_flag ? 0 : 1; i  <=  MaxLayersMinus1; i++ ) {
        vui->wpp_in_use_flag[ i ] = get_bits(reader,1);
      }
    }
    vui->single_layer_for_non_irap_flag = get_bits(reader,1);
    vui->higher_layer_irap_skip_flag = get_bits(reader,1);
    vui->ilp_restricted_ref_layers_flag = get_bits(reader,1);

    if (vui->ilp_restricted_ref_layers_flag) {
      for( int i = 1; i <= MaxLayersMinus1; i++ ) {
        for( int j = 0; j < NumDirectRefLayers[ vps_ext->layer_id_in_nuh[ i ] ]; j++ ) {
          if( vps->vps_base_layer_internal_flag  ||
              IdDirectRefLayer[ vps_ext->layer_id_in_nuh[ i ] ][ j ] > 0 ) {
            vui->min_spatial_segment_offset_plus1[i][j] = get_uvlc(reader);
            if( vui->min_spatial_segment_offset_plus1[ i ][ j ] > 0 ) {
              vui->ctu_based_offset_enabled_flag[ i ][ j ] = get_bits(reader,1);
              if( vui->ctu_based_offset_enabled_flag[ i ][ j ] ) {
                vui->min_horizontal_ctu_offset_plus1[ i ][ j ] = get_uvlc(reader);
              }
            }
          }
        }
      }
    }

    vui->vps_vui_bsp_hrd_present_flag = get_bits(reader,1);
    if (vui->vps_vui_bsp_hrd_present_flag) {
      
      // vps_vui_bsp_hrd_params( );
      vps_vui_bsp_hrd_params *bsp_hrd = &vui->bsp_hrd_params;

      bsp_hrd->vps_num_add_hrd_params = get_uvlc(reader);
      for( int i = vps->vps_num_hrd_parameters;
           i < vps->vps_num_hrd_parameters + bsp_hrd->vps_num_add_hrd_params; i++ ) {
        if (i > 0) {
          bsp_hrd->cprms_add_present_flag[ i ] = get_bits(reader,1);
        }
        bsp_hrd->num_sub_layer_hrd_minus1[i] = get_uvlc(reader);
        
        read_hrd_parameters(reader, &bsp_hrd->hrd_params[i], vps, bsp_hrd->cprms_add_present_flag[i], bsp_hrd->num_sub_layer_hrd_minus1[i]);
      }

      if( vps->vps_num_hrd_parameters + bsp_hrd->vps_num_add_hrd_params > 0 ) {
        for( int h = 1; h < NumOutputLayerSets; h++ ) {
          bsp_hrd->num_signalled_partitioning_schemes[ h ] = get_uvlc(reader);
          for( int j = 1; j < bsp_hrd->num_signalled_partitioning_schemes[ h ] + 1; j++ ) {
            bsp_hrd->num_partitions_in_scheme_minus1[ h ][ j ] = get_uvlc(reader);
            for( int k = 0; k <= bsp_hrd->num_partitions_in_scheme_minus1[ h ][ j ]; k++ ) {
              for( int r = 0; r < NumLayersInIdList[ OlsIdxToLsIdx[ h ] ]; r++ ) {
                bsp_hrd->layer_included_in_partition_flag[ h ][ j ][ k ][ r ] = get_bits(reader,1);
              }
            }
          }

          for( int i = 0; i < bsp_hrd->num_signalled_partitioning_schemes[ h ] + 1; i++ ) {
            for( int t = 0; t <= MaxSubLayersInLayerSetMinus1[ OlsIdxToLsIdx[ h ] ]; t++ ) {
              bsp_hrd->num_bsp_schedules_minus1[ h ][ i ][ t ] = get_uvlc(reader);
              for (int j = 0; j <= bsp_hrd->num_bsp_schedules_minus1[h][i][t]; j++) {
                for( int k = 0; k <= bsp_hrd->num_partitions_in_scheme_minus1[ h ][ i ]; k++ ) {
                  if( vps->vps_num_hrd_parameters + bsp_hrd->vps_num_add_hrd_params > 1 ) {
                    int nr_bits = ceil_log2( vps->vps_num_hrd_parameters + bsp_hrd->vps_num_add_hrd_params );
                    bsp_hrd->bsp_hrd_idx[ h ][ i ][ t ][ j ][ k ] = get_bits(reader,nr_bits);
                  }
                  bsp_hrd->bsp_sched_idx[ h ][ i ][ t ][ j ][ k ] = get_uvlc(reader);
                }
              }
            }
          }
        }
      }
    } // end vps_vui_bsp_hrd_params( );

    for (int i = 1; i <= MaxLayersMinus1; i++) {
      if( NumDirectRefLayers[ vps_ext->layer_id_in_nuh[ i ] ] == 0 ) {
        vui->base_layer_parameter_set_compatibility_flag[ i ] = get_bits(reader,1);
      }
    }
  } // end vps_vui()

  return DE265_OK;
}

de265_error read_profile_tier_level(bitreader* reader,
                             bool profile_present_flag,
                             struct profile_tier_level* hdr,
                             int max_sub_layers)
{
  if (profile_present_flag) {
    hdr->general_profile_space = get_bits(reader,2);
    hdr->general_tier_flag = get_bits(reader,1);
    hdr->general_profile_idc = get_bits(reader,5);

    for (int i=0; i<32; i++) {
      hdr->general_profile_compatibility_flag[i] = get_bits(reader,1);
    }

    hdr->general_progressive_source_flag = get_bits(reader,1);
    hdr->general_interlaced_source_flag  = get_bits(reader,1);
    hdr->general_non_packed_constraint_flag = get_bits(reader,1);
    hdr->general_frame_only_constraint_flag = get_bits(reader,1);
    skip_bits(reader,44);
  }

  hdr->general_level_idc = get_bits(reader,8);


  for (int i=0; i<max_sub_layers-1; i++)
    {
      hdr->profile[i].sub_layer_profile_present_flag = get_bits(reader,1);
      hdr->profile[i].sub_layer_level_present_flag   = get_bits(reader,1);
    }

  if (max_sub_layers > 1)
    {
      for (int i=max_sub_layers-1; i<8; i++)
        {
          skip_bits(reader,2);
        }
    }

  for (int i=0; i<max_sub_layers-1; i++)
    {
      if (hdr->profile[i].sub_layer_profile_present_flag)
        {
          hdr->profile[i].sub_layer_profile_space = get_bits(reader,2);
          hdr->profile[i].sub_layer_tier_flag = get_bits(reader,1);
          hdr->profile[i].sub_layer_profile_idc = get_bits(reader,5);

          for (int j=0; j<32; j++)
            {
              hdr->profile[i].sub_layer_profile_compatibility_flag[j] = get_bits(reader,1);
            }

          hdr->profile[i].sub_layer_progressive_source_flag = get_bits(reader,1);
          hdr->profile[i].sub_layer_interlaced_source_flag  = get_bits(reader,1);
          hdr->profile[i].sub_layer_non_packed_constraint_flag = get_bits(reader,1);
          hdr->profile[i].sub_layer_frame_only_constraint_flag = get_bits(reader,1);
          skip_bits(reader,44);
        }

      if (hdr->profile[i].sub_layer_level_present_flag)
        {
          hdr->profile[i].sub_layer_level_idc = get_bits(reader,8);
        }
    }

  return DE265_OK;
}


/*
void read_bit_rate_pic_rate_info(bitreader* reader,
                                 struct bit_rate_pic_rate_info* hdr,
                                 int TempLevelLow,
                                 int TempLevelHigh)
{
  for (int i=TempLevelLow; i<=TempLevelHigh; i++) {

    hdr->bit_rate_info_present_flag[i] = get_bits(reader,1);
    hdr->pic_rate_info_present_flag[i] = get_bits(reader,1);

    if (hdr->bit_rate_info_present_flag[i]) {
      hdr->avg_bit_rate[i] = get_bits(reader,16);
      hdr->max_bit_rate[i] = get_bits(reader,16);
    }

    if (hdr->pic_rate_info_present_flag[i]) {
      hdr->constant_pic_rate_idc[i] = get_bits(reader,2);
      hdr->avg_pic_rate[i] = get_bits(reader,16);
    }
  }
}
*/

#define LOG0(t) log2fh(fh, t)
#define LOG1(t,d) log2fh(fh, t,d)
#define LOG2(t,d1,d2) log2fh(fh, t,d1,d2)
#define LOG3(t,d1,d2,d3) log2fh(fh, t,d1,d2,d3)

void dump_vps(video_parameter_set* vps, int fd)
{
  FILE* fh;
  if (fd==1) fh=stdout;
  else if (fd==2) fh=stderr;
  else { return; }

  LOG0("----------------- VPS -----------------\n");
  LOG1("video_parameter_set_id                : %d\n", vps->video_parameter_set_id);
  LOG1("vps_base_layer_internal_flag          : %d\n", vps->vps_base_layer_internal_flag);
  LOG1("vps_base_layer_available_flag         : %d\n", vps->vps_base_layer_available_flag);
  LOG1("vps_max_layers                        : %d\n", vps->vps_max_layers);
  LOG1("vps_max_sub_layers                    : %d\n", vps->vps_max_sub_layers);
  LOG1("vps_temporal_id_nesting_flag          : %d\n", vps->vps_temporal_id_nesting_flag);

  dump_profile_tier_level(&vps->profile_tier_level, vps->vps_max_sub_layers, fh);
  //dump_bit_rate_pic_rate_info(&vps->bit_rate_pic_rate_info, 0, vps->vps_max_sub_layers-1);

  LOG1("vps_sub_layer_ordering_info_present_flag : %d\n",
       vps->vps_sub_layer_ordering_info_present_flag);

  if (vps->vps_sub_layer_ordering_info_present_flag) {
    for (int i=0;i<vps->vps_max_sub_layers;i++) {
      LOG2("layer %d: vps_max_dec_pic_buffering = %d\n",i,vps->layer[i].vps_max_dec_pic_buffering);
      LOG1("         vps_max_num_reorder_pics  = %d\n",vps->layer[i].vps_max_num_reorder_pics);
      LOG1("         vps_max_latency_increase  = %d\n",vps->layer[i].vps_max_latency_increase);
    }
  }
  else {
    LOG1("layer (all): vps_max_dec_pic_buffering = %d\n",vps->layer[0].vps_max_dec_pic_buffering);
    LOG1("             vps_max_num_reorder_pics  = %d\n",vps->layer[0].vps_max_num_reorder_pics);
    LOG1("             vps_max_latency_increase  = %d\n",vps->layer[0].vps_max_latency_increase);
  }


  LOG1("vps_max_layer_id   = %d\n", vps->vps_max_layer_id);
  LOG1("vps_num_layer_sets = %d\n", vps->vps_num_layer_sets);

  for (int i=1; i <= vps->vps_num_layer_sets-1; i++)
    for (int j=0; j <= vps->vps_max_layer_id; j++)
      {
        LOG3("layer_id_included_flag[%d][%d] = %d\n",i,j,
             vps->layer_id_included_flag[i][j]);
      }

  LOG1("vps_timing_info_present_flag = %d\n",
       vps->vps_timing_info_present_flag);

  if (vps->vps_timing_info_present_flag) {
    LOG1("vps_num_units_in_tick = %d\n", vps->vps_num_units_in_tick);
    LOG1("vps_time_scale        = %d\n", vps->vps_time_scale);
    LOG1("vps_poc_proportional_to_timing_flag = %d\n", vps->vps_poc_proportional_to_timing_flag);

    if (vps->vps_poc_proportional_to_timing_flag) {
      LOG1("vps_num_ticks_poc_diff_one = %d\n", vps->vps_num_ticks_poc_diff_one);
      LOG1("vps_num_hrd_parameters     = %d\n", vps->vps_num_hrd_parameters);

      for (int i=0; i<vps->vps_num_hrd_parameters; i++) {
        LOG2("hrd_layer_set_idx[%d] = %d\n", i, vps->hrd_layer_set_idx[i]);

        if (i > 0) {
          LOG2("cprms_present_flag[%d] = %d\n", i, vps->cprms_present_flag[i]);
        }

        //hrd_parameters(cprms_present_flag[i], vps_max_sub_layers_minus1)

        return; // TODO: decode hrd_parameters()
      }
    }
  }

  LOG1("vps_extension_flag = %d\n", vps->vps_extension_flag);
}


void dump_profile_tier_level(const struct profile_tier_level* hdr,
                             int max_sub_layers, FILE* fh)
{
  LOG1("  general_profile_space     : %d\n", hdr->general_profile_space);
  LOG1("  general_tier_flag         : %d\n", hdr->general_tier_flag);
  LOG1("  general_profile_idc       : %d\n", hdr->general_profile_idc);

  LOG0("  general_profile_compatibility_flags: ");
  for (int i=0; i<32; i++) {
    if (i) LOG0("*,");
    LOG1("*%d",hdr->general_profile_compatibility_flag[i]);
  }
  LOG0("*\n");

  LOG1("  general_level_idc         : %d\n", hdr->general_level_idc);

  for (int i=0; i<max_sub_layers-1; i++)
    {
      LOG1("  Profile/Tier/Level [Layer %d]\n",i);

      if (hdr->profile[i].sub_layer_profile_present_flag) {

        LOG1("    sub_layer_profile_space : %d\n",hdr->profile[i].sub_layer_profile_space);
        LOG1("    sub_layer_tier_flag     : %d\n",hdr->profile[i].sub_layer_tier_flag);
        LOG1("    sub_layer_profile_idc   : %d\n",hdr->profile[i].sub_layer_profile_idc);

        LOG0("    sub_layer_profile_compatibility_flags: ");
        for (int j=0; j<32; j++) {
          if (j) LOG0(",");
          LOG1("%d",hdr->profile[i].sub_layer_profile_compatibility_flag[j]);
        }
        LOG0("\n");

        LOG1("    sub_layer_progressive_source_flag : %d\n",hdr->profile[i].sub_layer_progressive_source_flag);
        LOG1("    sub_layer_interlaced_source_flag : %d\n",hdr->profile[i].sub_layer_interlaced_source_flag);
        LOG1("    sub_layer_non_packed_constraint_flag : %d\n",hdr->profile[i].sub_layer_non_packed_constraint_flag);
        LOG1("    sub_layer_frame_only_constraint_flag : %d\n",hdr->profile[i].sub_layer_frame_only_constraint_flag);
      }


      if (hdr->profile[i].sub_layer_level_present_flag) {
        LOG1("    sub_layer_level_idc   : %d\n", hdr->profile[i].sub_layer_level_idc);
      }
    }
}

#undef LOG0
#undef LOG1
#undef LOG2
#undef LOG3


/*
void dump_bit_rate_pic_rate_info(struct bit_rate_pic_rate_info* hdr,
                                 int TempLevelLow,
                                 int TempLevelHigh)
{
  for (int i=TempLevelLow; i<=TempLevelHigh; i++) {

    LOG("  Bitrate [Layer %d]\n", i);

    if (hdr->bit_rate_info_present_flag[i]) {
      LOG("    avg_bit_rate : %d\n", hdr->avg_bit_rate[i]);
      LOG("    max_bit_rate : %d\n", hdr->max_bit_rate[i]);
    }

    if (hdr->pic_rate_info_present_flag[i]) {
      LOG("    constant_pic_rate_idc : %d\n", hdr->constant_pic_rate_idc[i]);
      LOG("    avg_pic_rate[i]       : %d\n", hdr->avg_pic_rate[i]);
    }
  }
}
*/
