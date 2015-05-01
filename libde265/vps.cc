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


de265_error video_parameter_set::read_vps(decoder_context* ctx, bitreader* reader)
{
  int vlc;

  video_parameter_set_id = vlc = get_bits(reader, 4);
  if (vlc >= DE265_MAX_VPS_SETS) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  vps_base_layer_internal_flag = get_bits(reader, 1);
  vps_base_layer_available_flag = get_bits(reader, 1);

  vps_max_layers = vlc = get_bits(reader,6) +1;
  if (vlc > 63) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE; // vps_max_layers_minus1 (range 0...63)
  MaxLayersMinus1 = libde265_min(62, vps_max_layers-1);

  vps_max_sub_layers = vlc = get_bits(reader,3) +1;
  if (vlc >= MAX_TEMPORAL_SUBLAYERS) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  vps_temporal_id_nesting_flag = get_bits(reader,1);
  
  vlc = get_bits(reader,16);  // vps_reserved_0xffff_16bits
  if (vlc != 0xffff) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;

  profile_tier_level.read_profile_tier_level(reader, true, vps_max_sub_layers);

  /*
  read_bit_rate_pic_rate_info(reader, &bit_rate_pic_rate_info,
                              0, vps_max_sub_layers-1);
  */

  vps_sub_layer_ordering_info_present_flag = get_bits(reader,1);
  //assert(vps_max_sub_layers-1 < MAX_TEMPORAL_SUBLAYERS);

  int firstLayerRead = vps_sub_layer_ordering_info_present_flag ? 0 : (vps_max_sub_layers-1);

  for (int i=firstLayerRead;i<vps_max_sub_layers;i++) {
    layer[i].vps_max_dec_pic_buffering = get_uvlc(reader);
    layer[i].vps_max_num_reorder_pics  = get_uvlc(reader);
    layer[i].vps_max_latency_increase  = get_uvlc(reader);

    if (layer[i].vps_max_dec_pic_buffering == UVLC_ERROR ||
        layer[i].vps_max_num_reorder_pics  == UVLC_ERROR ||
        layer[i].vps_max_latency_increase  == UVLC_ERROR) {
      return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
  }

  if (!vps_sub_layer_ordering_info_present_flag) {
    assert(firstLayerRead < MAX_TEMPORAL_SUBLAYERS);

    for (int i=0;i<firstLayerRead;i++) {
      layer[i].vps_max_dec_pic_buffering = layer[firstLayerRead].vps_max_dec_pic_buffering;
      layer[i].vps_max_num_reorder_pics  = layer[firstLayerRead].vps_max_num_reorder_pics;
      layer[i].vps_max_latency_increase  = layer[firstLayerRead].vps_max_latency_increase;
    }
  }


  vps_max_layer_id = get_bits(reader,6);
  vps_num_layer_sets = get_uvlc(reader);
  if (vps_num_layer_sets==UVLC_ERROR) {
    return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
  }
  vps_num_layer_sets++;

  if (vps_num_layer_sets<0 ||
      vps_num_layer_sets>=1024) {
    ctx->add_warning(DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE, false);
    return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
  }

  for (int i=1; i <= vps_num_layer_sets-1; i++)
    for (int j=0; j <= vps_max_layer_id; j++)
      {
        layer_id_included_flag[i][j] = get_bits(reader,1);
      }

  vps_timing_info_present_flag = get_bits(reader,1);

  if (vps_timing_info_present_flag) {
    vps_num_units_in_tick = get_bits(reader,32);
    vps_time_scale        = get_bits(reader,32);
    vps_poc_proportional_to_timing_flag = get_bits(reader,1);

    if (vps_poc_proportional_to_timing_flag) {
      vps_num_ticks_poc_diff_one = get_uvlc(reader);
      if (vps_num_ticks_poc_diff_one == UVLC_ERROR) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
      vps_num_ticks_poc_diff_one++;

      vps_num_hrd_parameters     = get_uvlc(reader);
      if (vps_num_hrd_parameters == UVLC_ERROR) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
      if (vps_num_hrd_parameters >= 1024) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }

      for (int i=0; i<vps_num_hrd_parameters; i++) {
        hrd_layer_set_idx[i] = get_uvlc(reader);
        if (hrd_layer_set_idx[i] == UVLC_ERROR) {
          return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
        }

        if (i > 0) {
          cprms_present_flag[i] = get_bits(reader,1);
        }

        hrd_params[i].read_hrd_parameters(reader, cprms_present_flag[i], vps_max_sub_layers-1);
      }
    }
  }

  vps_extension_flag = get_bits(reader,1);

  if (vps_extension_flag) {
    // Parser the VPS extension
    vps_extension.read_vps_extension(ctx, reader, this);

    vps_extension2_flag = get_bits(reader,1);
    if (vps_extension2_flag) {
      // All remaining bits are vps_extension_data_flag
      return DE265_OK;
    }
  }
  else {
    //set_default_vps_extension(vps); ??
  }

  return DE265_OK;
}

de265_error rep_format::parse_rep_format( bitreader* reader)
{
  pic_width_vps_in_luma_samples = get_bits(reader,16);
  pic_height_vps_in_luma_samples = get_bits(reader,16);
  chroma_and_bit_depth_vps_present_flag = get_bits(reader,1);
  if( chroma_and_bit_depth_vps_present_flag ) {
    chroma_format_vps_idc = (de265_chroma)get_bits(reader,2);
    if (chroma_format_vps_idc == 3) {
      separate_colour_plane_vps_flag = get_bits(reader,1);
    }
    bit_depth_vps_luma_minus8 = get_bits(reader,4);
    bit_depth_vps_chroma_minus8 = get_bits(reader,4);
  }
  conformance_window_vps_flag = get_bits(reader,1);
  if( conformance_window_vps_flag ) {
    m_conformanceWindowVps.conf_win_vps_left_offset = get_uvlc(reader);
    m_conformanceWindowVps.conf_win_vps_right_offset = get_uvlc(reader);
    m_conformanceWindowVps.conf_win_vps_top_offset = get_uvlc(reader);
    m_conformanceWindowVps.conf_win_vps_bottom_offset = get_uvlc(reader);
  }

  return DE265_OK;
}

de265_error video_parameter_set_extension::read_vps_extension(decoder_context* ctx, bitreader* reader, video_parameter_set *vps)
{
  // Byte alignment (vps_extension_alignment_bit_equal_to_one)
  for (int nrBits = bits_to_byte_boundary(reader); nrBits > 0; nrBits--) {
    if (get_bits(reader, 1) != 1) {
      return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
  }

  if (vps->vps_max_layers - 1 > 0 && vps->vps_base_layer_internal_flag) {
    vps_ext_PTL[0].read_profile_tier_level(reader, false, vps->vps_max_sub_layers);
  }

  int vlc;
  splitting_flag = vlc = get_bits(reader,1);
  int NumScalabilityTypes = 0;
  for (int i = 0; i < 16; i++)
  {
    scalability_mask_flag[i] = vlc = get_bits(reader,1);
    NumScalabilityTypes += vlc;
  }

  for (int j = 0; j < NumScalabilityTypes - splitting_flag; j++)
  {
    dimension_id_len_minus1[j] = get_bits(reader,3);
  }

  if (splitting_flag) {
    int numBits = 0;
    for(int i = 0; i < NumScalabilityTypes - i; i++)
    {
      numBits += dimension_id_len_minus1[i] + 1;
    }
    if (numBits < 6) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    dimension_id_len_minus1[NumScalabilityTypes-1] = 6 - numBits;
    numBits = 6;
  }

  vps_nuh_layer_id_present_flag = get_bits(reader,1);

  int MaxLayersMinus1 = libde265_min(62, vps->vps_max_layers-1);

  layer_id_in_nuh[0] = 0;
  for (int i = 1; i <= MaxLayersMinus1; i++)
  {
    if (vps_nuh_layer_id_present_flag) {
      layer_id_in_nuh[i] = vlc = get_bits(reader,6);
      if (vlc > layer_id_in_nuh[i-1]) return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
    }
    else {
      layer_id_in_nuh[i] = i;
    }
    if (!splitting_flag) {
      for (int j = 0; j < NumScalabilityTypes; j++)
      {
        int nrBits = dimension_id_len_minus1[j] + 1;
        dimension_id[i][j] = get_bits(reader,nrBits);
      }
    }
  }

  // Standard:
  // For i from 0 to MaxLayersMinus1, inclusive, the variable LayerIdxInVps[ layer_id_in_nuh[ i ] ] is set equal to i.
  for (int i = 1; i <= MaxLayersMinus1; i++) {
    LayerIdxInVps[ layer_id_in_nuh[ i ] ] = 1;
  }

  view_id_len = vlc = get_bits(reader,4);
  if (vlc > 0) {
    
    // Standard F.7.4.3.1.1 (F-2) (JCTVC-R1008_v7)
    // The variable ScalabilityId[ i ][ smIdx ] specifying the identifier of the smIdx-th scalability dimension type of the i-th layer, and the variables ViewOrderIdx[ lId ], DependencyId[ lId ], and AuxId[ lId ] specifying the view order index, the spatial/quality scalability identifier, and the auxiliary identifier, respectively, of the layer with nuh_layer_id equal to lId  are derived as follows:
    int NumViews = 1;
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      int lId = layer_id_in_nuh[i];
      int_2d ScalabilityId;
      for (int smIdx = 0, j = 0; smIdx < 16; smIdx++)
      {
        if (scalability_mask_flag[smIdx])
          ScalabilityId[i][smIdx] = dimension_id[i][j++];
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
          if (ViewOrderIdx[ lId ] == ViewOrderIdx[ layer_id_in_nuh[j] ] )
            newViewFlag  = 0;
        }
        NumViews += newViewFlag ;
      }
    }
    
    for (int i = 0; i < NumViews; i++)
    {
      view_id_val[i] = get_bits(reader,vlc);
    }
  }

  for (int i = 1; i <= MaxLayersMinus1; i++) {
    for (int j = 0; j < i; j++)
    {
      direct_dependency_flag[i][j] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3.1.1 (F-3) (JCTVC-R1008_v7)
  // The variable DependencyFlag[ i ][ j ] is derived as follows:
  for( int i = 0; i <= MaxLayersMinus1; i++ ) {
    for( int j = 0; j <= MaxLayersMinus1; j++ ) {
      DependencyFlag[i][j] = direct_dependency_flag[i][j];
      for (int k = 0; k < i; k++) {
        if (direct_dependency_flag[i][k] && DependencyFlag[k][j])
          DependencyFlag[i][j] = true;
      }
    }
  }

  // Standard F.7.4.3.1.1 (F-4) (JCTVC-R1008_v7)
  // The variables NumDirectRefLayers[ iNuhLId ], IdDirectRefLayer[ iNuhLId ][ d ], NumRefLayers[ iNuhLId ], IdRefLayer[ iNuhLId ][ r ], NumPredictedLayers[ iNuhLId ], and IdPredictedLayer[ iNuhLId ][ p ] are derived as follows:
  for( int i = 0; i <= MaxLayersMinus1; i++ ) {
    int iNuhLId = layer_id_in_nuh[i];
    int d = 0; int r = 0; int p = 0;
    for( int j = 0; j <= MaxLayersMinus1; j++ ) {
      int jNuhLid = layer_id_in_nuh[ j ];
      if( direct_dependency_flag[i][j] )
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
    int iNuhLId = layer_id_in_nuh[i];
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
    num_add_layer_sets = get_uvlc(reader);
  }

  for( int i = 0; i < num_add_layer_sets; i++ ) {
    for( int j = 1; j < NumIndependentLayers; j++ ) {
      int nr_bits = ceil_log2(NumLayersInTreePartition[j] + 1);
      highest_layer_idx_plus1[ i ][ j ] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3 - layer_id_included_flag - (7-3) (JCTVC-R1013_v6)
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
  for( int i = 0; i < num_add_layer_sets; i++ ) {
    int layerNum = 0;
    int lsIdx = vps->vps_num_layer_sets + i;
    for (int treeIdx = 1; treeIdx < NumIndependentLayers; treeIdx++) {
      for( int layerCnt = 0; layerCnt < highest_layer_idx_plus1[ i ][ treeIdx ]; layerCnt++ ) {
        LayerSetLayerIdList[ lsIdx ][ layerNum++ ] = TreePartitionLayerIdList[ treeIdx ][ layerCnt ];
      }
    }
    NumLayersInIdList[ lsIdx ] = layerNum;
  }

  vps_sub_layers_max_minus1_present_flag = vlc = get_bits(reader,1);
  if (vlc) {
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      sub_layers_vps_max_minus1[i] = get_bits(reader,3);
    }
  }

  max_tid_ref_present_flag = vlc = get_bits(reader,1);
  if (vlc) {
    for (int i = 0; i <= MaxLayersMinus1; i++) {
      for( int j = i + 1; j <= MaxLayersMinus1; j++ ) {
        if (direct_dependency_flag[j][i])
          max_tid_il_ref_pics_plus1[i][j] = get_bits(reader,3);
      }
    }
  }

  default_ref_layers_active_flag = get_bits(reader,1);
  vps_num_profile_tier_level_minus1 = get_uvlc(reader);
  for( int i = vps->vps_base_layer_internal_flag ? 2 : 1;
             i <= vps_num_profile_tier_level_minus1; i++ ) {
    vps_profile_present_flag[i] = vlc = get_bits(reader,1);
    vps_ext_PTL[i].read_profile_tier_level(reader, vlc, vps->vps_max_sub_layers);
  }

  // Standard F.7.4.3.1.1 (F-6) (JCTVC-R1008_v7)
  // The variable NumLayerSets is derived as follows:
  NumLayerSets = vps->vps_num_layer_sets + num_add_layer_sets;
  if (NumLayerSets > 1) {
    num_add_olss = get_uvlc(reader);
    default_output_layer_idc = get_bits(reader,2);
  }

  output_layer_flag[0][0] = 1;
  OutputLayerFlag[0][0] = 1;

  NumOutputLayerSets = num_add_olss + NumLayerSets;
  layer_set_idx_for_ols_minus1[0] = -1;
  for( int i = 1; i < NumOutputLayerSets; i++ ) {
    if (NumLayerSets > 2 && i >= NumLayerSets) {
      int nr_bits = ceil_log2(NumLayerSets - 1);
      layer_set_idx_for_ols_minus1[i] = get_bits(reader,nr_bits);
    }
    else {
      layer_set_idx_for_ols_minus1[i] = i - 1;
    }
    int defaultOutputLayerIdc = libde265_min(default_output_layer_idc, 2);

    // Standard F.7.4.3.1.1 (F-10) (JCTVC-R1008_v7)
    // For i in the range of 0 to NumOutputLayerSets − 1, inclusive, the variable OlsIdxToLsIdx[ i ] is derived as specified in the following:
    OlsIdxToLsIdx[ i ] =  ( i < NumLayerSets ) ? i : ( layer_set_idx_for_ols_minus1[ i ] + 1 );

    if (i > vps->vps_num_layer_sets - 1 || defaultOutputLayerIdc == 2) {
      // i in the range of ( defaultOutputLayerIdc  = =  2 ) ? 0 : ( vps_num_layer_sets_minus1 + 1 ) to NumOutputLayerSets − 1, inclusive, 
      for (int j = 0; j < NumLayersInIdList[OlsIdxToLsIdx[i]]; j++) {
        output_layer_flag[ i ][ j ] = get_bits(reader,1);
        OutputLayerFlag[i][j] = output_layer_flag[ i ][ j ];
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
      if( NecessaryLayerFlag[ i ][ j ] && vps_num_profile_tier_level_minus1 > 0 ) {
        int nr_bits = ceil_log2(vps_num_profile_tier_level_minus1 + 1 );
        profile_tier_level_idx[ i ][ j ] = get_bits(reader,nr_bits);
      }
    }

    if (NumOutputLayersInOutputLayerSet[i] == 1 && NumDirectRefLayers[OlsHighestOutputLayerId[i]] > 0) {
      alt_output_layer_flag[ i ] = get_bits(reader,1);
    }
  }

  vps_num_rep_formats_minus1 = get_uvlc(reader);
  for (int i = 0; i <= vps_num_rep_formats_minus1; i++) {
    de265_error err = vps_ext_rep_format[i].parse_rep_format(reader);
    if ( err != DE265_OK) {
      return err;
    }
  }

  if (vps_num_rep_formats_minus1 > 0) {
    rep_format_idx_present_flag = get_bits(reader,1);
  }
  if (rep_format_idx_present_flag) {
    for( int i = vps->vps_base_layer_internal_flag ? 1 : 0; i <= MaxLayersMinus1; i++ ) {
      int nr_bits = ceil_log2( vps_num_rep_formats_minus1 + 1);
      vps_rep_format_idx[ i ] = get_bits(reader,nr_bits);
    }
  }
  max_one_active_ref_layer_flag = get_bits(reader,1);
  vps_poc_lsb_aligned_flag = get_bits(reader,1);

  for (int i = 1; i <= MaxLayersMinus1; i++) {
    if( NumDirectRefLayers[ layer_id_in_nuh[ i ] ] == 0 ) {
      poc_lsb_not_present_flag[ i ] = get_bits(reader,1);
    }
  }

  // Standard F.7.4.3.1.1 (F-9) (JCTVC-R1008_v7)
  // The variable MaxSubLayersInLayerSetMinus1[ i ] is derived as follows:
  MaxSubLayersInLayerSetMinus1;
  for( int i = 0; i < NumLayerSets; i++ ) {
    int maxSlMinus1 = 0;
    for( int k = 0; k < NumLayersInIdList[ i ]; k++ ) {
      int lId = LayerSetLayerIdList[ i ][ k ];
      maxSlMinus1 = libde265_max( maxSlMinus1, sub_layers_vps_max_minus1[ LayerIdxInVps[ lId ] ] );
    }
    MaxSubLayersInLayerSetMinus1[ i ] = maxSlMinus1;
  }

  // dpb_size() (F.7.3.2.1.3) (JCTVC-R1008_v7)
  dpb_size_table.read_decoded_picture_buffer_size_table(reader, vps);

  direct_dep_type_len_minus2 = get_uvlc(reader);
  direct_dependency_all_layers_flag = get_bits(reader,1);
  if (direct_dependency_all_layers_flag) {
    int nr_bits = direct_dep_type_len_minus2 + 2;
    direct_dependency_all_layers_type = get_bits(reader,nr_bits);
  }
  else {
    for (int i = vps->vps_base_layer_internal_flag ? 1 : 2; i <= MaxLayersMinus1; i++) {
      for( int j = vps->vps_base_layer_internal_flag ? 0 : 1; j < i; j++ ) {
        if( direct_dependency_flag[ i ][ j ] ) {
          int nr_bits = direct_dep_type_len_minus2 + 2;
          direct_dependency_type[ i ][ j ] = get_bits(reader,nr_bits);
        }
      }
    }
  }

  vps_non_vui_extension_length = get_uvlc(reader);
  for (int i = 1; i <= vps_non_vui_extension_length; i++) {
    int vps_non_vui_extension_data_byte = get_bits(reader, 8);
  }
	vps_vui_present_flag = get_bits(reader,1);

  if (vps_vui_present_flag) {
    // Byte alignment (vps_vui_alignment_bit_equal_to_one)
    for (int nrBits = bits_to_byte_boundary(reader); nrBits > 0; nrBits--) {
      if (get_bits(reader, 1) != 1) {
        return DE265_ERROR_CODED_PARAMETER_OUT_OF_RANGE;
      }
    }
    vui.read_vps_vui(reader, vps);
  }

  // Reading the VPS extension is done.
  // Check (calculate if not set) the output layer set index
  multilayer_decoder_parameters* ml_dec_param = ctx->get_multilayer_decoder_parameters();
  if (ml_dec_param) {
    if (!ml_dec_param->values_checked) {
      if (ml_dec_param->TargetOlsIdx == -1) {
        // Target output layer set ID not set yet
        if (ml_dec_param->TargetLayerId > vps->vps_max_layer_id) {
          // Target layer ID too high
          ml_dec_param->TargetLayerId = vps->vps_max_layer_id;
        }

        bool layerSetMatchFound = false;
        // Output layer set index not assigned.
        // Based on the value of targetLayerId, check if any of the output layer matches
        // Currently, the target layer ID in the encoder assumes that all the layers are decoded    
        // Check if any of the output layer sets match this description
        for(int i = 0; i < NumOutputLayerSets; i++)
        {
          bool layerSetMatchFlag = false;
          int layerSetIdx = layer_set_idx_for_ols_minus1[i] + 1;

          for(int j = 0; j < NumLayersInIdList[layerSetIdx]; j++)
          {
            if( LayerSetLayerIdList[layerSetIdx][j] == ml_dec_param->TargetLayerId )
            {
              layerSetMatchFlag = true;
              break;
            }
          }
      
          if( layerSetMatchFlag ) // Potential output layer set candidate found
          {
            // If target dec layer ID list is also included - check if they match
            if( !ml_dec_param->TargetDecLayerSetIdx.empty() )
            {
              for(int j = 0; j < NumLayersInIdList[layerSetIdx]; j++)
              {
                if (ml_dec_param->TargetDecLayerSetIdx[j] != layer_id_in_nuh[LayerSetLayerIdList[layerSetIdx][j]])
                {
                  layerSetMatchFlag = false;
                }
              }
            }
            if( layerSetMatchFlag ) // The target dec layer ID list also matches, if present
            {
              // Match found
              layerSetMatchFound = true;
              ml_dec_param->TargetOlsIdx = i;
              ml_dec_param->values_checked = true;
              break;
            }
          }
        }
        assert( layerSetMatchFound ); // No output layer set matched the value of either targetLayerId or targetdeclayerIdlist
  
      }
    }
    else {
      assert( ml_dec_param->TargetOlsIdx < NumOutputLayerSets );
      int layerSetIdx = layer_set_idx_for_ols_minus1[ ml_dec_param->TargetOlsIdx ] + 1;  // Index to the layer set
      // Check if the targetdeclayerIdlist matches the output layer set
      if( !ml_dec_param->TargetDecLayerSetIdx.empty() ) {
        for(int i = 0; i < NumLayersInIdList[layerSetIdx]; i++)
        {
          assert( ml_dec_param->TargetDecLayerSetIdx[i] == layer_id_in_nuh[LayerSetLayerIdList[layerSetIdx][i]]);
        }
      }
      ml_dec_param->values_checked = true;
    }
  }

  return DE265_OK;
}

de265_error decoded_picture_buffer_size_table::read_decoded_picture_buffer_size_table(bitreader* reader,
                                                                                      video_parameter_set *vps)
{
  video_parameter_set_extension *vps_ext = &vps->vps_extension;

  for( int i = 1; i < vps_ext->NumOutputLayerSets; i++ ) {
    int currLsIdx = vps_ext->OlsIdxToLsIdx[ i ];
    sub_layer_flag_info_present_flag[ i ] = get_bits(reader,1);
    for( int j = 0; j <= vps_ext->MaxSubLayersInLayerSetMinus1[ currLsIdx ]; j++ ) {
      if (j > 0 && sub_layer_flag_info_present_flag[i]) {
        sub_layer_dpb_info_present_flag[ i ][ j ] = get_bits(reader,1);
      }
      sub_layer_dpb_info_present_flag[i][0] = true;
      if( sub_layer_dpb_info_present_flag[ i ][ j ] ) {
        for (int k = 0; k < vps_ext->NumLayersInIdList[currLsIdx]; k++) {
          if( vps_ext->NecessaryLayerFlag[ i ][ k ]  &&  ( vps->vps_base_layer_internal_flag || 
            ( vps_ext->LayerSetLayerIdList[ currLsIdx ][ k ] != 0 ) ) ) {
            max_vps_dec_pic_buffering_minus1[i][k][j] = get_uvlc(reader);
          }
        }
        max_vps_num_reorder_pics[ i ][ j ] = get_uvlc(reader);
        max_vps_latency_increase_plus1[ i ][ j ] = get_uvlc(reader);
      }
    }
  }

  return DE265_OK;
}

de265_error profile_tier_level::read_profile_tier_level( bitreader* reader,
                                                         bool profile_present_flag,
                                                         int max_sub_layers)
{
  if (profile_present_flag) {
    general_profile_space = get_bits(reader,2);
    general_tier_flag = get_bits(reader,1);
    general_profile_idc = get_bits(reader,5);

    for (int i=0; i<32; i++) {
      general_profile_compatibility_flag[i] = get_bits(reader,1);
    }

    general_progressive_source_flag = get_bits(reader,1);
    general_interlaced_source_flag  = get_bits(reader,1);
    general_non_packed_constraint_flag = get_bits(reader,1);
    general_frame_only_constraint_flag = get_bits(reader,1);
    skip_bits(reader,44);
  }

  general_level_idc = get_bits(reader,8);


  for (int i=0; i<max_sub_layers-1; i++)
    {
      profile[i].sub_layer_profile_present_flag = get_bits(reader,1);
      profile[i].sub_layer_level_present_flag   = get_bits(reader,1);
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
      if (profile[i].sub_layer_profile_present_flag)
        {
          profile[i].sub_layer_profile_space = get_bits(reader,2);
          profile[i].sub_layer_tier_flag = get_bits(reader,1);
          profile[i].sub_layer_profile_idc = get_bits(reader,5);

          for (int j=0; j<32; j++)
            {
              profile[i].sub_layer_profile_compatibility_flag[j] = get_bits(reader,1);
            }

          profile[i].sub_layer_progressive_source_flag = get_bits(reader,1);
          profile[i].sub_layer_interlaced_source_flag  = get_bits(reader,1);
          profile[i].sub_layer_non_packed_constraint_flag = get_bits(reader,1);
          profile[i].sub_layer_frame_only_constraint_flag = get_bits(reader,1);
          skip_bits(reader,44);
        }

      if (profile[i].sub_layer_level_present_flag)
        {
          profile[i].sub_layer_level_idc = get_bits(reader,8);
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

void video_parameter_set::dump_vps(int fd)
{
  FILE* fh;
  if (fd==1) fh=stdout;
  else if (fd==2) fh=stderr;
  else { return; }

  LOG0("----------------- VPS -----------------\n");
  LOG1("video_parameter_set_id                : %d\n", video_parameter_set_id);
  LOG1("vps_base_layer_internal_flag          : %d\n", vps_base_layer_internal_flag);
  LOG1("vps_base_layer_available_flag         : %d\n", vps_base_layer_available_flag);
  LOG1("vps_max_layers                        : %d\n", vps_max_layers);
  LOG1("vps_max_sub_layers                    : %d\n", vps_max_sub_layers);
  LOG1("vps_temporal_id_nesting_flag          : %d\n", vps_temporal_id_nesting_flag);

  profile_tier_level.dump_profile_tier_level(vps_max_sub_layers, fh);
  //dump_bit_rate_pic_rate_info(&bit_rate_pic_rate_info, 0, vps_max_sub_layers-1);

  LOG1("vps_sub_layer_ordering_info_present_flag : %d\n",
       vps_sub_layer_ordering_info_present_flag);

  if (vps_sub_layer_ordering_info_present_flag) {
    for (int i=0;i<vps_max_sub_layers;i++) {
      LOG2("layer %d: vps_max_dec_pic_buffering = %d\n",i,layer[i].vps_max_dec_pic_buffering);
      LOG1("         vps_max_num_reorder_pics  = %d\n",layer[i].vps_max_num_reorder_pics);
      LOG1("         vps_max_latency_increase  = %d\n",layer[i].vps_max_latency_increase);
    }
  }
  else {
    LOG1("layer (all): vps_max_dec_pic_buffering = %d\n",layer[0].vps_max_dec_pic_buffering);
    LOG1("             vps_max_num_reorder_pics  = %d\n",layer[0].vps_max_num_reorder_pics);
    LOG1("             vps_max_latency_increase  = %d\n",layer[0].vps_max_latency_increase);
  }


  LOG1("vps_max_layer_id   = %d\n", vps_max_layer_id);
  LOG1("vps_num_layer_sets = %d\n", vps_num_layer_sets);

  for (int i=1; i <= vps_num_layer_sets-1; i++)
    for (int j=0; j <= vps_max_layer_id; j++)
      {
        LOG3("layer_id_included_flag[%d][%d] = %d\n",i,j,
             layer_id_included_flag[i][j]);
      }

  LOG1("vps_timing_info_present_flag = %d\n",
       vps_timing_info_present_flag);

  if (vps_timing_info_present_flag) {
    LOG1("vps_num_units_in_tick = %d\n", vps_num_units_in_tick);
    LOG1("vps_time_scale        = %d\n", vps_time_scale);
    LOG1("vps_poc_proportional_to_timing_flag = %d\n", vps_poc_proportional_to_timing_flag);

    if (vps_poc_proportional_to_timing_flag) {
      LOG1("vps_num_ticks_poc_diff_one = %d\n", vps_num_ticks_poc_diff_one);
      LOG1("vps_num_hrd_parameters     = %d\n", vps_num_hrd_parameters);

      for (int i=0; i<vps_num_hrd_parameters; i++) {
        LOG2("hrd_layer_set_idx[%d] = %d\n", i, hrd_layer_set_idx[i]);

        if (i > 0) {
          LOG2("cprms_present_flag[%d] = %d\n", i, cprms_present_flag[i]);
        }

        //hrd_parameters(cprms_present_flag[i], vps_max_sub_layers_minus1)

        return; // TODO: decode hrd_parameters()
      }
    }
  }

  LOG1("vps_extension_flag = %d\n", vps_extension_flag);
}


void profile_tier_level::dump_profile_tier_level(int max_sub_layers, FILE* fh)
{
  LOG1("  general_profile_space     : %d\n", general_profile_space);
  LOG1("  general_tier_flag         : %d\n", general_tier_flag);
  LOG1("  general_profile_idc       : %d\n", general_profile_idc);

  LOG0("  general_profile_compatibility_flags: ");
  for (int i=0; i<32; i++) {
    if (i) LOG0("*,");
    LOG1("*%d",general_profile_compatibility_flag[i]);
  }
  LOG0("*\n");

  LOG1("  general_level_idc         : %d\n", general_level_idc);

  for (int i=0; i<max_sub_layers-1; i++)
    {
      LOG1("  Profile/Tier/Level [Layer %d]\n",i);

      if (profile[i].sub_layer_profile_present_flag) {

        LOG1("    sub_layer_profile_space : %d\n",profile[i].sub_layer_profile_space);
        LOG1("    sub_layer_tier_flag     : %d\n",profile[i].sub_layer_tier_flag);
        LOG1("    sub_layer_profile_idc   : %d\n",profile[i].sub_layer_profile_idc);

        LOG0("    sub_layer_profile_compatibility_flags: ");
        for (int j=0; j<32; j++) {
          if (j) LOG0(",");
          LOG1("%d",profile[i].sub_layer_profile_compatibility_flag[j]);
        }
        LOG0("\n");

        LOG1("    sub_layer_progressive_source_flag : %d\n",profile[i].sub_layer_progressive_source_flag);
        LOG1("    sub_layer_interlaced_source_flag : %d\n",profile[i].sub_layer_interlaced_source_flag);
        LOG1("    sub_layer_non_packed_constraint_flag : %d\n",profile[i].sub_layer_non_packed_constraint_flag);
        LOG1("    sub_layer_frame_only_constraint_flag : %d\n",profile[i].sub_layer_frame_only_constraint_flag);
      }


      if (profile[i].sub_layer_level_present_flag) {
        LOG1("    sub_layer_level_idc   : %d\n", profile[i].sub_layer_level_idc);
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
