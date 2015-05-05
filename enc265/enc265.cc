/*
 * H.265 video codec.
 * Copyright (c) 2013-2014 struktur AG, Dirk Farin <farin@struktur.de>
 *
 * This file is part of libde265.
 *
 * libde265 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libde265 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with libde265.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "libde265/en265.h" //coder-context.h"

//#include "libde265/image-io.h"
//#include "libde265/encoder/analyze.h"
#include "libde265/util.h"

#include <getopt.h>
#include <stdlib.h>


#if HAVE_VIDEOGFX
#include <libvideogfx.hh>
using namespace videogfx;

void debug_show_image_libvideogfx(const de265_image* input, int slot)
{
    static X11Win debugwin;
    static bool opened=false;
    int w = input->get_width();
    int h = input->get_height();
    if (!opened) {
      opened=true;
      debugwin.Create(w,h, "debug");
    }

    Image<Pixel> img;
    img.Create(w,h,Colorspace_YUV, Chroma_420);

    for (int y=0;y<h;y++)
      memcpy(img.AskFrameY()[y], input->get_image_plane_at_pos(0,0,y), w);

    for (int y=0;y<h/2;y++) {
      memcpy(img.AskFrameU()[y], input->get_image_plane_at_pos(1,0,y), w/2);
      memcpy(img.AskFrameV()[y], input->get_image_plane_at_pos(2,0,y), w/2);
    }

    debugwin.Display(img);
    //debugwin.WaitForKeypress();
}
#endif



int show_help=false;
int verbosity=0;
const char* output_bitstream = "out.str";
const char* input_YUV = NULL;
const char* rec_YUV = NULL;
int inp_width = -1;
int inp_height = -1;
int frame_start = 0;
int frame_num = -1;

static struct option long_options[] = {
  {"help",            no_argument,       &show_help, 1 },
  {"verbose",         no_argument,       0, 'v' },
  {"bitstream",       required_argument, 0, 'b' },
  {"input",           required_argument, 0, 'i' },
  {"reconstruction",  required_argument, 0, 'r' },
  {"width",           required_argument, 0, 'w' },
  {"height",          required_argument, 0, 'h' },
  {"frame_start",     required_argument, 0, 's' },
  {"frame_num",       required_argument, 0, 'n' },
  {0,            0,                 0,  0 }
};

void test_parameters_API(en265_encoder_context* ectx)
{
  const char** param = en265_list_parameters(ectx);
  if (param) {
    for (int i=0; param[i]; i++) {
      printf("|%s| ",param[i]);

      enum en265_parameter_type type = en265_get_parameter_type(ectx, param[i]);
      const char* type_name="unknown";
      switch (type) {
      case en265_parameter_int: type_name="int"; break;
      case en265_parameter_bool: type_name="bool"; break;
      case en265_parameter_string: type_name="string"; break;
      case en265_parameter_choice: type_name="choice"; break;
      }

      printf("(%s)",type_name);

      if (type==en265_parameter_choice) {
        const char** choices = en265_list_parameter_choices(ectx, param[i]);
        if (choices) {
          for (int k=0; choices[k]; k++) {
            printf(" %s",choices[k]);
          }
        }
      }

      printf("\n");
    }
  }

  // en265_set_parameter_int(ectx, "min-tb-size", 8);
}


extern int skipTBSplit, noskipTBSplit;
extern int zeroBlockCorrelation[6][2][5];

int main(int argc, char** argv)
{
  de265_init();

  en265_encoder_context* ectx = en265_new_encoder();

  bool cmdline_errors = false;

  // --- in/out parameters ---



  // --- read encoder parameters ---

  if (en265_parse_command_line_parameters(ectx, &argc, argv) != DE265_OK) {
    cmdline_errors = true;
  }



  while (1) {
    int option_index = 0;

    int c = getopt_long(argc, argv, "vb:i:r:w:h:s:n:"
                        , long_options, &option_index);
    if (c == -1)
      break;

    switch (c) {
      case 'v': verbosity++; break;
      case 'b' : output_bitstream = optarg; break;
      case 'i' : input_YUV = optarg; break;
      case 'r' : rec_YUV = optarg; break;
      case 'w' : inp_width = atoi(optarg); break;
      case 'h' : inp_height = atoi(optarg); break;
      case 's' : frame_start = atoi(optarg); break;
      case 'n' : frame_num = atoi(optarg); break;
    }
  }

  if (!show_help) {
    // These options are required
    if (input_YUV == NULL) {
      fprintf(stderr,"No input YUV file provided.\n");
      cmdline_errors = true;
    }
    if (inp_width == -1) {
      fprintf(stderr,"Width of the input YUV file not prvided.\n");
      cmdline_errors = true;
    }
    if (inp_height == -1) {
      fprintf(stderr,"Height of the input YUV file not prvided.\n");
      cmdline_errors = true;
    }
  }

  // --- show usage information ---

  if (optind != argc || cmdline_errors || show_help) {
    fprintf(stderr," enc265  v%s\n", de265_get_version());
    fprintf(stderr,"--------------\n");
    fprintf(stderr,"usage: enc265 [options] -i inp.yuv -w 352 -h 288\n");
    fprintf(stderr,"The video file must be a raw YUV file\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"options:\n");
    fprintf(stderr,"      --help              show help\n");
    fprintf(stderr,"  -v, --verbose           increase verbosity level (up to 3 times)\n");
    fprintf(stderr,"  -b, --bitstream         the output bitstream file (default: out.str)\n");
    fprintf(stderr,"  -i, --input             the raw input YUV file\n");
    fprintf(stderr,"  -r, --reconstruction    the reconstruction YUV file\n");
    fprintf(stderr,"  -w, --width             the width of the input YUV file\n");
    fprintf(stderr,"  -h, --height            the height of the input YUV file\n");
    fprintf(stderr,"  -s, --frame_start       the frame number of frames to skip in the input YUV file\n");
    fprintf(stderr,"  -n, --frame_num         the number of frames to encode\n");

    //inout_param_config.print_params();
    fprintf(stderr,"\n");
    en265_show_parameters(ectx);

    exit(show_help ? 0 : 5);
  }

  de265_set_verbosity(verbosity);
#if HAVE_VIDEOGFX
  //debug_set_image_output(debug_show_image_libvideogfx);
#endif

  //test_parameters_API(ectx);


  //ImageSink_YUV reconstruction_sink;
  //if (strlen(inout_params.reconstruction_yuv.get().c_str()) != 0) {
  //  reconstruction_sink.set_filename(inout_params.reconstruction_yuv.get().c_str());
  //  //ectx.reconstruction_sink = &reconstruction_sink;
  //}

  //ImageSource_YUV image_source;
  //image_source.set_input_file(inout_params.input_yuv.get().c_str(),
  //                            inout_params.input_width,
  //                            inout_params.input_height);

  //PacketSink_File packet_sink;
  //packet_sink.set_filename(inout_params.output_filename.get().c_str());


  // --- run encoder ---

  //image_source.skip_frames( inout_params.first_frame );

  //en265_start_encoder(ectx, 0);

  //int maxPoc = INT_MAX;
  //if (inout_params.max_number_of_frames.is_defined()) {
  //  maxPoc = inout_params.max_number_of_frames;
  //}

  //bool eof = false;
  //for (int poc=0; poc<maxPoc && !eof ;poc++)
  //  {
  //    // push one image into the encoder

  //    de265_image* input_image = image_source.get_image();
  //    if (input_image==NULL) {
  //      en265_push_eof(ectx);
  //      eof=true;
  //    }
  //    else {
  //      en265_push_image(ectx, input_image);
  //    }



  //    // encode images while more are available

  //    en265_encode(ectx);


  //    // write all pending packets

  //    for (;;) {
  //      en265_packet* pck = en265_get_packet(ectx,0);
  //      if (pck==NULL)
  //        break;

  //      packet_sink.send_packet(pck->data, pck->length);

  //      en265_free_packet(ectx,pck);
  //    }
  //  }



  // --- print statistics ---

  //en265_print_logging((encoder_context*)ectx, "tb-split", NULL);


  en265_free_encoder(ectx);

  de265_free();

  return 0;
}
