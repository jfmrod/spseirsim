#include <iostream>

extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/avutil.h>
#include <libavutil/time.h>
#include <libavutil/opt.h>
#include <libswscale/swscale.h>
}



AVFrame *videoFrame = nullptr;
AVFrame *audioFrame = nullptr;
AVCodecContext *cctx = nullptr;
SwsContext *swsCtx = nullptr;
AVCodecContext *acctx = nullptr;
int frameCounter = 0;
AVFormatContext *ofctx = nullptr;
AVOutputFormat *oformat = nullptr;

int fps = 5;
//int width = 1920;
//int height = 1080;
//int bitrate = 2000;

void audioPushFrame(uint8_t *data){
  int err;
  if (!audioFrame) {
    audioFrame = av_frame_alloc();
    audioFrame->format = AV_SAMPLE_FMT_U8;
//    audioFrame->nChannels = 1;
//    audioFrame->nSamplesPerSec = 22000;
    
//    videoFrame->width = cctx->width;
//    videoFrame->height = cctx->height;
    if ((err = av_frame_get_buffer(audioFrame, 32)) < 0) {
      std::cout <<  "Failed to allocate picture" << err << std::endl;
      return;
    }
  }

//  sws_scale(swsCtx, (const uint8_t * const *)&data, inLinesize, 0, cctx->height, videoFrame->data, videoFrame->linesize);

  audioFrame->pts = (1.0/fps)*90000.0*(frameCounter-1);
//  std::cout << videoFrame->pts <<" " << cctx->time_base.num << " " << cctx->time_base.den << " " << frameCounter<< std::endl;
  if ((err = avcodec_send_frame(acctx, audioFrame)) < 0) {
    std::cout << "Failed to send audio frame" << err <<std::endl;
    return;
  }

  AVPacket pkt;
  av_init_packet(&pkt);
  pkt.data = NULL;
  pkt.size = 0;
  pkt.stream_index = 0;
  pkt.flags |= AV_PKT_FLAG_KEY;
  if (avcodec_receive_packet(acctx, &pkt) == 0) {
    static int counter = 0;
//    std::cout << "pkt key: " << (pkt.flags & AV_PKT_FLAG_KEY) <<" " << pkt.size << " " << (counter++) << std::endl;
    uint8_t *size = ((uint8_t*)pkt.data);
//    std::cout << "first: " << (int)size[0] << " " << (int)size[1] << " " << (int)size[2] << " " << (int)size[3] <<" "  << (int)size[4] << " " << (int)size[5] << " " << (int)size[6] << " " << (int)size[7] << std::endl;
    av_interleaved_write_frame(ofctx, &pkt);
    av_packet_unref(&pkt);
  }

}

void videoPushFrame(uint8_t *data){
  int err;
  if (!videoFrame) {
    videoFrame = av_frame_alloc();
    videoFrame->format = AV_PIX_FMT_YUV420P;
    videoFrame->width = cctx->width;
    videoFrame->height = cctx->height;
    if ((err = av_frame_get_buffer(videoFrame, 32)) < 0) {
      std::cout <<  "Failed to allocate picture" << err << std::endl;
      return;
    }
  }
  if (!swsCtx)
    swsCtx = sws_getContext(cctx->width, cctx->height, AV_PIX_FMT_RGBA, cctx->width, cctx->height, AV_PIX_FMT_YUV420P, SWS_BICUBIC, 0, 0, 0);

  int inLinesize[1] = { 4 * cctx->width };

  // From RGB to YUV
  sws_scale(swsCtx, (const uint8_t * const *)&data, inLinesize, 0, cctx->height, videoFrame->data, videoFrame->linesize);
  videoFrame->pts = (1.0/fps)*90000.0*(frameCounter++);
//  std::cout << videoFrame->pts <<" " << cctx->time_base.num << " " << cctx->time_base.den << " " << frameCounter<< std::endl;
  if ((err = avcodec_send_frame(cctx, videoFrame)) < 0) {
    std::cout << "Failed to send frame" << err <<std::endl;
    return;
  }
//  AV_TIME_BASE;

  AVPacket pkt;
  av_init_packet(&pkt);
  pkt.data = NULL;
  pkt.size = 0;
  pkt.flags |= AV_PKT_FLAG_KEY;
  if (avcodec_receive_packet(cctx, &pkt) == 0) {
    static int counter = 0;
    if (counter == 0){
      FILE *fp = fopen("dump_first_frame1.dat", "wb");
      fwrite(pkt.data, pkt.size,1,fp);
      fclose(fp);
    }
//    std::cout << "pkt key: " << (pkt.flags & AV_PKT_FLAG_KEY) <<" " << pkt.size << " " << (counter++) << std::endl;
    uint8_t *size = ((uint8_t*)pkt.data);
//    std::cout << "first: " << (int)size[0] << " " << (int)size[1] << " " << (int)size[2] << " " << (int)size[3] <<" "  << (int)size[4] << " " << (int)size[5] << " " << (int)size[6] << " " << (int)size[7] << std::endl;
    av_interleaved_write_frame(ofctx, &pkt);
    av_packet_unref(&pkt);
  }
}

static void finish() {
  //DELAYED FRAMES
  AVPacket pkt;
  av_init_packet(&pkt);
  pkt.data = NULL;
  pkt.size = 0;

  for (;;) {
    avcodec_send_frame(cctx, NULL);
    if (avcodec_receive_packet(cctx, &pkt) == 0) {
        av_interleaved_write_frame(ofctx, &pkt);        
        av_packet_unref(&pkt);
    }
    else {
        break;
    }
  }

  av_write_trailer(ofctx);
  if (!(oformat->flags & AVFMT_NOFILE)) {
    int err = avio_close(ofctx->pb);
    if (err < 0) {
        std::cout << "Failed to close file" << err <<std::endl;
    }
  }
}

static void free(){
  if (videoFrame) {
    av_frame_free(&videoFrame);
  }
  if (cctx) {
    avcodec_free_context(&cctx);
  }
  if (ofctx) {
    avformat_free_context(ofctx);
  }
  if (swsCtx) {
    sws_freeContext(swsCtx);
  }
}

int videoOpen(int width,int height,int _fps,int bitrate)
{
  fps=_fps;
  av_register_all();
  avcodec_register_all();

  oformat = av_guess_format(nullptr, "test.mp4", nullptr);
  if (!oformat) {
    std::cout << "can't create output format" << std::endl;
    return -1;
  }
//  oformat->video_codec = AV_CODEC_ID_H265;

  int err = avformat_alloc_output_context2(&ofctx, oformat, nullptr, "test.mp4");
  if (err){
    std::cout << "can't create output context" << std::endl;
    return -1;
  }

  AVCodec *codec = avcodec_find_encoder(oformat->video_codec);
  if (!codec) {
    std::cout << "can't create codec" << std::endl;
    return -1;
  }

  AVStream *stream = avformat_new_stream(ofctx, codec);
  if (!stream) {
    std::cout << "can't find format" << std::endl;
    return -1;
  }

  cctx = avcodec_alloc_context3(codec);
  if (!cctx)  {
    std::cout << "can't create codec context" << std::endl;
    return -1;
  }

  stream->codecpar->codec_id = oformat->video_codec;
  stream->codecpar->codec_type = AVMEDIA_TYPE_VIDEO;
  stream->codecpar->width = width;
  stream->codecpar->height = height;
//  stream->codecpar->aspect_ratio = {width,height};
  stream->codecpar->format = AV_PIX_FMT_YUV420P;
  stream->codecpar->bit_rate = bitrate * 1000;
  avcodec_parameters_to_context(cctx, stream->codecpar);
//  cctx->time_base = (AVRational){ 1, 1 };
//  cctx->time_base = (AVRational){ 1, fps };

  cctx->width = width;
  cctx->height = height;
//  cctx->aspect_ratio = {width,height};
  cctx->pix_fmt = AV_PIX_FMT_YUV420P;


  cctx->pix_fmt = AV_PIX_FMT_YUV420P;
  cctx->time_base = (AVRational){ 1000, 1 };
  cctx->max_b_frames = 0;
//  cctx->max_b_frames = 2;
  cctx->gop_size = 12;
//  cctx->gop_size = fps;
  cctx->framerate= (AVRational){fps, 1};

//  cctx->codec_tag = MKTAG('h', 'v', 'c', '1');
/*
  cctx->bit_rate_tolerance = 0;
  cctx->rc_max_rate = 0;
  cctx->rc_buffer_size = 0;
  cctx->gop_size = 40;
  cctx->max_b_frames = 3;
  cctx->me_cmp = 1;
  cctx->me_range = 16;
  cctx->qmin = 10;
  cctx->qmax = 51;
  cctx->flags |= CODEC_FLAG_LOOP_FILTER;
  cctx->me_subpel_quality = 5;
  cctx->i_quant_factor = 0.71;
  cctx->qcompress = 0.6;
  cctx->max_qdiff = 4;

  cctx->flags2 |= CODEC_FLAG2_FAST;
//  cctx->flags2 |= CODEC_FLAG2_FASTPSKIP;
*/
//must remove the following
//  cctx->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

//  cctx->level=30;
  cctx->level=32;
  av_opt_set(cctx->priv_data, "profile", "baseline", AV_OPT_SEARCH_CHILDREN);

/*
  if (stream->codecpar->codec_id == AV_CODEC_ID_H264)
    av_opt_set(cctx, "preset", "ultrafast", 0);
  else if (stream->codecpar->codec_id == AV_CODEC_ID_H265)
    av_opt_set(cctx, "preset", "ultrafast", 0);
*/
//  if (oformat->flags & AVFMT_GLOBALHEADER)
//    cctx->flags |= CODEC_FLAG_GLOBAL_HEADER;

  avcodec_parameters_from_context(stream->codecpar, cctx);
  if ((err = avcodec_open2(cctx, codec, NULL)) < 0) {
    std::cout << "Failed to open codec" << err << std::endl;
    return -1;
  }


/*
  // audio codec

  AVCodec *acodec = avcodec_find_encoder(oformat->audio_codec);
  if (!acodec) {
    std::cout << "can't create audio codec" << std::endl;
    return -1;
  }

  AVStream *astream = avformat_new_stream(ofctx, acodec);
  if (!stream) {
    std::cout << "can't find audio format" << std::endl;
    return -1;
  }

  acctx = avcodec_alloc_context3(acodec);
  if (!acctx)  {
    std::cout << "can't create audio codec context" << std::endl;
    return -1;
  }

  astream->codecpar->codec_id = oformat->audio_codec;
  astream->codecpar->codec_type = AVMEDIA_TYPE_AUDIO;
  avcodec_parameters_to_context(acctx, astream->codecpar);

  acctx->time_base = (AVRational){ 1000, 1 };
  acctx->framerate= (AVRational){fps, 1};


  acctx->sample_rate = 22050;
  acctx->channel_layout = 1;
  acctx->channels = av_get_channel_layout_nb_channels(acctx->channel_layout);
  // take first format from list of supported formats
//  acctx->sample_fmt = acodec->sample_fmts[0];
  acctx->sample_fmt = AV_SAMPLE_FMT_U8;
  acctx->time_base = (AVRational){1, acctx->sample_rate};

  avcodec_parameters_from_context(astream->codecpar, acctx);
  if ((err = avcodec_open2(acctx, acodec, NULL)) < 0) {
    std::cout << "Failed to open audio codec" << err << std::endl;
    return -1;
  }
*/

  if (!(oformat->flags & AVFMT_NOFILE)) {
    if ((err = avio_open(&ofctx->pb, "test.mp4", AVIO_FLAG_WRITE)) < 0) {
        std::cout << "Failed to open file" << err << std::endl;
        return -1;
    }
  }

  AVDictionary *fmtOptions = nullptr;
  av_dict_set(&fmtOptions, "movflags", "faststart", 0);
  av_dict_set(&fmtOptions, "brand", "mp42", 0);

//  if ((err = avformat_write_header(ofctx, NULL)) < 0) {
  if ((err = avformat_write_header(ofctx, &fmtOptions)) < 0) {
    std::cout << "Failed to write header" << err << std::endl;
    return -1;
  }

//  av_dump_format(ofctx, 0, "test.mp4", 1);
  return 0;
}

/*
  uint8_t *frameraw = new uint8_t[1920*1080*4];
  memset(frameraw, 222, 1920*1080*4);
  for(int i =0;i<180;++i)
    pushFrame(frameraw);
*/

void videoClose(){
//  delete [] frameraw;
  finish();
  free();
}
