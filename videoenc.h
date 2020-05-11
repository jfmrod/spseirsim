#ifndef VIDEOENC_H
#define VIDEOENC_H

#include <eutils/estr.h>

//int videoOpen();
int videoOpen(const estr& fname,int width,int height,int fps=30,int bitrate=2000);
void videoClose();
void videoPushFrame(uint8_t *data);

#endif

