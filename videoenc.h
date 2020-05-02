#ifndef VIDEOENC_H
#define VIDEOENC_H

//int videoOpen();
int videoOpen(int width,int height,int fps=30,int bitrate=2000);
void videoClose();
void videoPushFrame(uint8_t *data);

#endif

