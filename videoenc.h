#ifndef VIDEOENC_H
#define VIDEOENC_H

int videoOpen();
void videoClose();
void videoPushFrame(uint8_t *data);

#endif

