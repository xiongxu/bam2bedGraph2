#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "sds.h"

static inline sds int2sds(uint32_t pos);
static inline uint32_t sds2int(const sds string);

sds int2sds(uint32_t pos){
	short codeLen=1;
	if (pos>16777215){
		codeLen=4;
	}else if(pos>65535){
		codeLen=3;
	}else if(pos>255){
		codeLen=2;
	}
	sds value=sdsnewlen(NULL,codeLen);
	short i;
	for (i=codeLen-1;i>=0 ;i-- )
		value[i]= pos>>(i*8) & 0xFF;
	return value;
}

uint32_t sds2int(const sds string){
	uint32_t pos=0;
	int8_t i;
	for (i=sdslen(string)-1; i>=0 ;i-- ){
		pos=pos << 8 | (unsigned char)string[i];
	}
	return pos;
}