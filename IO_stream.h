#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <err.h>
#include <fcntl.h>
#include <unistd.h>
/*
typedef pFILE FILE *;
#define OPEN_STREAM(TYPE,fileno,mode)	\
TYPE open_stream_##TYPE(char* filename,void(*f)(int , const char *)) {		\
	int fd ;																\
	if (strncmp(filename,"-", 1)==0) {										\
		fd = fileno;														\
	} else {																\
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );			\
		if (fd==-1) err(1, "Failed to create output file (%s)", filename);	\
	}																		\
	TYPE out = f(fd,mode);													\
	return out;																\
}

#define CREATE_OUTFILE(TYPE)	\
TYPE creat_outfile_##TYPE(char *outfile,char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	if (!strncmp(outfile,"-",1)){
		TYPE fo=open_stream_##TYPE(out,fdopen);
	}
	else{
		TYPE fo=open_stream_##TYPE(out,gzdopen);
	}
	return fo;
}

OPEN_STREAM(pFILE,STDOUT_FILENO,"wb")
OPEN_STREAM(gzFile,STDOUT_FILENO,"wb")

CREATE_OUTFILE(gzFile)
CREATE_OUTFILE(pFILE)
*/

FILE *fopen_output_stream(const char *filename) {
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDOUT_FILENO;
	} else {
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );
		if (fd==-1) err(1, "Failed to create output file (%s)", filename);
	}
	FILE *out = fdopen(fd,"wb");
	return out;
}

FILE *fcreat_outfile(const char *outfile,const char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	FILE *fo=fopen_output_stream(out);
	return fo;
}

gzFile open_output_stream(char* filename) {
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDOUT_FILENO;
	} else {
		fd = open(filename, O_CREAT | O_WRONLY | O_TRUNC, 0666 );
		if (fd==-1) err(1, "Failed to create output file (%s)", filename);
	}
	gzFile out = gzdopen(fd,"wb");
	return out;
}

gzFile creat_outfile(const char *outfile,const char *suffix) {
	char *out=(char *)calloc(strlen(outfile)+strlen(suffix)+2,sizeof(char));
	sprintf(out,"%s%s",outfile,suffix);
	gzFile fo=open_output_stream(out);
	return fo;
}

gzFile open_input_stream(const char *filename){
	int fd ;
	if (strncmp(filename,"-", 1)==0 || !strcmp(filename,"")) {
		fd = STDIN_FILENO;
	} else {
		fd = open(filename, O_CREAT | O_RDONLY , 0666 );
		if (fd==-1) err(1, "Failed to create input file (%s)", filename);
	}
	gzFile in = gzdopen(fd,"rb");
	return in;
}
