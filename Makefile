CC=			gcc
CFLAGS=		-O3 -Wall -W -Wstrict-prototypes -Wwrite-strings -g -ggdb -fPIC 
LDFLAGS=	-lhiredis -lz -lpthread -lm -lbam
INCLUDES=	-I$(CURDIR)/hiredis -I$(CURDIR)/zlib-1.2.11/include -I$(CURDIR)/samtools-0.1.19
LIBPATH=    -L$(CURDIR)/hiredis -L$(CURDIR)/zlib-1.2.11/lib -L$(CURDIR)/samtools-0.1.19

.PHONY:all clean
all:bam2bedGraph
#SUBDIRS = `find -maxdepth 1 -type d | sed "1d"`
$(CURDIR)/hiredis/libhiredis.a:
	wdir=`pwd`; \
	cd  $(CURDIR)/hiredis;\
	$(MAKE)  ;\
	cd $$wdir

$(CURDIR)/samtools-0.1.19/libbam.a:
	wdir=`pwd`; \
	cd  $(CURDIR)/samtools-0.1.19;\
	$(MAKE)  INCLUDES="-I. -I$(CURDIR)/zlib-1.2.11/include";\
	cd $$wdir

$(CURDIR)/zlib-1.2.11:
	wget http://www.zlib.net/zlib-1.2.11.tar.gz; \
	tar -zxvf zlib-1.2.11.tar.gz && mv zlib-1.2.11 zlib;\
	cd zlib;\
	test -d $(CURDIR)/zlib-1.2.11 || mkdir -p $(CURDIR)/zlib-1.2.11;\
	./configure --prefix=$(CURDIR)/zlib-1.2.11;\
	make;\
	make install;\
	cd $(CURDIR) && rm -rf zlib

bam2bedGraph:bam2bedGraph.o
	$(CC) $(CFLAGS) $< -o $@ $(LIBPATH) $(LDFLAGS) 

bam2bedGraph.o:bam2bedGraph.c $(CURDIR)/zlib-1.2.11 $(CURDIR)/hiredis/libhiredis.a $(CURDIR)/samtools-0.1.19/libbam.a
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDES)

clean:
	rm -rf zlib-1.2.11 zlib-1.2.11.tar.gz bam2bedGraph  *.o;\
	wdir=`pwd`; \
	cd ./hiredis;\
	$(MAKE) clean;\
	cd $$wdir; \
	cd  $(CURDIR)/samtools-0.1.19;\
	$(MAKE) clean ;\
	cd $$wdir
