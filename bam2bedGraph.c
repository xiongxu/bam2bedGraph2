// gcc -g -O3 -Wall bam2wig.c  -o bam2wig -I./samtools-0.1.19 -I./zlib_1.2.8/include -L./samtools-0.1.19 -L./zlib_1.2.8/lib -lpthread -lm -lz -lbam -lhiredis
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <inttypes.h>
#include <sys/time.h>
#include "IO_stream.h"
#include "dict.h"
#include "sds.h"
#include "int2sds.h"
#include "sam.h"

#define initdictType(mydictType)	\
do{													\
	(mydictType)->hashFunction=myhashFunction;		\
	(mydictType)->keyDup=mykeyDup;					\
	(mydictType)->valDup=NULL;						\
	(mydictType)->keyCompare=mykeyCompare;			\
	(mydictType)->keyDestructor=mykeyDestructor;	\
	(mydictType)->valDestructor=myvalDestructor;	\
}while (0)

typedef struct _data_buf_ {
	samfile_t *fp;
	dict *Start;
	dict *End;
	uint32_t *chrRC;
	uint64_t *chrReadLen;
}DataBuf;

struct globalArgs_t {
	char **infiles;
	short numInfiles;
	const char *outfile;
	const char *region;
	int strand;
} globalArgs;

static inline long long usec(void);
void display_usage(char * argv[]);
void region_sam2wig(DataBuf *databuf ,bam_index_t *idx,const char *region,FILE *bedGraph,dictType *mydictType);
void hash2BedGraph(DataBuf *databuf,int ref,FILE *stream);
sds *union_hashed_keys(dict *hashtbl_a,dict *hashtbl_b);
int bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf);
bam_index_t *load_bam_index(char *filename);
static inline int fetch_func(const bam1_t *b, void *data) ;
static inline int insertInt2Hash(dict *Start,uint32_t pos);

static inline void *mykeyDup(void *privdata, const void *key);
static inline int mykeyCompare(void *privdata, const void *key1, const void *key2);
static inline void mykeyDestructor(void *privdata, void *key);
static inline void myvalDestructor(void *privdata, void *obj);
static inline unsigned int myhashFunction(const void *key);

unsigned int myhashFunction(const void *key){
	return dictGenHashFunction((unsigned char *)key,sdslen((sds) key));
}

void *mykeyDup(void *privdata, const void *key){
	return (void *)sdsdup((sds) key);
}

int mykeyCompare(void *privdata, const void *key1, const void *key2){
	return sdscmp((sds) key1,(sds) key2);
}

void mykeyDestructor(void *privdata, void *key){
	sdsfree((sds) key);
}

void myvalDestructor(void *privdata, void *obj){
	free(obj);
}

long long usec(void) {
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (((long long)tv.tv_sec)*1000000)+tv.tv_usec;
}

void display_usage(char * argv[]){
	char *buffer=(char* )malloc(8092*sizeof(char));
	const char* usage=
"\nCopyright (c) 2012-2013\n" \
"Contact: XiongXu <xuxiong19880610@163.com> <xiongxu@me.com> \n" \
"Usage: %s [-o OUTFILE] [-r chr1:1-2000000] [-h] bamFile1 bamFile2 ..\n" \
"Discription:\n  This program is used for converting a bam file to wig file format, the bamfile should be samtools indexed. If you don't specify the -r parameter ,the program will convert the whole .bam file by chromosomes\n" \
"Example1:\n  %s -o out -r chr1:1-100000 /share/work1/staff/xuxiong/twins/batch2/test_20131224/RT144FD_L1_I001.R1.clean.fastq.gz_1224/accepted_hits.bam\n" \
"Example2:\n  %s -o out /share/work1/staff/xuxiong/twins/batch2/test_20131224/RT144FD_L1_I001.R1.clean.fastq.gz_1224/accepted_hits.bam\n" \
"\n" \
"   [-o OUTPUT_FILE]  = OUTPUT file.                                   [required]\n" \
"   [-r]              = region, default is whole genome.(chr1:1-20000) [option]\n" \
"   [-h]              = This helpful help screen.                      [option]\n" \
"   Infiles           = bam format input file(s),at least 1 bam file.  [required]\n" \
"\n";
	sprintf(buffer,usage,argv[0],argv[0],argv[0]);
	fprintf(stderr,"%s",buffer);
	exit(1);
}

int insertInt2Hash(dict *Start,uint32_t pos) {
	sds value=int2sds(pos);
	dictEntry *val_data=dictFind(Start,value);
	if (val_data==NULL) {
		int *init=(int *)malloc(sizeof(int));
		*init=1;
		dictAdd(Start,value,init);
	}
	else{
		(*(int *)val_data->val)++;
	}
	sdsfree(value);
	return 0;
}

int fetch_func(const bam1_t *b, void *data) {
	DataBuf *databuf=(DataBuf *)data;
	uint32_t *cigar = bam1_cigar(b);
	const bam1_core_t *c = &b->core;
	if (c->flag&BAM_DEF_MASK || b->core.tid < 0) return 0;
	int i, l;
	unsigned int temp_start= c->pos;
	for (i = l = 0; i < c->n_cigar; ++i) {
		int op = bam_cigar_op(cigar[i]);
		if (op == BAM_CINS) continue;
		if ( op == BAM_CDEL || op == BAM_CREF_SKIP){
			l = bam_cigar_oplen(cigar[i]);
			temp_start+=l;
		}
		else if (op == BAM_CMATCH ) {
			l = bam_cigar_oplen(cigar[i]);
			insertInt2Hash(databuf->Start,temp_start);
			temp_start+=l;
			insertInt2Hash(databuf->End,temp_start);
		}
	}
	databuf->chrRC[c->tid]++;
	databuf->chrReadLen[c->tid]+=l;
	return 0;
}

bam_index_t *load_bam_index(char *filename){
	bam_index_t *idx;
	if ((idx = bam_index_load(filename)) == 0) {
		fprintf(stderr, "bam2bed: BAM indexing file is not available.\n");
		exit(1);
	}
	return idx;
}

int bam_fetch_chr(bam_index_t *idx,const char *region,DataBuf *databuf){
	int ref, beg, end;
	bam_parse_region(databuf->fp->header,region, &ref, &beg, &end);
	if (ref < 0) {
		fprintf(stderr, "bam2bed: Invalid region %s\n", globalArgs.region);
		exit(1);
	}
//	fprintf(stderr,"%s\t%d\t%d\n",databuf->fp->header->target_name[ref],beg,end);
	bam_fetch(databuf->fp->x.bam, idx, ref, beg, end, databuf, fetch_func);
	return ref;
}

int cmpInt(const void *a, const void *b){
	return sds2int(*(sds const*)a)-sds2int(*(sds const*)b);
}

sds *union_hashed_keys(dict *hashtbl_a,dict *hashtbl_b){
	sds *keys=(sds *)calloc((hashtbl_a->used+hashtbl_b->used),sizeof(sds));
	int j=0;
	dictIterator *iter=dictGetIterator(hashtbl_a);
	dictEntry *entry=dictNext(iter);
	for (;entry;entry=dictNext(iter) ){
		keys[j++]=(sds)(entry->key);
	}
	dictReleaseIterator(iter);

	iter=dictGetIterator(hashtbl_b);
	entry=dictNext(iter);
	for (;entry;entry=dictNext(iter) ){
		keys[j++]=(sds)(entry->key);
	}
	dictReleaseIterator(iter);
	qsort(keys,j,sizeof(sds),cmpInt);
	return keys;
}

void hash2BedGraph(DataBuf *databuf,int ref,FILE *stream){
	unsigned long all_keys_count = databuf->Start->used + databuf->End->used;
	if (!all_keys_count) return;
	uint32_t i=0,prevkey=0,last_start=0,last_end=0,last_depth=0,Count=0,pos=0;
	sds *all_keys=union_hashed_keys(databuf->Start,databuf->End);
	for (i = 0; i < all_keys_count; ++i) {
		pos=sds2int(all_keys[i]);
		if (pos && pos==prevkey) continue;
		if (prevkey) {
			if (last_depth==Count) {
				prevkey=last_start;
			}
			else{
				if (last_depth) fprintf(stream,"%s\t%u\t%u\t%u\n",databuf->fp->header->target_name[ref],last_start,last_end,last_depth);
			}
		}
		last_start=prevkey;
		last_end=pos;
		last_depth=Count;
		dictEntry *data=dictFind(databuf->Start,all_keys[i]);
		if (data!=NULL) Count+= *((int*)data->val);
		data=dictFind(databuf->End,all_keys[i]);
		if (data!=NULL) Count-= *((int*)data->val);
		prevkey=pos;
	}
	if (last_depth){
		fprintf(stream,"%s\t%u\t%u\t%u\n",databuf->fp->header->target_name[ref],last_start,last_end,last_depth);
	}
	fprintf(stdout,"%s\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu64 "\n",databuf->fp->header->target_name[ref],databuf->fp->header->target_len[ref],databuf->chrRC[ref],databuf->chrReadLen[ref]);
	free(all_keys);
}

void region_sam2wig(DataBuf *databuf ,bam_index_t *idx,const char *region,FILE *bedGraph,dictType *mydictType){
	databuf->Start=dictCreate(mydictType,NULL);
	databuf->End=dictCreate(mydictType,NULL);
	int ref=bam_fetch_chr(idx,region,databuf);
	hash2BedGraph(databuf,ref,bedGraph);
	dictRelease(databuf->Start);
	dictRelease(databuf->End);
}

int main(int argc, char *argv[])
{
	int opt = 0;
	globalArgs.infiles=NULL;
	globalArgs.numInfiles=0;
	globalArgs.outfile="-";
	globalArgs.region="-";
	const char *optString = "o:r:h?";
	if (argc<2) display_usage(argv);
	opt = getopt( argc, argv, optString );
	while( opt != -1 ) {
		switch( opt ) {
			case 'o':
				globalArgs.outfile = optarg;
				break;
			case 'r':
				globalArgs.region = optarg;
				break;
			case '?':	/* fall-through is intentional */
				break;
			case 'h':
				display_usage(argv);
				break;
			default:
				fprintf(stderr,"error parameter!\n");
				break;
		}
		opt = getopt( argc, argv, optString );
	}
	globalArgs.infiles = argv + optind;
	globalArgs.numInfiles = argc - optind;
	long long begin;
	begin=usec();

	DataBuf *databuf=(DataBuf *)calloc(globalArgs.numInfiles,sizeof(DataBuf));
	FILE *bedGraph=fcreat_outfile(globalArgs.outfile,".bedGraph");
	dictType *mydictType=(dictType *)malloc(sizeof(dictType));
	initdictType(mydictType);
	int32_t i=0,j=0;
	for (i=0;i<globalArgs.numInfiles ;i++ ) {
		if (((databuf+i)->fp = samopen(*(globalArgs.infiles+i), "rb", 0)) == 0) {
			err(1, "bam2bed: Fail to open BAM file %s\n", *(globalArgs.infiles+i));
		}
		bam_index_t *idx=load_bam_index(*(globalArgs.infiles+i));
		(databuf+i)->chrRC = (uint32_t *)calloc((databuf+i)->fp->header->n_targets,sizeof(uint32_t));//databuf->fp->header->target_len
		(databuf+i)->chrReadLen =(uint64_t *)calloc((databuf+i)->fp->header->n_targets,sizeof(uint64_t));
		if (strncmp(globalArgs.region,"-",6)) { /* if a region is not specified */
			region_sam2wig(databuf+i ,idx,globalArgs.region,bedGraph,mydictType);
		}else{
			for (j=0;j<(databuf+i)->fp->header->n_targets ;j++ ) {
				region_sam2wig(databuf+i ,idx,(databuf+i)->fp->header->target_name[j],bedGraph,mydictType);
				fprintf(stderr,"%s at %.3f s\n",(databuf+i)->fp->header->target_name[j],(double)(usec()-begin)/CLOCKS_PER_SEC);
			}
		}
		free((databuf+i)->chrRC);
		free((databuf+i)->chrReadLen);
		bam_index_destroy(idx);
		samclose((databuf+i)->fp);
		fprintf(stderr,"Converted %s to wig format at %.3f s\n",*(globalArgs.infiles+i),(double)(usec()-begin)/CLOCKS_PER_SEC);
	}
	free(mydictType);
	free(databuf);
	fclose(bedGraph);
	return 0;
}
