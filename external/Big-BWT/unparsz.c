#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <fcntl.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <time.h>
#include <sys/mman.h>
#include "utils.h" 

// -------------------------------------------------------------
// struct containing command line parameters and other globals
typedef struct {
   char *basename;
   char *outname;
} Args;


static void print_help(char *name)
{
  printf("Usage: %s <basename> [options]\n\n", name);
  puts("Restore the original file given a quasi prefix-free parse");
  puts("ie files ."EXTDICZ" ."EXTDZLEN" and ."EXTPARSE);
  puts("  Options:");
  puts("\t-o outfile   output file (def. <basename>.out)");
  puts("\t-w wsize     window size (def. 10)");
  puts("\t-h           show help and exit");
  exit(1);
}

static void parseArgs(int argc, char** argv, Args *arg ) {
  extern int optind, opterr, optopt;
  extern char *optarg;  
  int c;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("\n");

  arg->outname = NULL;
  while ((c = getopt( argc, argv, "ho:") ) != -1) {
    switch(c) {
      case 'o':
      arg->outname = strdup(optarg); break;
      case 'h':
         print_help(argv[0]); exit(1);
      case '?':
      puts("Unknown option. Use -h for help.");
      exit(1);
    }
  }
  // read base name as the only non-option parameter   
  if (argc!=optind+1) 
    print_help(argv[0]);
  arg->basename = strdup(argv[optind]);
  // if not given create output file name with .out extension
  if(arg->outname==NULL) {
    int e = asprintf(&arg->outname,"%s.out",arg->basename);
    if(e<0) die("Error creating output file name");
  }
}

void *mmap_fd(int fd, size_t *n)
{
  off_t off = lseek(fd,0,SEEK_END);
  if(n<0) die("seek error on dictionary file");
  void *a = mmap(NULL,off,PROT_READ|PROT_WRITE,MAP_PRIVATE,fd,0);
  if(a==MAP_FAILED) die("mmap error on dictionary file");
  *n =off;
  return a;
}

int main(int argc, char *argv[])
{
  Args arg;        // command line arguments
  char *Dict;      // dictionary 
  size_t n;        // length of dictionary
  uint32_t *Wstart; // starting positions of words inside dictionary 
  
  // read arguments 
  parseArgs(argc,argv,&arg);
  // start measuring wall clock time 
  time_t start_wc = time(NULL);
  
  // mmap dictionary file to memory
  int dict_fd = fd_open_aux_file(arg.basename,EXTDICZ,O_RDONLY);
  Dict = mmap_fd(dict_fd,&n);
  if(close(dict_fd)!=0) die("error closing dictionary file");

  // compute # words, and save starting points
  int wlen_fd = fd_open_aux_file(arg.basename,EXTDZLEN,O_RDONLY);
  long size = lseek(wlen_fd, 0, SEEK_END); // move to end of gile and get position
  assert(size%4==0); 
  long words = size/4;
  Wstart = malloc((words+1)*sizeof(*Wstart));
  if(Wstart==NULL) die("Allocation error");
  Wstart[0] = 0; // first word starts at Dict[0]
  lseek(wlen_fd,0,SEEK_SET);     // rewind 
  ssize_t r = read(wlen_fd,Wstart+1,4*words); // assume word lengths are 4 bytes
  if(r!=4*words) die("error reading from wlen file"); 
  for(long i=2;i<=words;i++)    // Wstart[i] is starting position of word i
    Wstart[i] += Wstart[i-1];   // Wstart[i+1]-1 is ending position of word i
  
  
  // === check starting from here 
  
  fprintf(stderr,"Found %ld dictionary words\n",words);
  fprintf(stderr,"Recovering file %s\n",arg.outname);
  // create output file reading word id's from the parse
  FILE *f = fopen(arg.outname,"wb");
  if(f==NULL) die("Cannot open output file");
  FILE *parse = open_aux_file(arg.basename,EXTPARSE,"rb");
  if(parse==NULL) die("Cannot open parse file");
  while(true) {
    uint32_t w; char *s;
    int e = fread(&w,4,1,parse);
    if(e==0 && feof(parse)) break; // done
    if(e!=1) die("Error reading parse file");
    if(w==0 || w>words) die("Invalid word ID in the parse file");
    s = Dict + Wstart[w-1]; // dictionary word starting position 
    e = fwrite(s,1,Wstart[w]-Wstart[w-1],f);
    if(e!=Wstart[w]-Wstart[w-1]) die("Error writing to the output file"); 
  }
  fclose(parse);
  fclose(f);
  free(Wstart);
  munmap(Dict,n);
  free(arg.outname);
  free(arg.basename);  
  printf("==== Elapsed time: %.0lf wall clock seconds\n", difftime(time(NULL),start_wc));  
  return 0;
}

