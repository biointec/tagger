/* *********************************************************************
   Check that file argv[2] contains the BWT of file argv[1]
   Optionally check also a suffix array and a lcp array files.
   
   The checking is done computing the SA of argv[1] and, to test a lcp file, 
   also the LCP array. The SA is computed with O(n\log n) divsufsort so the 
   only requirement is that the text consist of 1 byte symbol (no distinct EOS
   required as for the linear time SACA and SAIS algorithms). 
   
   The LCP array is computed from the SA using the lcp9 algorithm from 
   SWAT '04: not the fastest available but a good balance between simplicity,
   space usage and speed. 
   Works for files of size up to 2^40 since the SA and LCP array are
   stored using 40 bits per entry, stored together in 64+16 bits
   for a total space usage of 11n bytes.  
      
   Compiling with -DVERY_BIG uses 48 bits per each SA and LCP entry
   stored together in 64+32 bits: the total space usage becomes 
   13n bytes but the input files can have size op to 2^48 (not tested). 
   
   The SA and LCP entries in the files must be stored using a fixed
   number of bytes per entry between 4 and 8 (controlled by option -i)
   The program was tested only for 5 bytes x entry   
   
   A note on the running time for SA construction. Although divsufsort 
   runs Th(n \log n) in the worst case it appears to be faster in 
   practive than the linea time SACA-K algorithm. For a 3GB collection
   for human chr19 variants SA construction with divsufsort took 460 
   seconds vs 1025 seconds for SACA-K. For einstein.en.txt divsufsort 
   took 49 seconds vs 140 seconds for SACA-K 
  *********************************************************************** */
#define _GNU_SOURCE  
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <time.h>
#include "sa64/divsufsort64.h"


#define ALPHABET_SIZE 256


#ifndef VERY_BIG
// 5 byte per entry version
typedef uint16_t ulcpL_t;   
#define LOW_BITS 16
#define ALL_BITS 40
#define HIGH_BITS 24
#define MASKSA 0xFFFFFFFFFFL     // used to read a SA value 
#define MASKHIGH 0xFFFFFF0000L   // used to read the lcp high bits
#else
#error "The 6 byte version is not currently used" 
// 6 byte per entry version 
typedef uint32_t ulcpL_t; 
#define LOW_BITS 32
#define ALL_BITS 48
#define HIGH_BITS 16
#define MASKSA 0xFFFFFFFFFFFFL
#define MASKHIGH 0xFFFF00000000L
#endif

// macro to read a single LCP value combing SA and LCP values  
#define ReadLCP(i) (((SA[(i)]>>HIGH_BITS)&MASKHIGH)|LCP[(i)])


// vebosity level
int Verbose;

// local prototypes
static void print_help(const char *progname, int status);
static long get_file_size(FILE *f, const char *fname);
static void check_int(FILE *f,long i, int bytes, long pos, char *name);
static void check_char(FILE *f,int c, long pos);
static void sa2lcp_9n(sauchar_t *t, saidx64_t n, saidx64_t *sa, ulcpL_t *lcp);
static saidx64_t sa2ranknext(sauchar_t *t, saidx64_t n, saidx64_t *sa, ulcpL_t *lcp);



int main(int argc, char *argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;
  FILE *fp, *fbwt;
  const char *fname, *bwtname;
  sauchar_t *T;
  saidx64_t *SA;
  ulcpL_t *LCP=NULL;
  saidx64_t n;
  clock_t start0, start, finish;

  
  start0 = clock();
  Verbose=0;  // default values
  char *saname = NULL, *lcpname=NULL;
  int c, ibytes = 5;
  while ((c=getopt(argc, argv, "vS:L:i:")) != -1) {
    switch (c) 
      {
      case 'v':
        Verbose++; break;
      case 'S':
        saname=optarg; break;
      case 'L':
        lcpname=optarg; break;        
      case 'i':
        ibytes=atoi(optarg); break;        
      case '?':
        print_help(argv[0],EXIT_FAILURE);;
      }
  }
  if(optind+2==argc) {
    fname=argv[optind];
    bwtname=argv[optind+1];
  }
  else  print_help(argv[0],EXIT_FAILURE);
  // --- check ibytes
  if(ibytes<4 || ibytes>8) {
    fprintf(stderr, "%s: Invalid ibytes parameter (%d); supported range i 4-8\n", argv[0], ibytes);
    exit(EXIT_FAILURE);
  }
    

  // ---------- make sure bwt file exists
  if((fbwt = fopen(bwtname, "rb")) == NULL) {
    fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], bwtname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  // ----- make sure SA and LCP file exist if requested 
  FILE *fsa = NULL;
  if(saname!=NULL && (fsa = fopen(saname, "rb")) == NULL) {
    fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], saname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  FILE *flcp = NULL;
  if(lcpname!=NULL && (flcp = fopen(lcpname, "rb")) == NULL) {
    fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], lcpname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  
  //--------- store text in memory 
  /* Open a file for reading. */
  if((fp = fopen(fname, "rb")) == NULL) {
    fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], fname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  /* Get text file size. */
  n = (saidx64_t) get_file_size(fp,fname);
  assert(n>0);
  if(n>MASKSA) {
    fprintf(stderr, "%s: input file %s is too big: max size is 2^%d-1\n", argv[0],fname,ALL_BITS);
    exit(EXIT_FAILURE);
  }

  // --- check that bwt/SA/LCP files are of the right size
  long nbwt = get_file_size(fbwt,bwtname);
  if(n+1!=nbwt) {
    fprintf(stderr, "%s: bwt file has the wrong size: %ld(file) vs %ld(expected)\n", argv[0],nbwt,n+1);
    exit(EXIT_FAILURE);
  }
  if(saname!=NULL) {
    long nsa = get_file_size(fsa,saname);
    if(n*ibytes!=nsa) {
      fprintf(stderr, "%s: sa file has the wrong size: %ld(file) vs %ld(expected)\n", argv[0],nsa,ibytes*n);
      exit(EXIT_FAILURE);
    }
  } 
  if(lcpname!=NULL) {
    long nlcp = get_file_size(flcp,lcpname);
    if(n*ibytes!=nlcp) {
      fprintf(stderr, "%s: lcp file has the wrong size: %ld vs %ld\n", argv[0],nlcp,ibytes*n);
      exit(EXIT_FAILURE);
    }
  } 


  /* Allocate Space for Text and SA */
  SA = (saidx64_t *)malloc((size_t)n * sizeof(saidx64_t));
  T = (sauchar_t *)malloc((size_t)n * sizeof(sauchar_t));
  // if using SA-1 later on is a problem use something like:
  // saidx64_t *sa1 = (saidx64_t *)malloc((size_t)n * sizeof(saidx64_t)+1);
  // SA = sa1+1;
  if((T == NULL) || (SA == NULL)) {
    fprintf(stderr, "%s: Cannot allocate memory for Text and SA.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  /* Read n bytes of data. */
  if(fread(T, sizeof(sauchar_t), (size_t)n, fp) != (size_t)n) {
    fprintf(stderr, "%s: %s `%s': ",
      argv[0],
      (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
      fname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
  fclose(fp);

  /* Construct the suffix array. */
  start = clock();
  if(divsufsort64(T, SA, (saidx64_t)n) != 0) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  finish = clock();
  if(Verbose)
    fprintf(stderr, "SA computation: %.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);

  /* Check the BWT and suffix array. */
  int bwtc = T[n-1];         // first bwt char
  check_char(fbwt,bwtc,-1);  // check first bwt char 
  for(long i=0; i<n; i++) {
    saidx64_t s = SA[i];
    if(fsa) check_int(fsa,s,ibytes,i,"Suffix array");
    if(s>0) {
      bwtc = T[s-1];
      check_char(fbwt,bwtc,i);
    }
    else { // eos positon 
      int eos = fgetc(fbwt);
      if(Verbose>0) fprintf(stderr,"Eos character in file %s: %d\n",bwtname, eos);
    }
  }
  // check we are really at the end of teh file  
  assert(fgetc(fbwt)==EOF && feof(fbwt));  
  if(Verbose)
    fprintf(stderr, "Bwt file %s is correct\n",bwtname);
  if(fsa) {
    assert(fgetc(fsa)==EOF && feof(fsa));  
    if(Verbose)
      fprintf(stderr, "Suffix array file %s is correct\n",saname);
  }


  // compute the lcp
  if(flcp!=NULL) {
    LCP = (ulcpL_t *)malloc((size_t)n * sizeof(ulcpL_t));
    if((LCP==NULL)) {
      fprintf(stderr, "%s: Cannot allocate memory for LCP array.\n", argv[0]);
      exit(EXIT_FAILURE);
    }
    start = clock();
    sa2lcp_9n(T,n,SA-1,LCP-1);
    finish = clock();
    if(Verbose)
      fprintf(stderr, "LCP computation: %.4f sec\n", (double)(finish - start) / (double)CLOCKS_PER_SEC);
    // --- first LCP value is not real a LCP; it is the maxlcp in sa2lcp_9n
    saidx64_t maxlcp = ReadLCP(0);
    if(Verbose)
      fprintf(stderr, "Largest LCP value: %ld\n", maxlcp);
    fseek(flcp,ibytes,SEEK_SET); // skip first LCP file entry     
    for(saidx64_t i=1;i<n;i++) {
      saidx64_t lcpi = ReadLCP(i);
      assert(lcpi<=maxlcp); 
      check_int(flcp,lcpi,ibytes,i,"LCP array");
    }
    // check we are really at the end of the file  
    assert(fgetc(flcp)==EOF && feof(flcp));  
    if(Verbose)
      fprintf(stderr, "LCP file %s is correct\n",lcpname);
  }

  finish = clock();
  fprintf(stderr, "Total elapsed time: %.4f sec\n", (double)(finish - start0) / (double)CLOCKS_PER_SEC);  

  /* Deallocate memory. */
  free(T);
  free(SA);
  if(LCP) free(LCP);

  return 0;
}



/* ***********************************************************************
   WARNING sa, lcp, rank_next are 1 based: sa[n] is ok sa[0] is segfault  
   space economical computation of the lcp array
   This procedure is similar to the algorithm of kasai et al.  but we compute 
   the rank of i using the rank_next map (stored in the rn array) instead of
   the rank array. The advantage is that as soon as we have read rn[k] that
   position is no longer needed and we can use it to store the lcp.
   Thus, rn[] and lcp[] share the same memory.
   It is assumed that n < 2^40 (2^48 if VERY_BIG is defined) so each 
   sa/lcp/rn entry is stored in 5 (6) bytes 
   The overall space requirements of the procedure is
      n + 5n (sa) + 5n (lcp/rn) = 11n
      (or n + 6n (sa) + 6n (lcp/rn) = 13n if VERY_BIG)
   The name 9n is historical and denotes the algorithm as defined in the SWAT 04 paper 
    input
      t[0,n-1] input text 
      sa[1,n]  suffix array each entry taking at most 40 bits, the remaining 24 bits are 0
      lcp[1,n] uninitialized array of uint16 (or uint32 if VERY_BIG)
    return
      nothing
      the lcp[2,n] entries are stored in the high 24 bit of each sa[] entry
                   and in the 16 bit of each lcp[] entry 
                   sa[j][40-63] contains bits 16-39, while lcp[j] bits 0-15
      lcp[i] is lcp between t[sa[i-1]...] and t[sa[i]...]
    no additional space in addition to t[], sa[], lcp[] is used 
   *********************************************************************** */
static void sa2lcp_9n(sauchar_t *t, saidx64_t n, saidx64_t *sa, ulcpL_t *lcp)
{
  saidx64_t i,h,j,k,nextk=-1;
  ulcpL_t *rn;
  saidx64_t maxlcp = 0;
  
  // rn and lcp are pointers to the same array
  rn =  lcp;
  // compute rank_next map
  k = sa2ranknext(t, n, sa, rn); // k is the rank of t[0...]
  if(Verbose>1)
    fprintf(stderr,"ranknext computation completed\n");
  for(h=i=0; i<n; i++,k=nextk) {
    assert(k>0 && i==(sa[k]&MASKSA)); 
    nextk=((sa[k]>>HIGH_BITS)&MASKHIGH)|rn[k];     // read nextk before it is overwritten
    if(k>1) {
      // compute lcp between suffixes of rank k and k-1 (recall i==sa[k])
      j = sa[k-1] & MASKSA;
      while(i+h<n && j+h<n && t[i+h]==t[j+h])
        h++;
      // check lcp fits in 40/48 bits   
      if(h>MASKSA) {
        fprintf(stderr,">%d bits lcp value: %ld. Exiting\n",ALL_BITS, (long) h); 
        exit(EXIT_FAILURE);
      }
      if(h>maxlcp) maxlcp=h;           // compute also maximum LCP 
      lcp[k]=h;        // save bits 0:15
      sa[k] = (h>>LOW_BITS)<<ALL_BITS | (sa[k]&MASKSA); // save bits 16:39   
    }
    if(h>0) h--;
  }
  assert(nextk==0);        // we have reached the last suffix s[n-1]
  lcp[1] = maxlcp;         // write maxlcp to lcp[1]
  sa[1] = ((maxlcp>>LOW_BITS)<<ALL_BITS) | (sa[1]&MASKSA);
  return;
}


/* *******************************************************************
   WARNING sa, lcp, rank_next are 1 based 
   input 
     n size of text
     t[0,n-1] input text
     sa[1,n] suffix array (low 32 bits)
  output
    rank_next[1,n] is written to rn[1,n] (32 bits) and in the high 16 bits of sa[1,n] 
    note: ranks are in the range [1,n]
  return 
    r0 = rank of t[0,n-1] (it is not a rank_next value), 
    sa[r0]=0 (this is an invalid  rank)
  space: inplace except for the occ and count arrays (occ can be avoided)
  ******************************************************************* */
static saidx64_t sa2ranknext(sauchar_t *t, saidx64_t n, saidx64_t *sa, ulcpL_t *rn)
{
  saidx64_t i,j,sai,c,r0=0;
  saidx64_t count[ALPHABET_SIZE], occ[ALPHABET_SIZE];

  // compute # occurrences of each char 
  for(i=0;i<ALPHABET_SIZE;i++) occ[i]=0;
  for(i=0;i<n;i++) occ[t[i]]++;
  // --- occ -> count
  count[0]=0;
  for(i=1;i<ALPHABET_SIZE;i++) {
    count[i] = count[i-1] + occ[i-1];
  }
  
  // --- sa+t -> rank_next
  j = ++count[t[n-1]];       // this is bwt[0]
  assert(j>0 && j<=n);
  sa[j] |= 0; // here there was rank_next[j]=0, note 0 is an invalid rank value;
  for(i=1;i<=n;i++) {
    sai = (sa[i] & MASKSA);
    assert(sai>=0 && sai<n);
    if(sai == 0)
      r0 = i;
    else {
      c = t[sai-1];
      j = ++count[c];  // rank_next[j] = i
      assert(j>0 && j<=n); 
      rn[j] = i;                // save low bits 0..15
      if(i>>LOW_BITS)           // is there any high bit set?
        sa[j] |= ((i>>LOW_BITS)<<ALL_BITS); // save bits 16..39 of i in sa[j][40..63]
    }
  }
  assert(r0>0);
  return r0;
}

static void check_int(FILE *f,long i, int bytes, long pos, char *name) {
  long fi=0;
  int e = fread(&fi,bytes,1,f);
  if(e!=1) { 
    fprintf(stderr,"Cannot read %d-byte %s integer\n",bytes,name); 
    perror(NULL);
    exit(EXIT_FAILURE);
  }    
  if(fi!=i) {
    fprintf(stderr,"%s mismatch at sa position %ld, %ld(file) vs %ld(real)\n",name, pos,fi,i);
    exit(EXIT_FAILURE);
  }
}

static void check_char(FILE *f,int c, long pos) {
  int fc = fgetc(f);
  if(fc!=c) {
    fprintf(stderr,"bwt char mismatch at sa position %ld, %d(file) vs %d(real)\n",pos,fc,c);
    exit(EXIT_FAILURE);
  }
}


// get file size leaving the file in a rewind status  
static 
long get_file_size(FILE *f, const char *fname) {
  if(fseek(f, 0, SEEK_END) == 0) {
    long nx = ftell(f);
    rewind(f);
    if(nx < 0) {
      fprintf(stderr, "Cannot ftell `%s': ",fname);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
    return nx;
  }
  else {
    fprintf(stderr, "Cannot fseek `%s': ", fname);
    perror(NULL);
    exit(EXIT_FAILURE);
  }
}



static void
print_help(const char *progname, int status) {
  fprintf(stderr, "Usage:\n\t %s textfile bwtfile [options]\n\n", progname);
  fprintf(stderr, "Compare BWT(textfile) with the content of bwtfile. Optionally also check\n");
  fprintf(stderr, "the Suffix/LCP Array assuming they are written using 'size' bytes x entry.\n");
  fprintf(stderr, "The program terminates as soon as a mismatch is encountered\n\n"); 
  fprintf(stderr, "Options:\n"); 
  fprintf(stderr, "\t-S safile     compare safile with suffix array of textfile\n");
  fprintf(stderr, "\t-L lcpfile    compare lcpfile with the lcp array of textfile\n");
  fprintf(stderr, "\t-i size       bytes x entry in sa/lcp files (def. 5)\n");
  fprintf(stderr, "\t-v verbose output\n\n");
  exit(status);
}


