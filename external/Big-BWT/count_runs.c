#define _GNU_SOURCE 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include <errno.h>   


// quick and dirty tool to count the number of runs in a binary file

// Input:  file_name and number of bytes x symbol
// Output: # runs, maxrun, etc on stdout



void termina(const char *messaggio)
{
  if(errno==0) fprintf(stderr,"%s\n",messaggio);
  else         perror(messaggio);
  exit(1);
}
 

int main(int argc, char *argv[])
{
  if(argc!=3) {
    printf("Usage: %s filename bytes_x_symbol\n",argv[0]);
    return 1;
  }
  int bytes = atoi(argv[2]);
  if(bytes<1) 
    termina("Number of bytes x symbol must be positive");
  
  // open file
  FILE *f = fopen(argv[1],"r");
  if(f==NULL) termina("Unable to open file");

  long sym=0;
  int e = fread(&sym,1,bytes,f);
  if(e!=bytes) termina("Unable to read the first symbol");
  long runs = 1;
  long totsym = 1;
  long curlen = 1;
  long maxlen = 0;  
  while(true) {
    long newsym = 0;
    e = fread(&newsym,1,bytes,f);
    if(e==0) break;
    if(e<bytes) {
      fprintf(stderr,"Unable to read a full symbol. Bytes x symbol is wrong?\n");
      break;
    }
    // new symbol is valid
    totsym++;
    if(newsym!=sym) {
      if(curlen>maxlen) maxlen = curlen;
      curlen = 1;
      runs +=1;
      sym = newsym;
    }
    else curlen++;
  }
  printf("File: %s   [%d byte(s) x symbol]\n",argv[1],bytes);
  printf("Number of symbols: %ld\n",totsym);
  printf("Number of runs: %ld\n",runs);
  printf("Longest run: %ld (ave length: %lf)\n",maxlen,1.0*totsym/runs);
  return 0;
}
