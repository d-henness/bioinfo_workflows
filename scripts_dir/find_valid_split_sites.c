#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char *argv[]){
  int opt;
  int start = 0, finish = 0;
  char *filenm;

  // i = input file, x = smallest position, y = largest position
  while((opt = getopt(argc, argv, "i:x:y:")) != -1){
    switch(opt){
      case 'i':
        filenm = optarg;
        break;
      case 'x':
        start = strtol(optarg, NULL, 10);
        break;
      case 'y':
        finish = strtol(optarg, NULL, 10);
        break;
      case '?':
        return EXIT_FAILURE;
        break;
      default:
        printf("could not parse %s\n", optarg);
        break;
    }
  }
  if ((start == 0) || (finish == 0) || (finish < start) || (filenm == NULL)){
    fprintf(stderr, "Invalid arguments\n");
    return EXIT_FAILURE;
  }
  printf("file %s start %d finish %d\n", filenm, start, finish);

  FILE *gtf_file;
  if ((gtf_file = fopen(filenm, "r")) == NULL){
    fprintf(stderr, "Could not open file %s\n", filenm);
    return EXIT_FAILURE;
  }

  int *invalid_sites;
  invalid_sites = calloc(finish - start + 1, sizeof(int));
  int small = 0, large = 0;
  while(fscanf(gtf_file, "%*s %*s %*s %d %d %*[^\n]", &small, &large) != -1){
    for (int i = small - start; i < large - start; i++){
      *(invalid_sites + i) = 1;
    }
  }

  for(int i = 0; i < finish - start; i++){
    if (!(*(invalid_sites + i))){
      printf("%d\n", i + start);
    }
  }
  free(invalid_sites);

  return 0;
}
