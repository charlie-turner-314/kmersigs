#include "uthash.h"
#include <chrono>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ISAAC-rand.h"

typedef unsigned char byte;

#define SIGNATURE_LEN 64

int DENSITY = 21;
int PARTITION_SIZE;

int WORDLEN;
FILE *sig_file;

typedef struct {
  char term[100];
  short sig[SIGNATURE_LEN];
  UT_hash_handle hh;
} hash_term;

typedef struct Signature {
  // int for document number
  int doc;
  // signature len bits for sig
  byte sig[SIGNATURE_LEN / 8];
} Signature;

hash_term *vocab = NULL;
#pragma omp threadprivate(vocab)

#define min(a, b) ((a) < (b) ? (a) : (b))

int non_zero = SIGNATURE_LEN * DENSITY / 100;

struct randctx ctx;
// make ctx thread private
#pragma omp threadprivate(ctx)

short *compute_new_term_sig(char *term, short *term_sig) {
  seed_random(term, WORDLEN, &ctx);

  int positive = 0;
  while (positive < non_zero / 2) {
    short pos = random_num(SIGNATURE_LEN, &ctx);
    if (term_sig[pos] == 0) {
      term_sig[pos] = 1;
      positive++;
    }
  }

  int negative = 0;
  while (negative < non_zero / 2) {
    short pos = random_num(SIGNATURE_LEN, &ctx);
    if (term_sig[pos] == 0) {
      term_sig[pos] = -1;
      negative++;
    }
  }
  return term_sig;
}

short *find_sig(char *term) {
//  hash_term *entry;
//  HASH_FIND(hh, vocab, term, WORDLEN, entry);
//  if (entry == NULL) {
//    entry = (hash_term *)malloc(sizeof(hash_term));
//    strncpy(entry->term, term, WORDLEN);
//    memset(entry->sig, 0, sizeof(entry->sig));
//    compute_new_term_sig(term, entry->sig);
//    HASH_ADD(hh, vocab, term, WORDLEN, entry);
//  }
  short *sig = (short *)malloc(sizeof(short) * SIGNATURE_LEN);
  compute_new_term_sig(term, sig);

  return sig;
}

void signature_add(char *term, int *doc_sig) {
  short *term_sig = find_sig(term);
#pragma omp simd
  for (int i = 0; i < SIGNATURE_LEN; i++)
    doc_sig[i] += term_sig[i];
  free(term_sig);
}

void compute_signature(char *sequence, int length, Signature *sig) {
  int doc_sig[SIGNATURE_LEN] = {0};

  for (int i = 0; i < length - WORDLEN + 1; i++)
    signature_add(sequence + i, doc_sig);

  // flatten and output to sig
  for (int i = 0; i < SIGNATURE_LEN; i += 8) {
    byte c = 0;
    for (int j = 0; j < 8; j++)
      c |= (doc_sig[i + j] > 0) << (7 - j);
    sig->sig[i / 8] = c;
  }
}

typedef struct Stat {
  long parts;
  long lines;
  long size;
} Stat;

Stat num_parts(FILE *file) {
  // basically do what partition does, except just count the number of them
  // and return that
  char buffer[10000];
  Stat stats = {0, 0};
  while (fgets(buffer, 10000, file) != NULL) {
    fgets(buffer, 10000, file);
    stats.lines++;
    int length = strlen(buffer);
    int i = 0;
    do {
      stats.parts++;
      i += PARTITION_SIZE / 2;
    } while (i + PARTITION_SIZE / 2 < length);
  }
  fseek(file, 0, SEEK_END);
  stats.size = ftell(file);
  rewind(file);
  return stats;
}

void process(char *file, long num_lines, Signature *signatures) {
  // go through each partition, and create an omp task for computing the
  // signature line number is document number
#pragma omp parallel
  {
#pragma omp single
    {
      int part_num = 0;
      for (int doc = 0; doc < num_lines; doc++) {
        // skip one line:
        char *line = strsep(&file, "\n");
        line = strsep(&file, "\n");
        int length = strlen(line);
        // we have a line, now need to partition and compute signature
        int j = 0;
        do {
// create a task for each partition
#pragma omp task
          {
            signatures[part_num].doc = doc;
            compute_signature(line + j, min(PARTITION_SIZE, length - j),
                              &signatures[part_num]);
          }
          j += PARTITION_SIZE / 2;
          part_num++;
        } while (j + PARTITION_SIZE / 2 < length);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <filename> <wordlen> <partition_size>\n",
            argv[0]);
    return 1;
  }

  const char *filename = argv[1];
  WORDLEN = atoi(argv[2]);
  PARTITION_SIZE = atoi(argv[3]);

  auto start = std::chrono::high_resolution_clock::now();

  FILE *file;
  file = fopen(filename, "r");
  if (file == NULL) {
    fprintf(stderr, "Error: failed to open file %s\n", filename);
    return 1;
  }

  char outfile[256];
  snprintf(outfile, 256, "%s.part%d_sigs%02d_%d_omp", filename, PARTITION_SIZE,
           WORDLEN, SIGNATURE_LEN);
  sig_file = fopen(outfile, "w");

  Stat stats = num_parts(file);

  // read in the file
  char *buffer = (char *)malloc(sizeof(char) * stats.size);
  size_t result = fread(buffer, 1, stats.size, file);
  if (result != stats.size) {
    fprintf(stderr, "Error: failed to read file %s\n", filename);
    return 1;
  }
  // need to have enough memory to store all the signatures
  Signature *signatures = (Signature *)malloc(sizeof(Signature) * stats.parts);

  // now we can process the file
  process(buffer, stats.lines, signatures);

  // write out the signatures to the file
  // want to write out the document number and the signature
  for (int i = 0; i < stats.parts; i++) {
    fwrite(&signatures[i].doc, sizeof(int), 1, sig_file);
    fwrite(signatures[i].sig, sizeof(byte), SIGNATURE_LEN / 8, sig_file);
  }

  fclose(file);

  fclose(sig_file);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  printf("%s %f seconds\n", filename, duration.count());

  return 0;
}
