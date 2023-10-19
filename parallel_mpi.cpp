#include "uthash.h"
#include <chrono>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ISAAC-rand.h"

typedef unsigned char byte;

#define SIGNATURE_LEN (64 /* in bits */)

int rank, size;

int DENSITY = 21;   // Density of the signature (non-zero) (percent)
int WORDLEN;        // k-mer length (characters)
int PARTITION_SIZE; // size to partition sequences into (characters)
char out_filename[256];
// int doc_sig[SIGNATURE_LEN];

typedef struct {
  char term[100];
  short sig[SIGNATURE_LEN];
  UT_hash_handle hh;
} hash_term;

hash_term *vocab = NULL;

typedef struct Signature {
  // int for document number
  int doc;
  // signature len bits for sig
  byte sig[SIGNATURE_LEN / 8];
} Signature;

MPI_Comm WORKER_COM;

typedef struct Stat {
  long parts;
  long lines;
  long size;
} Stat;

typedef struct WorkInfo {
  long start;
  long end;
  int doc_start;
  int doc_end;
  int num_parts;
} WorkInfo;

Stat file_stats(FILE *file) {
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

int distribute(char *filename, WorkInfo *work_infos) {
  // rank 0 opens file, runs through it, and fills out the above
  // open file
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    MPI_Finalize();
    return 1;
  }
  Stat stats = file_stats(fp);
  // fill out work_infos, with approximately equal parts in each
  char buffer[10000];
  int doc = 0;
  int num_parts = 0;
  int parts_per_proc = stats.parts / size;
  int current_proc = 0;
  work_infos[current_proc].start = 0;
  work_infos[current_proc].doc_start = doc;
  int length;
  int i;
  while (fgets(buffer, 10000, fp) != NULL) {
    fgets(buffer, 10000, fp);
    length = strlen(buffer);
    i = 0;
    do {
      num_parts++;
      i += PARTITION_SIZE / 2;
    } while (i + PARTITION_SIZE / 2 < length);
    int bytes_read = ftell(fp);
    if (num_parts >= parts_per_proc) {
      // this is the last line for this process
      work_infos[current_proc].end = bytes_read;
      work_infos[current_proc].doc_end = doc;
      work_infos[current_proc].num_parts = num_parts;
      current_proc++;
      work_infos[current_proc].start = ftell(fp) + 1;
      work_infos[current_proc].doc_start = doc + 1;
      num_parts = 0;
    }
    doc++;
  }
  // if more processes than lines, set the rest to -1
  for (int i = current_proc + 1; i < size; i++) {
    work_infos[i].start = -1;
  }

  return 0;
}

short *compute_new_term_sig(char *term, short *term_sig) {
  seed_random(term, WORDLEN);

  int non_zero = SIGNATURE_LEN * DENSITY / 100;

  int positive = 0;
  while (positive < non_zero / 2) {
    short pos = random_num(SIGNATURE_LEN);
    if (term_sig[pos] == 0) {
      term_sig[pos] = 1;
      positive++;
    }
  }

  int negative = 0;
  while (negative < non_zero / 2) {
    short pos = random_num(SIGNATURE_LEN);
    if (term_sig[pos] == 0) {
      term_sig[pos] = -1;
      negative++;
    }
  }
  return term_sig;
}

short *find_sig(char *term) {
  // hash_term *entry;
  // HASH_FIND(hh, vocab, term, WORDLEN, entry);
  // if (entry == NULL) {
  //   entry = (hash_term *)malloc(sizeof(hash_term));
  //   strncpy(entry->term, term, WORDLEN);
  //   memset(entry->sig, 0, sizeof(entry->sig));
  //   compute_new_term_sig(term, entry->sig);
  //   HASH_ADD(hh, vocab, term, WORDLEN, entry);
  // }
  // return entry->sig;

  short *sig = (short *)malloc(sizeof(short) * SIGNATURE_LEN);

  compute_new_term_sig(term, sig);

  return sig;
}

void signature_add(char *term, int *doc_sig) {
  short *term_sig = find_sig(term);
  for (int i = 0; i < SIGNATURE_LEN; i++)
    doc_sig[i] += term_sig[i];
  free(term_sig);
}

int doc = 0;

char *sig_buffer;     // buffer to write signatures to
long sig_buffer_size; // in bytes
long offset = 0;      // offset into sig_buffer (in bytes)

void compute_signature(char *sequence, int length, Signature *sig) {
  int *doc_sig = (int *)calloc(SIGNATURE_LEN, sizeof(int));

  for (int i = 0; i < length - WORDLEN + 1; i++) {
    signature_add(sequence + i, doc_sig);
  }

  // flatten doc_sig and write to sig
  for (int i = 0; i < SIGNATURE_LEN; i += 8) {
    byte c = 0;
    for (int j = 0; j < 8; j++)
      c |= (doc_sig[i + j] > 0) << (7 - j);
    sig->sig[i / 8] = c;
  }
  offset += ((SIGNATURE_LEN / 8) + sizeof(int));
}

int partitions = 0;
Signature *signatures;

#define min(a, b) ((a) < (b) ? (a) : (b))
void partition(char *sequence, int length, int doc) {
  int i = 0;
  do {
    signatures[partitions].doc = doc;
    compute_signature(sequence + i, min(PARTITION_SIZE, length - i),
                      &signatures[partitions]);
    i += PARTITION_SIZE / 2;
    partitions++;
  } while (i + PARTITION_SIZE / 2 < length);
}

// Process a portion of the file (from start_locs[rank] to end_locs[rank])
int process(MPI_File *fasta, const int rank, WorkInfo *work) {
  if (work->start == -1) {
    printf("Rank %d has no work\n", rank);
    return 0;
  }
  int chunk_size = work->end - work->start + 1;
  char *chunk = (char *)malloc(chunk_size * sizeof(char));
  signatures = (Signature *)malloc(work->num_parts * sizeof(Signature));
  MPI_File_read_at(*fasta, work->start, chunk, chunk_size, MPI_CHAR,
                   MPI_STATUS_IGNORE);
  int doc = work->doc_start;
  // send each line to partition
  char *line = strsep(&chunk, "\n");
  while (line != NULL) {
    // dispose of documentation line
    line = strsep(&chunk, "\n");
    if (line == NULL)
      break;
    // send to partition
    partition(line, strlen(line), doc);
    // get next line
    line = strsep(&chunk, "\n");
    doc++;
  }
  free(chunk);

  int all_offsets[size];
  // synchronise all_offsets
  MPI_Allgather(&offset, 1, MPI_INT, all_offsets, 1, MPI_INT, WORKER_COM);
  // my write offset will be the sum of all offsets before me
  // e.g if rank 0 has offset of 27, will still write at 0
  // rank 1 will then write at 27
  int write_offset = 0;
  for (int i = 0; i < rank; i++) {
    write_offset += all_offsets[i];
  }
  // write to file
  MPI_File sig_file;
  MPI_File_open(WORKER_COM, out_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &sig_file);
  MPI_File_write_at_all(sig_file, write_offset, signatures,
                        work->num_parts * sizeof(Signature), MPI_CHAR,
                        MPI_STATUS_IGNORE);
  MPI_File_close(&sig_file);
  return 0;
}
int main(int argc, char *argv[]) {
  std::chrono::time_point<std::chrono::high_resolution_clock> start;

  int err;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {
    start = std::chrono::high_resolution_clock::now();
  }

  WorkInfo work_info;

  // check input
  if (argc != 4) {
    if (rank == 0)
      fprintf(stderr,
              "Usage: mpirun -np <num_procs> %s <input_fasta> <WORD_LEN> "
              "<PARTITION_SIZE>\n",
              argv[0]);
    MPI_Finalize();
    return 1;
  }

  // get input
  char *filename = argv[1];
  WORDLEN = atoi(argv[2]);
  PARTITION_SIZE = atoi(argv[3]);

  // create output filename and delete contents if it exists
  snprintf(out_filename, 256, "%s.part%d_sigs%02d_%d_mpi", filename,
           PARTITION_SIZE, WORDLEN, SIGNATURE_LEN);
  if (rank == 0) {
    FILE *fp = fopen(out_filename, "w");
    fclose(fp);
  }

  WorkInfo *work_infos = (WorkInfo *)malloc(sizeof(WorkInfo) * size);
  if (rank == 0) {
    err = distribute(filename, work_infos);
    if (err) {
      MPI_Finalize();
      return 1;
    }
  }
  // scatter work_infos into each process so that each process has its own in
  // work_info
  MPI_Scatter(work_infos, sizeof(WorkInfo), MPI_CHAR, &work_info,
              sizeof(WorkInfo), MPI_CHAR, 0, MPI_COMM_WORLD);
  free(work_infos);

  // split into excess communicators
  int color = work_info.start == -1 ? MPI_UNDEFINED : 0;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &WORKER_COM);

  // open the file (MPI)
  MPI_File fasta;
  err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL,
                      &fasta);
  if (err) {
    fprintf(stderr, "Error opening file %s\n", filename);
    MPI_Finalize();
    return 1;
  }
  // process the file
  if (color == 0)
    process(&fasta, rank, &work_info);

  if (rank == 0) {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    printf("%s %f seconds\n", out_filename, duration.count());
  }
  MPI_Finalize();
}
