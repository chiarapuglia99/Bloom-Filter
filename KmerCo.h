#ifndef KMERCO_H
#define KMERCO_H

#include <stdlib.h>
#include "mask8.h"
#include "murmur.h"

extern int seed[8];
extern unsigned long long nc;
extern unsigned long long bc;


void dim(int p, int q);
unsigned long long _test_(unsigned long long **a,int kmer_len,char *kmer, int k);
unsigned long long **allocate(void);
void _set_(unsigned long long **a, int kmer_len, char *kmer, int k);
unsigned long long long_test_(unsigned long long **a, int kmer_len, char *kmer, int k);
int _set_canonical_(unsigned long long **a, int kmer_len, char *kmer, char *rev_kmer, int k);
unsigned long long long_test_canonical_(unsigned long long **a, int kmer_len, char *kmer, char *rev_kmer, int k, int *result);
void _free_(unsigned long long **a);

#endif