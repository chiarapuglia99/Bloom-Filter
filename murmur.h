#ifndef MURMUR_H
#define MURMUR_H

unsigned int murmur2(const void *key, int len, unsigned int seed);
unsigned int murmur_backup2(const void *key, int len, unsigned int seed);

#endif