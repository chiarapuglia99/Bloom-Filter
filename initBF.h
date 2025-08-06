#ifndef INITBF_H
#define INITBF_H

extern int x;
extern int y;

extern unsigned long long size;
void dim(int, int);
unsigned long long **allocate();
void _free_(unsigned long long **);


unsigned long long selectPrime(unsigned long long k);
double error(unsigned long long m, unsigned long long n);
unsigned long long memory(unsigned long long n, double err);
unsigned long long number(unsigned long long m, double err);
void setDim(unsigned long long m);

#endif