#include <stdio.h>
#include <math.h>
#include "prime.h"
#include "KmerCo.h"
#include "initBF.h"

unsigned long long size=0;
int x = 0;
int y = 0;

void dim(int p, int q)
{
	x=p;
	y=q;
}

unsigned long long **allocate()
{
    int i, j;
    printf("Allocating matrix of size %d x %d\n", x, y);

    unsigned long long **a = (unsigned long long**)malloc(x * sizeof(unsigned long long*));
    if (a == NULL) {
        perror("malloc failed for rows");
        return NULL;
    }

    for (i = 0; i < x; i++) {
        a[i] = (unsigned long long*)malloc(y * sizeof(unsigned long long));
        if (a[i] == NULL) {
            perror("malloc failed for columns");
            return NULL;
        }
        if (i % 100 == 0)   // stampa ogni 100 righe allocate
            printf("Row %d allocated\n", i);
    }

    printf("Initialization...\n");
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++)
            a[i][j] = 0;
        if (i % 100 == 0)   // stampa ogni 100 righe inizializzate
            printf("Row %d initialized\n", i);
    }

    printf("Allocated and Initialized 2DBF Successfully...\n");
    return a;
}

//Insertion function of KmerCo with only K-mer as input
void _set_(unsigned long long**a,int kmer_len,char *kmer, int k)
{
	unsigned long long h;
	int i,j,pos,loop;
	unsigned long long c; // Declare 'c' here
	for(loop=0;loop<k;loop++){
		h=murmur2(kmer,kmer_len,seed[loop]);
		i=h%x;
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		c=c>>(bc*(unsigned long int)pos);
		c=c+1UL;		
		if(c==0xFF)
			return;
		c=c<<(bc*pos); 
		a[i][j]=a[i][j]&rm[pos];
		a[i][j]=a[i][j] | c;
	}
}


//if( (c==0x1F && bc==5) || (c==0xFF && bc==8) ||  (c==0x3F && bc==6) || (c==0x7F && bc==7) ||  (c==0x1FF && bc==9)|| (c==0x3FF && bc==10) || (c==0xFFF && bc==12) || (c==0x3FFF && bc==14) ||  (c==0xFFFF && bc==16) ) //Write the condition of current experiment at the beginning for quicker time
//Query function of KmerCo with only K-mer as input
unsigned long long long_test_(unsigned long long**a,int kmer_len,char *kmer, int k)
{
	unsigned long long h;
	int i,j,pos,loop;
	unsigned long long c, count[k], min;
	for(loop=0;loop<k;loop++){
		h=murmur2(kmer,kmer_len,seed[loop]);
		i=h%x; 
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		count[loop]=c>>(bc*(unsigned long int)pos);
		if(count[loop]==0)
			return 0;
	}
	switch(k){
		case 1:
			return count[0];
		case 2:
			if(count[0]<count[1])
				return count[0];
			else
				return count[1];
		default:	
			min=count[0];
			for(loop=1;loop<k;loop++){
				if (min>count[loop])
					min=count[loop];
			}
			return min;
	}
}

//Insertion function of KmerCo with K-mer and reverse complement of K-mer as input
int _set_canonical_(unsigned long long**a,int kmer_len,char *kmer, char *rev_kmer,int k) //h(kmer)<h(rev_kmer)?return 0: return 1
{
	unsigned long long h, h1;
	int i,j,pos,loop,result=0;
	unsigned long longc;
	h=murmur2(kmer,kmer_len,seed[0]); 
	h1=murmur2(rev_kmer,kmer_len,seed[0]);
	if (h1<h){
		h=h1;
		result=1;
		kmer=rev_kmer;
	}
	for(loop=0;loop<k;loop++){
		unsigned long long c;
		i=h%x;
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		c=c>>(bc*(unsigned long int)pos);
		c=c+1UL;		
		if(c==0xFF)
			return result;
		c=c<<(bc*pos); 
		a[i][j]=a[i][j]&rm[pos];
		a[i][j]=a[i][j] | c;
		if (k>1)
			h=murmur2(kmer,kmer_len,seed[loop+1]); 
	}
	return result;	
}

//if( (c==0x1F && bc==5) || (c==0xFF && bc==8) ||   (c==0x3F && bc==6) || (c==0x7F && bc==7) ||  (c==0x1FF && bc==9)|| (c==0x3FF && bc==10)|| (c==0xFFF && bc==12) || (c==0x3FFF && bc==14) || (c==0xFFFF && bc==16)) //Write the condition of current experiment at the beginning for quicker time
//Query function of KmerCo with K-mer and reverse complement of K-mer as input
unsigned long long long_test_canonical_(unsigned long long**a,int kmer_len,char *kmer, char *rev_kmer, int k, int *result)
{
	unsigned long long h, h1;
	int i,j,pos,loop;	
	unsigned long long c, count[k], min;
	*result=0;

	h=murmur2(kmer,kmer_len,seed[0]); 
		h1=murmur2(rev_kmer,kmer_len,seed[0]);
	if (h1<h){
		h=h1;
		*result=1;
		kmer=rev_kmer;
	}		
	for(loop=0;loop<k;loop++){
		i=h%x; 
		j=h%y;
		pos=h%nc; 
		c=a[i][j]&em[pos];
		count[loop]=c>>(bc*(unsigned long int)pos);
		if(count[loop]==0)
			return 0;
		if (k>1)
			h=murmur2(kmer,kmer_len,seed[loop+1]); 
	}
	switch(k){
		case 1:
			return count[0];
		case 2:
			if(count[0]<count[1])
				return count[0];
			else
				return count[1];
		default:	
			min=count[0];
			for(loop=1;loop<k;loop++){
				if (min>count[loop])
					min=count[loop];
			}
			return min;
	}
}

void _free_(unsigned long long**a)
{
	free(a);
	printf("\nMemory freed successfully...\n");
}

unsigned long long selectPrime(unsigned long long k)
{
	unsigned long long i;
	for(i=1;i<total_prime;i++)
	{
		if(prime[i]>k)
			return i;
	}
}

double error(unsigned long long m, unsigned long long n)
{
	return pow((1-exp(-2*n/m)),2);
}
unsigned long long memory(unsigned long long n, double err)
{
	return (unsigned long int)(-(n*log(err))/pow(log(2),2));
}
unsigned long long number(unsigned long long m,double err)
{
	return (unsigned long long)(-(m*pow(log(2),2))/log(err));
}

void setDim(unsigned long long m)
{
	unsigned long long k=m/256; 
	int a,b,c,d,e,f;
	unsigned long long i;
	f=sqrt(k);
	i=selectPrime(f);
	a=prime[i+3];
	b=prime[i];
	//c=prime[i/2-3];
	//d=prime[i/2+3];
	dim(a,b);
	//dim2(c,d);
	printf("2DBF dimensions for aBF: \n%d  %d\n",a,b);
	//printf("2DBF dimensions for bBF: \n%d  %d\n",c,d);
	size=a*b*64;

}
