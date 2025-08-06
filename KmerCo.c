#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "prime.h"
#include "murmur.h"
#include "KmerCo.h"
#include "initBF.h"
#ifdef _WIN32
    #include <windows.h>
#else
     #include <unistd.h>
#endif
#include <fcntl.h>
#include <sys/types.h>

// Variabili globali
unsigned long long TP = 0, TN = 0, FP = 0;
int seed[8] = { 123, 456, 789, 1011, 1213, 1415, 1617, 1819 };
unsigned long long nc = 8;
unsigned long long bc = 8;
unsigned long long **aBF;

// Funzione di query KmerCo con solo K-mer come input
unsigned long long _test_(unsigned long long **a, int kmer_len, char *kmer, int k)
{
    unsigned long long h;
    int i, j, pos, loop;
    unsigned long long c, count[8], min; // 8 è il massimo valore di k previsto
    for(loop = 0; loop < k; loop++) {
        h = murmur2(kmer, kmer_len, seed[loop]);
        i = h % x;
        j = h % y;
        pos = h % nc;
        c = a[i][j] & em[pos];
        count[loop] = c >> (bc * (unsigned long long)pos);
        if(count[loop] == 0)
            return 0;
    }
    switch(k) {
        case 1:
            return count[0];
        case 2:
            return (count[0] < count[1]) ? count[0] : count[1];
        default:
            min = count[0];
            for(loop = 1; loop < k; loop++) {
                if(min > count[loop])
                    min = count[loop];
            }
            return min;
    }
}

// Funzione di query KmerCo con K-mer e reverse complement come input
unsigned long long _test_canonical_(unsigned long long **a, int kmer_len, char *kmer, char *rev_kmer, int k, int *result)
{
    unsigned long long h, h1;
    int i, j, pos, loop;
    unsigned long long c, count[8], min;
    *result = 0;

    h = murmur2(kmer, kmer_len, seed[0]);
    h1 = murmur2(rev_kmer, kmer_len, seed[0]);
    if(h1 < h) {
        h = h1;
        *result = 1;
        kmer = rev_kmer;
    }
    for(loop = 0; loop < k; loop++) {
        i = h % x;
        j = h % y;
        pos = h % nc;
        c = a[i][j] & em[pos];
        count[loop] = c >> (bc * (unsigned long long)pos);
        if(count[loop] == 0)
            return 0;
        if(k > 1)
            h = murmur2(kmer, kmer_len, seed[loop + 1]);
    }
    switch(k) {
        case 1:
            return count[0];
        case 2:
            return (count[0] < count[1]) ? count[0] : count[1];
        default:
            min = count[0];
            for(loop = 1; loop < k; loop++) {
                if(min > count[loop])
                    min = count[loop];
            }
            return min;
    }
}

// Funzione di memoria
unsigned long long memory(unsigned long long n, double err);

//#### Insertion senza file writing (Canonical) ##########
void insertion_canonical_without_filewrite(char fname[6][100], int kmer_len, int threshold, int k) {
    int result, kcount = 0, rcount = 0;
    char ch, kmer[kmer_len + 1], rev_kmer[kmer_len + 1];
    double err = 0.001;
    unsigned long long m;
    clock_t start, end;
    long long int count, fileLength, total_kmers;
    FILE *fkmer, *fres;

    printf("Insertion process\n");
    int f = open(fname[0], O_RDONLY);
    if (f == -1)
        printf("Error in opening file\n");

    char cmd[150] = "wc -c <", *s1 = " > linecount.txt";
    strcat(cmd, fname[0]);
    strcat(cmd, s1);
    system(cmd);

    fkmer = fopen("linecount.txt", "r");
    if (fkmer == NULL) {
        printf("File can't be created\n");
        exit(0);
    }
    fscanf(fkmer, "%lld", &fileLength);

    printf("File length: %lld\n", (long long)fileLength);
    total_kmers = fileLength - kmer_len;
    printf("No of kmers: %lld\n", (long long)total_kmers);
    fclose(fkmer);

    printf("File initiated!\n");
    int mem = (int)(total_kmers);
    m = memory(mem, err);
    printf("Memory initiated!\n");
    setDim(m);
    printf("Dimensions initiated!\n");
    aBF = allocate();
    printf("Filters are created!\n");

#ifdef _WIN32
    HANDLE hMap = CreateFileMapping((HANDLE)_get_osfhandle(f), NULL, PAGE_READONLY, 0, 0, NULL);
    if (hMap == NULL) {
        printf("CreateFileMapping failed: %lu\n", GetLastError());
        exit(1);
    }
    char* buff = (char*)MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, fileLength);
    if (buff == NULL) {
        printf("MapViewOfFile failed: %lu\n", GetLastError());
        CloseHandle(hMap);
        exit(1);
    }
#else
    char* buff = mmap(NULL, fileLength, PROT_READ, MAP_SHARED, f, 0);
    if (buff == MAP_FAILED) {
        perror("mmap failed");
        exit(1);
    }
#endif

    off_t bindex = 0;
    count = 1;
    kmer[kmer_len] = '\0';
    rev_kmer[kmer_len] = '\0';
    start = clock();
    for (kcount = 0, rcount = kmer_len - 1; kcount < kmer_len; kcount++, rcount--) {
        ch = buff[bindex++];
        kmer[kcount] = ch;
        switch (ch) {
            case 'A': rev_kmer[rcount] = 'T'; break;
            case 'C': rev_kmer[rcount] = 'G'; break;
            case 'G': rev_kmer[rcount] = 'C'; break;
            case 'T': rev_kmer[rcount] = 'A'; break;
            default: rev_kmer[rcount] = ch;
        }
    }
    result = _set_canonical_(aBF, kmer_len, kmer, rev_kmer, k);

    while (count < total_kmers) {
        for (kcount = 1, rcount = kmer_len - 1; kcount < kmer_len; kcount++, rcount--) {
            kmer[kcount - 1] = kmer[kcount];
            rev_kmer[rcount] = rev_kmer[rcount - 1];
        }
        ch = buff[bindex++];
        kmer[--kcount] = ch;
        switch (ch) {
            case 'A': rev_kmer[0] = 'T'; break;
            case 'C': rev_kmer[0] = 'G'; break;
            case 'G': rev_kmer[0] = 'C'; break;
            case 'T': rev_kmer[0] = 'A'; break;
            default: rev_kmer[0] = ch;
        }
        result = _set_canonical_(aBF, kmer_len, kmer, rev_kmer, k);
        count++;
    }
    end = clock();

#ifdef _WIN32
    UnmapViewOfFile(buff);
    CloseHandle(hMap);
    _close(f);
#else
    munmap(buff, fileLength);
    close(f);
#endif

    printf("Insertion complete!\n\n");
    printf("Total insertion:%llu\n", count);
    printf("Elapsed Time of insertion:%f\n\n", (double)(end - start) / CLOCKS_PER_SEC);
    double mb = 8 * 1024 * 1024.0;
    printf("\nRequired memory size in bits: %lu\n", size);
    printf("\nTotal memory size in MB: %lf\n", (double)(size) / mb);
    printf("\nRequired memory size in bits: %lu\n", m);
    printf("\nRequired memory size in MB: %lf\n", (double)(m) / mb);

    fres = fopen(fname[4], "a+");
    if (fres == NULL) {
        printf("Result File can't be created/Opened\n");
        exit(0);
    }
    fprintf(fres, "\n\n########### Results of dataset %s without write (Canonical) ############\n", fname[0]);
    fprintf(fres, "Kmer length: %d\n", kmer_len);
    fprintf(fres, "Total kmers: %llu\n", total_kmers);
    fprintf(fres, "Total insertion: %llu\n", count);
    fprintf(fres, "Elapsed Time of insertion: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    fprintf(fres, "Required memory size in bits and MB: %lu\t %lf\n", size, (double)(size) / mb);
    fprintf(fres, "Required memory size in bits: %lu\t %lf\n", m, (double)(m) / mb);
    fclose(fres);

    _free_(aBF);
}

//#### Insertion con file writing (Canonical) ##########
void insertion_canonical_with_filewrite(char fname[6][100], int kmer_len, int threshold, int k) {
    int result, kcount = 0, rcount = 0;
    char ch, kmer[kmer_len + 1], rev_kmer[kmer_len + 1];
    double err = 0.001;
    unsigned long long m;
    clock_t start, end;
    long long int count, fileLength, total_kmers;
    FILE *fkmer, *ferror, *ftrust, *fres;

    printf("Insertion process\n");
    int f = open(fname[0], O_RDONLY);
    if (f == -1)
        printf("Error in opening file\n");

    char cmd[150] = "wc -c <", *s1 = " > linecount.txt";
    strcat(cmd, fname[0]);
    strcat(cmd, s1);
    system(cmd);

    fkmer = fopen("linecount.txt", "r");
    if (fkmer == NULL) {
        printf("File can't be created\n");
        exit(0);
    }
    fscanf(fkmer, "%lld", &fileLength);

    printf("File length: %lld\n", (long long)fileLength);
    total_kmers = fileLength - kmer_len;
    printf("No of kmers: %lld\n", (long long)total_kmers);
    fclose(fkmer);

    printf("File initiated!\n");
    int mem = (int)(total_kmers);
    m = memory(mem, err);
    printf("Memory initiated!\n");
    setDim(m);
    printf("Dimensions initiated!\n");
    aBF = allocate();
    printf("Filters are created!\n");

#ifdef _WIN32
    HANDLE hMap = CreateFileMapping((HANDLE)_get_osfhandle(f), NULL, PAGE_READONLY, 0, 0, NULL);
    if (hMap == NULL) {
        printf("CreateFileMapping failed: %lu\n", GetLastError());
        exit(1);
    }
    char* buff = (char*)MapViewOfFile(hMap, FILE_MAP_READ, 0, 0, fileLength);
    if (buff == NULL) {
        printf("MapViewOfFile failed: %lu\n", GetLastError());
        CloseHandle(hMap);
        exit(1);
    }
#else
    char* buff = mmap(NULL, fileLength, PROT_READ, MAP_SHARED, f, 0);
    if (buff == MAP_FAILED) {
        perror("mmap failed");
        exit(1);
    }
#endif

    fkmer = fopen(fname[1], "w");
    if (fkmer == NULL) {
        printf("Distinct Kmer File can't be created\n");
        exit(0);
    }

    off_t bindex = 0;
    kmer[kmer_len] = '\0';
    rev_kmer[kmer_len] = '\0';
    start = clock();
    for (kcount = 0, rcount = kmer_len - 1; kcount < kmer_len; kcount++, rcount--) {
        ch = buff[bindex++];
        kmer[kcount] = ch;
        switch (ch) {
            case 'A': rev_kmer[rcount] = 'T'; break;
            case 'C': rev_kmer[rcount] = 'G'; break;
            case 'G': rev_kmer[rcount] = 'C'; break;
            case 'T': rev_kmer[rcount] = 'A'; break;
            default: rev_kmer[rcount] = ch;
        }
    }
    result = 0;
    unsigned long long test_res = _test_canonical_(aBF, kmer_len, kmer, rev_kmer, k, &result);
    if (test_res == 0) {
        if (result == 0) {
            fprintf(fkmer, "%s\n", kmer);
            _set_(aBF, kmer_len, kmer, k);
        } else {
            fprintf(fkmer, "%s\n", rev_kmer);
            _set_(aBF, kmer_len, rev_kmer, k);
        }
    } else {
        if (result == 0)
            _set_(aBF, kmer_len, kmer, k);
        else
            _set_(aBF, kmer_len, rev_kmer, k);
    }
    count = 1;

    while (count < total_kmers) {
        for (kcount = 1, rcount = kmer_len - 1; kcount < kmer_len; kcount++, rcount--) {
            kmer[kcount - 1] = kmer[kcount];
            rev_kmer[rcount] = rev_kmer[rcount - 1];
        }
        ch = buff[bindex++];
        kmer[--kcount] = ch;
        switch (ch) {
            case 'A': rev_kmer[0] = 'T'; break;
            case 'C': rev_kmer[0] = 'G'; break;
            case 'G': rev_kmer[0] = 'C'; break;
            case 'T': rev_kmer[0] = 'A'; break;
            default: rev_kmer[0] = ch;
        }
        result = 0;
        test_res = _test_canonical_(aBF, kmer_len, kmer, rev_kmer, k, &result);
        if (test_res == 0) {
            if (result == 0) {
                fprintf(fkmer, "%s\n", kmer);
                _set_(aBF, kmer_len, kmer, k);
            } else {
                fprintf(fkmer, "%s\n", rev_kmer);
                _set_(aBF, kmer_len, rev_kmer, k);
            }
        } else {
            if (result == 0)
                _set_(aBF, kmer_len, kmer, k);
            else
                _set_(aBF, kmer_len, rev_kmer, k);
        }
        count++;
    }
    end = clock();

#ifdef _WIN32
    UnmapViewOfFile(buff);
    CloseHandle(hMap);
    _close(f);
#else
    munmap(buff, fileLength);
    close(f);
#endif

    fclose(fkmer);

    printf("Insertion complete!\n\n");
    printf("Total insertion:%llu\n", count);
    printf("Elapsed Time of insertion:%f\n\n", (double)(end - start) / CLOCKS_PER_SEC);
    TP = 0; TN = 0; FP = 0;
    double mb = 8 * 1024 * 1024.0;
    printf("\nRequired memory size in bits: %lu\n", size);
    printf("\nTotal memory size in MB: %lf\n", (double)(size) / mb);
    printf("\nRequired memory size in bits: %lu\n", m);
    printf("\nRequired memory size in MB: %lf\n", (double)(m) / mb);

    fres = fopen(fname[4], "a+");
    if (fres == NULL) {
        printf("Result File can't be created/Opened\n");
        exit(0);
    }
    fprintf(fres, "\n\n########### Results of dataset %s with write (Canonical) ############\n", fname[0]);
    fprintf(fres, "Kmer length: %d\n", kmer_len);
    fprintf(fres, "Total kmers: %llu\n", total_kmers);
    fprintf(fres, "Total insertion: %llu\n", count);
    fprintf(fres, "Elapsed Time of insertion: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    fprintf(fres, "Required memory size in bits and MB: %lu\t %lf\n", size, (double)(size) / mb);
    fprintf(fres, "Required memory size in bits and MB: %lu\t %lf\n", m, (double)(m) / mb);

    //########################## Classification ##########################################
    fkmer = fopen(fname[1], "r");
    ferror = fopen(fname[2], "w");
    ftrust = fopen(fname[3], "w");
    if (fkmer == NULL || ferror == NULL || ftrust == NULL) {
        printf("Distinct file or erroneous kmers file or trustworthy kmer file can't be created\n");
        exit(0);
    }

    fscanf(fkmer, "%s", kmer);
    count = 1;
    start = clock();
    while (!feof(fkmer)) {
        if (_test_(aBF, kmer_len, kmer, k) > threshold)
            fprintf(ftrust, "%s\n", kmer);
        else
            fprintf(ferror, "%s\n", kmer);
        fscanf(fkmer, "%s", kmer);
        count++;
    }
    end = clock();
    fclose(fkmer);
    fclose(ftrust);
    fclose(ferror);

    printf("Elapsed Time of classification: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Total queries: %llu\n", count - 1);
    fprintf(fres, "\nElapsed Time of query: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    fprintf(fres, "Total queries: %llu\n", count);

    char *s = "wc -l < ";
    for (int i = 1; i < 4; i++) {
        strcpy(cmd, s);
        strcat(cmd, fname[i]);
        strcat(cmd, s1);
        system(cmd);

        ferror = fopen("linecount.txt", "r");
        if (ferror == NULL) {
            printf("File can't be created\n");
            exit(0);
        }
        fscanf(ferror, "%lld", &fileLength);
        fclose(ferror);
        if (i == 1)
            fprintf(fres, "Total number of distinct kmers: %lld\n", (long long)fileLength + 1);
        else if (i == 2)
            fprintf(fres, "Total number of erroneous kmers: %lld\n", (long long)fileLength + 1);
        else
            fprintf(fres, "Total number of trustworthy kmers: %lld\n", (long long)fileLength + 1);
    }
    fclose(fres);
}

int main(int argc, char* argv[])
{
    int kmer_len, threshold, k;
    char fname[6][100] = {"0", "Distinct.txt", "Errorneous.txt", "Trustworthy.txt", "Result.txt", ""};
    char* parameter[] = {"-K", "-eta", "-h"};

    if (argc == 2) {
        printf("Check the command format in README file!\n");
        printf("Considering the default values: K=28, threshold=5 and number of hash function: 1\n\n");
        kmer_len = 28;
        threshold = 5;
        k = 1;
        strcpy(fname[0], argv[1]);
    }
    else if (argc == 8) {
        for (int i = 1; i < 6; i = i + 2) {
            if (strcmp(argv[i], parameter[0]) == 0)
                kmer_len = atoi(argv[i + 1]);
            else if (strcmp(argv[i], parameter[1]) == 0)
                threshold = atoi(argv[i + 1]);
            else if (strcmp(argv[i], parameter[2]) == 0)
                k = atoi(argv[i + 1]);
            else {
                printf("Check the command format in README file!\n");  exit(0);
            }
        }
        strcpy(fname[0], argv[7]);
    }
    else {
        printf("Check the command format in README file!\n");
        printf("At least provide genomic dataset in fasta format!\n");
        exit(0);
    }

    insertion_canonical_without_filewrite(fname, kmer_len, threshold, k);
    insertion_canonical_with_filewrite(fname, kmer_len, threshold, k);

    return 0;
}