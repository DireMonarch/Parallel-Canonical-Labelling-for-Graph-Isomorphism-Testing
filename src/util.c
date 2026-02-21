/**
 * Copyright 2025 Jim Haslett
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "util.h"
#include <time.h>
#include <limits.h>


double wtime() {
    /**
     * wtime function used to generate a wall time, in seconds, similar to MPI_wtime()
     * 
     * Based entirely on the MPI_wtime implementation:
     * https://github.com/open-mpi/ompi/blob/main/ompi/mpi/c/wtime.c
     */

    double wtime;
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
    wtime  = (double)tp.tv_nsec/1.0e+9;
    wtime += tp.tv_sec;
    return wtime;
}


// Assumes little endian
void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;
    
    for (i = size-1; i >= 0; i--) {
        for (j = 7; j >= 0; j--) {
            byte = (b[i] >> j) & 1;
            printf(" %u", byte);
        }
    }
    // puts("");
}

void fprintUnsignedLongBits(FILE * f, unsigned long num) {
    // Determine the number of bits in a long
    int num_bits = sizeof(unsigned long) * CHAR_BIT;

    // Iterate from the most significant bit to the least significant bit
    for (int i = num_bits - 1; i >= 0; i--) {
        // Check if the i-th bit is set
        if ((num >> i) & 1) {
            fprintf(f, "1");
        } else {
            fprintf(f, "0");
        }
    }
}


void UnsignedLongtoBitsCharStr(unsigned long num, char *buff) {
    // Determine the number of bits in a long
    int num_bits = sizeof(unsigned long) * CHAR_BIT;
    int pos = 0;

    // Iterate from the most significant bit to the least significant bit
    for (int i = num_bits - 1; i >= 0; i--) {
        // Check if the i-th bit is set
        if ((num >> i) & 1) {
            buff[pos++] = '1';
        } else {
            buff[pos++] = '0';
        }
    }
    buff[pos++] = 0;
}


int calculate_numcells(int *ptn, int n) {
    int calc_numcells = 0;
    for (int i = 0; i < n; ++i) {
        if (ptn[i] < NAUTY_INFINITY) ++calc_numcells;
    }
    return calc_numcells;
}


void dump_automorphisms(int rank, set *workspace, set *fmptr, int m) {
    /** visualize automorphisms */
    int pos = 0;
    char filename[128];
    sprintf(filename, "auto_%02d.txt", rank);
    FILE *f = fopen(filename, "w");

    const int wordlen = sizeof(unsigned long) * CHAR_BIT;
    char fix_buff[m * wordlen + 1];
    char mcr_buff[m * wordlen + 1];
    while (workspace + pos < fmptr)
    {
        for(int i = 0; i < m; ++i) {
            UnsignedLongtoBitsCharStr(workspace[pos++], fix_buff + (i * wordlen));
            // printf(" ");
        }
        for(int i = 0; i < m; ++i) {
            UnsignedLongtoBitsCharStr(workspace[pos++], mcr_buff + (i * wordlen));
            // printf(" ");
        }
        fprintf(f, "MPI: Process %d: Automorphism: %s    %s\n", rank, fix_buff, mcr_buff);
    }       
    fclose(f);
}

boolean is_fix_not_ok(set *fmptr, int n, int m) {
    boolean firstone = FALSE;
    boolean zero = FALSE;
    boolean secondone = FALSE;

    int num_bits = sizeof(unsigned long) * CHAR_BIT;
    char bits[num_bits * m + 1];
    for (int i = 0; i < m; ++i) {
        UnsignedLongtoBitsCharStr(fmptr[i], bits+(num_bits*i));
    }

    /**
     * I think the rule is going to be, once we see a 1, then zeros, we see ones again, if we see another zero, that's bad
     */
    for (int i = 0; i < n; ++i) {
        if (!firstone && bits[i] == '1') firstone = TRUE;
        if (firstone && !zero && bits[i] == '0') zero = TRUE;
        if (zero && !secondone && bits[i] == '1') secondone = TRUE;
        if (secondone && bits[i] == '0') return TRUE;
    }
    return FALSE;
}

void recalculate_processes_running(int *processes_running, int num_processes) {
    processes_running[0] = 0;
    for (int i = 1; i < num_processes; ++i) processes_running[0] += processes_running[i];
}