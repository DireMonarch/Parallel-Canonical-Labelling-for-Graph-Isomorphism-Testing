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

#include "nauty.h"
#include <stddef.h>
#include <stdio.h>
#ifndef _UTIL_H_
#define _UTIL_H_

double wtime();
void printBits(size_t const size, void const * const ptr);
void fprintUnsignedLongBits(FILE * f, unsigned long num);
void UnsignedLongtoBitsCharStr(unsigned long num, char *buff);
int calculate_numcells(int *ptn, int n);
void dump_automorphisms(int rank, set *workspace, set *fmptr, int m);
boolean is_fix_not_ok(set *fmptr, int n, int m);
void recalculate_processes_running(int *processes_running, int num_processes);

#endif /* _UTIL_H_ */