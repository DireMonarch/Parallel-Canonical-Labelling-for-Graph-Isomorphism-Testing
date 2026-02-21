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

#include "path.h"

// level: 4   ccode:  0 132 31388 2432 27196 32767 0 0 0 0[0, 132, 31388, 2432] code: 27196

int path_find_equiv_level_old(Path *a, Path *b) {
    int i = 0;
    for (i = 0; i < a->sz && i < b->sz; ++i) {
        if (a->data[i] != b->data[i]) return i+1;  /* while i is the first node that doesn't match, 
                                                    we don't do -1 because of the 1 indexed nature 
                                                    of nauty's levels */
    }
    return i; /* while i is the first node that doesn't match, 
                 we don't do -1 because of the 1 indexed nature 
                 of nauty's levels */
}
int path_find_equiv_level(Path *a, short *b) {
    int i = 0;
    for (i = 0; i < a->sz; ++i) {
        if (a->data[i] != b[i]) return i-1;  
    }
    return i-1; 
}


Path* path_make_child_path(Path *src) {
    Path *dst;
    DYNALLOCPATH(dst, src->sz+1, "Path Make Child Path"); 

    for (int i = 0; i < src->sz; ++i) {
        dst->data[i] = src->data[i];
    }
    dst->data[dst->sz-1] = -1;
    return dst;
}

Path* path_deep_copy(Path *src) {
    Path *dst;
    DYNALLOCPATH(dst, src->sz, "Path Deep Copy"); 

    for (int i = 0; i < src->sz; ++i) {
        dst->data[i] = src->data[i];
    }
    return dst;
}

int path_greatest_common_ancestor(Path *a, Path *b) {
    int i;
    for( i = 0; i < a->sz && i < b->sz; ++i) {
        if (a->data[i] != b->data[i]) return i + 1;  /* i+1 because level is one indexed
                                                    but data is zero indexed, and the
                                                    [0] entry is level 2, not level 1.
                                                    The array is actuall two indexed */
    }
    if (a->sz > b->sz) return b->sz + 1;
    return a->sz + 1;
}

boolean path_compare_to_level(Path *a, Path *b, int level) {
    if (a->sz < level - 1 || b->sz < level - 1) return FALSE;

    int i;
    for( i = 0; i < (level - 1); ++i) {
        if (a->data[i] != b->data[i]) return FALSE;
    }
    return TRUE;
}

void path_visualize(Path *path) {
    printf("[");
    for (int i = 0; i < path->sz; ++i) {
        if (i > 0) printf(", ");
        printf("%d", path->data[i]);
    }
    printf("]");    
}


void fpath_visualize(FILE *f, Path *path) {
    fprintf(f, "[");
    for (int i = 0; i < path->sz; ++i) {
        if (i > 0) fprintf(f, ", ");
        fprintf(f, "%d", path->data[i]);
    }
    fprintf(f, "]");    
}