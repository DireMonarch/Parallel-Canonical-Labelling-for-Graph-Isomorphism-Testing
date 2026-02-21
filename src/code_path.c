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

#include "code_path.h"

/**
 * I'm sussing out how this should work based on the original implementation
 * 
 * It seems that when comparing the code paths of the current node, from the 
 * root to the code path of the canon node, that the first deviation in that path
 * is one level higher than equivelant level (eqlev).  Also, for the purpose
 * of comparing which is greater, it's important to track if that first deviation
 * is greater or less than the canon code at the same level.ACCESSX_MAX_TABLESIZE
 * 
 * This function walks the two code paths (node and canon).  Once it finds
 * values in those paths that are not equal (they deviate) it retuns the last level
 * in which the two paths were equal via the refernce perameter eqlev, and returns
 * -1 if the node value is less than the canon value, 1 if the node value is greater
 * than the canon value, or zero if the paths are equal.
 * 
 * Parameters:
 *      *node   pointer to the current node's CodePath struct
 *      *canon  pointer to the canon node's CodePath struct
 *      *eqlev  pointer to an integer used to return the level at which the two CodePaths are equivelent
 * 
 * Returns:
 *      -1 if the node value is less than the canon value, 1 if the node value is greater
 *          than the canon value, or zero if the paths are equal
 */
int codepath_find_eqlev_level(CodePath *node, CodePath *canon, int *eqlev) {
    // printf("Codepath Compare:\n");
    // codepath_visualize(node);printf("\n");
    // codepath_visualize(canon);printf("\n");
    *eqlev = 0;
    for (*eqlev = 0; *eqlev < node->sz && *eqlev < canon->sz; ++*eqlev) {
        if (node->data[*eqlev] != canon->data[*eqlev])
        {   
            --*eqlev;
            if (node->data[*eqlev+1] < canon->data[*eqlev+1]) return -1;
            return 1;
        }
    }
    --*eqlev;
    return 0; 
}


CodePath* codepath_make_child_path(CodePath *src) {
    CodePath *dst;
    DYNALLOCCODEPATH(dst, src->sz+1, "CodePath Make Child CodePath"); 

    for (int i = 0; i < src->sz; ++i) {
        dst->data[i] = src->data[i];
    }
    dst->data[dst->sz-1] = -1;
    return dst;
}

CodePath* codepath_deep_copy(CodePath *src) {
    CodePath *dst;
    DYNALLOCCODEPATH(dst, src->sz, "CodePath Deep Copy"); 

    for (int i = 0; i < src->sz; ++i) {
        dst->data[i] = src->data[i];
    }
    return dst;
}

void codepath_visualize(CodePath *path) {
    printf("[");
    for (int i = 0; i < path->sz; ++i) {
        if (i > 0) printf(", ");
        printf("%d", path->data[i]);
    }
    printf("]");    
}


CodePath* codepath_resize_to_new(CodePath *path, int new_size) {
    if (new_size <= path->sz) {
        printf("Error:  new size (%d) is not larger than current size (%d)\n", new_size, path->sz);
        exit(1);
    }

    CodePath *dst;
    DYNALLOCCODEPATH(dst, new_size, "CodePath resize"); 

    int i;
    for(i = 0; i < path->sz; ++i) {
        dst->data[i] = path->data[i];
    }
    for(;i < dst->sz; ++i) {
        dst->data[i] = 0;
    }

    FREECODEPATH(path);
    return dst;
}
