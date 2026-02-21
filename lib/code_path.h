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

#ifndef _CODE_PATH_H_
#define _CODE_PATH_H_

#include <stddef.h>
#include "nauty.h"

typedef struct {
    short *data;  /* array of values in the path */
    int sz;     /* size of the path array */
} CodePath;

int codepath_find_eqlev_level(CodePath *node, CodePath *canon, int *eqlev);
CodePath* codepath_make_child_path(CodePath *);
CodePath* codepath_deep_copy(CodePath *);
void codepath_visualize(CodePath *);
CodePath* codepath_resize_to_new(CodePath *path, int new_size);

#define DYNALLOCCODEPATH(name, size, msg) \
    if ((name= (CodePath*)malloc(sizeof(CodePath))) == NULL) {alloc_error(msg);}; \
    if ((name->data = (short*)malloc((size) * sizeof(short))) == NULL) {alloc_error(msg);}; \
    name->sz = size;


#define FREECODEPATH(name) \
    if (name) { \
        FREES(name->data); \
        FREES(name); \
    } \
    name = NULL; 

#endif /* _CODE_PATH_H_ */