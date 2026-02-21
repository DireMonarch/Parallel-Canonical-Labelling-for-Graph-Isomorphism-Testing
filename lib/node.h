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


#ifndef _NODE_H_
#define _NODE_H_

#include <stddef.h>
#include "nauty.h"
#include "path.h"
#include "code_path.h"

typedef struct Node Node;
struct Node
{
    Path *path;         /* Path in the tree of the current node */
    int *lab, *ptn;     /* Define the partition nest of this node */
    CodePath *codes;    /* Code path from root to this node at each level.  used to calculate eqlev_canon (and possibly eqlev_first) */
    set *fixedpts;      /* points which were explicitly fixed to get current node */
    set *active;        /* used to contain index to cells now active for refinement purposes */
    set *tcell;         /* nodes active in the current target cell */
    int level;          /* Level of the tree this node is at */
    int target_cell;    /* position of target cell in lab */
    int target_vertex;  /* target vertex within the target cell that produces this child */
    int numcells;       /* Number of cells in the partition nest (lab,ptn) */
    int gca_first;      /* level of greatest common ancestor of current node and first leaf */
    int gca_canon;      /* ditto for current node and bsf leaf */
    int eqlev_first;    /* level to which codes for this node match those for first leaf */
    int eqlev_canon;    /* level to which codes for this node match those for the bsf leaf. */    
    int comp_canon;     /* -1,0,1 according as code at eqlev_canon+1 is <,==,> that for bsf leaf.  Also used for similar purpose during leaf processing */
    int cosetindex;     /* the point being fixed at level gca_first */
    int noncheaplevel;  /* level of greatest ancestor for which cheapautom==FALSE */
    boolean fpchild;    /* TRUE if node is immediate child of a firstpathnode.  i.e. generated in firstpathnode not othernode */
};
 

Node* node_make_child(Node *, set *tcell, int, int);
void node_visualize(Node *, int);
void fnode_visualize(FILE *f, Node *node, int n);
void node_serialize(Node *, int, int, int *, int *, int *, unsigned long *);
Node *node_deserialize(int n, int m, int *pos, int *buffer, int *ul_pos, unsigned long *ul_buffer);
void node_get_serialization_size(Node *, int *, int *, int, int);

#define DYNALLOCNODE(name,msg) \
    if ((name= (Node*)malloc(sizeof(Node))) == NULL) {alloc_error(msg);}; 


#define FREENODE(name) \
    if (name) { \
        FREES(name->lab); \
        FREES(name->ptn); \
        FREECODEPATH(name->codes); \
        FREEPATH(name->path); \
        FREES(name->fixedpts); \
        FREES(name->active); \
        FREES(name->tcell); \
        FREES(name); \
        name = NULL; \
    }


 #endif /* _NODE_H_ */