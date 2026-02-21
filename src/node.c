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

#include "node.h"
#include "nauty.h"

Node* node_make_child(Node *src, set *tcell, int n, int m) {
    Node *dst;
    int lab_sz=0, ptn_sz=0;

    DYNALLOCNODE(dst, "node_deep_copy");
    DYNALLOC1(int,dst->lab,lab_sz,n,"node_deep_copy");
    DYNALLOC1(int,dst->ptn,ptn_sz,n,"node_deep_copy");



    dst->comp_canon = src->comp_canon;
    dst->cosetindex = src->cosetindex;
    dst->eqlev_canon = src->eqlev_canon;
    dst->eqlev_first = src->eqlev_first;
    dst->gca_canon = src->gca_canon;
    dst->gca_first = src->gca_first;
    dst->level = src->level+1;
    dst->numcells = src->numcells;
    dst->target_cell = src->target_cell;
    dst->target_vertex = src->target_vertex;
    dst->noncheaplevel = src->noncheaplevel;

    dst->path = path_make_child_path(src->path);
    dst->codes = codepath_make_child_path(src->codes);

    dst->fixedpts = (set*)malloc(m * sizeof(set));
    dst->active = (set*)malloc(m * sizeof(set));
    dst->tcell = (set*)malloc(m * sizeof(set));

    for (int i = 0; i < m; ++i) {
        dst->fixedpts[i] = src->fixedpts[i];
        dst->active[i] = src->active[i];
        dst->tcell[i] = tcell[i];
    }

    for (int i = 0; i < n; ++i) {
        dst->lab[i] = src->lab[i];
        dst->ptn[i] = src->ptn[i];
    }

    return dst;
}


void node_visualize(Node *node, int n) {

    /* print LAB */
    for (int i = 0; i < n; ++i) {
        // if (i > 0) printf(", ");
        printf("%-4d", node->lab[i]);
    }
    printf("\n");

    /* print PTN */
    for (int i = 0; i < n; ++i) {
        // if (i > 0) printf(", ");
        if (node->ptn[i] >= NAUTY_INFINITY)
            printf("-   ");
        else
            printf("%-4d", node->ptn[i]);
    }        
    printf("\nLV: %d  TC: %d  TV: %d\n", node->level, node->target_cell, node->target_vertex);
    printf("EQL_F: %d   COMP_C: %d   EQL_C: %d\n", node->eqlev_first, node->comp_canon, node->eqlev_canon);
}

void fnode_visualize(FILE *f, Node *node, int n) {

    /* print LAB */
    fprintf(f, "\tlab: ");
    for (int i = 0; i < n; ++i) {
        // if (i > 0) printf(", ");
        fprintf(f, "%-4d", node->lab[i]);
    }
    fprintf(f,"\n");

    /* print PTN */
    fprintf(f, "\tptn: ");
    for (int i = 0; i < n; ++i) {
        // if (i > 0) printf(", ");
        if (node->ptn[i] >= NAUTY_INFINITY)
            fprintf(f,"-   ");
        else
            fprintf(f,"%-4d", node->ptn[i]);
    }        
    fprintf(f,"\n");
}

void node_serialize(Node *node, int n, int m, int *pos, int *buffer, int *ul_pos, unsigned long *ul_buffer) {
    /* add the path to the buffer */
    buffer[(*pos)++] = node->path->sz;
    for (int i = 0; i < node->path->sz; ++i) {
        buffer[(*pos)++] = node->path->data[i];
    }

    /* add lab to the buffer */
    for (int i = 0; i < n; ++i) {
        buffer[(*pos)++] = node->lab[i];
    }

    /* add ptn to the buffer */
    for (int i = 0; i < n; ++i) {
        buffer[(*pos)++] = node->ptn[i];
    }

    /* add codes to the buffer */
    buffer[(*pos)++] = (int)(node->codes->sz);
    for (int i = 0; i < node->codes->sz; ++i) {
        buffer[(*pos)++] = (int)(node->codes->data[i]);
    }


    /* add the rest of the node data */
    buffer[(*pos)++] = node->level;          
    buffer[(*pos)++] = node->target_cell;    
    buffer[(*pos)++] = node->target_vertex;  
    buffer[(*pos)++] = node->numcells;       
    buffer[(*pos)++] = node->gca_first;      
    buffer[(*pos)++] = node->gca_canon;                      
    buffer[(*pos)++] = node->eqlev_first;    
    buffer[(*pos)++] = node->eqlev_canon;       
    buffer[(*pos)++] = node->comp_canon;     
    buffer[(*pos)++] = node->cosetindex;     
    buffer[(*pos)++] = node->noncheaplevel;  
    buffer[(*pos)++] = node->fpchild;        

    /** Handle the unsigned long stuff here */
    /* add fixedpts to the buffer */
    for (int i = 0; i < m; ++i) {
        ul_buffer[(*ul_pos)++] = node->fixedpts[i];
    }
    /* add active to the buffer */
    for (int i = 0; i < m; ++i) {
        ul_buffer[(*ul_pos)++] = node->active[i];
    }
    /* add tcell to the buffer */
    for (int i = 0; i < m; ++i) {
        ul_buffer[(*ul_pos)++] = node->tcell[i];
    }
    /**  */
}

Node *node_deserialize(int n, int m, int *pos, int *buffer, int *ul_pos, unsigned long *ul_buffer) {
    Node *node;
    DYNALLOCNODE(node, "node_deserialize");


    /* deserialize path */
    DYNALLOCPATH(node->path, buffer[*pos], "node_deserialize");
    ++(*pos);  // move past the path size
    for (int i = 0; i < node->path->sz; ++i) {
        node->path->data[i] = buffer[(*pos)++];
    }

    int lab_sz = 0, ptn_sz = 0;

    /* deserialize lab */
    DYNALLOC1(int, node->lab, lab_sz, n, "node_deserialize");
    for (int i = 0; i < n; ++i) {
        node->lab[i] = buffer[(*pos)++];
    }

    /* deserialize ptn */
    DYNALLOC1(int, node->ptn, ptn_sz, n, "node_deserialize");
    for (int i = 0; i < n; ++i) {
        node->ptn[i] = buffer[(*pos)++];
    }

    /* deserialize codes */
    DYNALLOCCODEPATH(node->codes, buffer[*pos], "node_deserialize");
    ++(*pos);  // move past the codes size
    for (int i = 0; i < node->codes->sz; ++i) {
        node->codes->data[i] = buffer[(*pos)++];
    }

    /* deserialize the rest of the node data */
    node->level = buffer[(*pos)++];
    node->target_cell = buffer[(*pos)++];
    node->target_vertex = buffer[(*pos)++];
    node->numcells = buffer[(*pos)++];
    node->gca_first = buffer[(*pos)++];
    node->gca_canon = buffer[(*pos)++];
    node->eqlev_first = buffer[(*pos)++];
    node->eqlev_canon = buffer[(*pos)++];
    node->comp_canon = buffer[(*pos)++];
    node->cosetindex = buffer[(*pos)++];
    node->noncheaplevel = buffer[(*pos)++];
    node->fpchild = buffer[(*pos)++] != 0;


    /** Deserialize the unsigned long items */
    /* deserialize fixedpts */
    node->fixedpts = (set*)malloc(m * sizeof(set));
    for (int i = 0; i < m; ++i) {
        node->fixedpts[i] = ul_buffer[(*ul_pos)++];
    }

    /* deserialize active */
    node->active = (set*)malloc(m * sizeof(set));
    for (int i = 0; i < m; ++i) {
        node->active[i] = ul_buffer[(*ul_pos)++];
    }

    /* deserialize tcell */
    node->tcell = (set*)malloc(m * sizeof(set));
    for (int i = 0; i < m; ++i) {
        node->tcell[i] = ul_buffer[(*ul_pos)++];
    }    
    
    return node;
}

    void node_get_serialization_size(Node *node, int *sz, int *ul_sz, int n, int m) {

    *sz += 1;  /* for the path size */
    *sz += node->path->sz;  /* for the path data */

    *sz += n;  /* for the lab */
    *sz += n;  /* for the ptn */

    *sz += 1;  /* for the codes size */
    *sz += node->codes->sz;  /* for the codes data */

    *ul_sz += m;  /* for fixedpts */
    *ul_sz += m;  /* for active */
    *ul_sz += m;  /* for tcell */

    *sz += 12; /* for the rest of the node data */
}
