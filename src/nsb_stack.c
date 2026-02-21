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

#include "nsb_stack.h"
#include "node.h"
#include "path.h"


/* put node on top of stack (using this makes the stack a LIFO stack) */ 
void stack_push(Stack *stack, Node *node) {
    StackNode *new_stacknode;
    DYNALLOCSTACKNODE(new_stacknode, "stack_push")
    new_stacknode->node = node;
    new_stacknode->next = stack->top;
    stack->top = new_stacknode;
    if (stack->bottom == NULL) stack->bottom = new_stacknode; /* basically, if the stack was empty, point stack->bottom to the new node, it's the bottom */
    
    stack->size++;
}


/* put node at bottom of stack (using this makes the stack a FIFO queue) */
void stack_enqueue(Stack *stack, Node *node) {
    StackNode *new_stacknode;
    DYNALLOCSTACKNODE(new_stacknode, "stack_enqueue")
    new_stacknode->node = node;
    new_stacknode->next = NULL;
    if (stack->bottom != NULL) stack->bottom->next = new_stacknode;
    stack->bottom = new_stacknode;
    if (stack->top == NULL) stack->top = new_stacknode;
    
    stack->size++;
}


Node* stack_pop(Stack *stack) {
    if (stack->top == NULL)
        return NULL;

    StackNode *top = stack->top;
    stack->top = top->next;
    Node *ret = top->node;
    FREESTACKNODE(top)
    stack->size--;
    
    return ret;
}

Node* stack_peek(Stack *stack) {
    if (stack->top == NULL)
        return NULL;

    return stack->top->node;
}

int stack_size(Stack *stack) {
    return stack->size;
}

Node* stack_peek_at(Stack *stack, int idx) {
    if (idx >= stack->size) {printf("Index out of bounds peeking at stack location %d, stack is only %d\n", idx, stack->size);  exit(0);}
    StackNode *ptr = stack->top;
    
    for (int i = 0; i < idx; ++i, ptr = ptr->next);                 /* this line pulls from top of stack */   
    // for (int i = 1; i < stack->size - idx; ++i, ptr = ptr->next);   /* this line pulls from bottom of stack */
    return ptr->node;
}


Node* stack_pop_at(Stack *stack, int idx) {

    if (idx >= stack->size) {printf("Index out of bounds popping at stack location %d, stack is only %d\n", idx, stack->size);  exit(0);}
    StackNode *ptr = stack->top;
    StackNode *parent = NULL;

    for (int i = 0; i < idx; ++i) {
        parent = ptr;
        ptr = ptr->next;
    }

    if (parent == NULL) {
        stack->top = ptr->next;
    } else {
        parent->next = ptr->next;
    }

    --stack->size;
    
    Node *curr = ptr->node;
    FREESTACKNODE(ptr);
    return curr;
}



void stack_visualize(Stack *stack, int n) {
    StackNode *curr = stack->top;
    Node *node;
    while (curr != NULL) {
        node = curr->node;

        printf("\t"); path_visualize(node->path); printf("\n");

        // node_visualize(node, n);
        // printf("\n");

        /* get next StackNode */
        curr = curr->next;
    }
}


void fstack_visualize(FILE *f, Stack *stack, int n, char *prefix) {
    StackNode *curr = stack->top;
    Node *node;
    int idx = 0;
    while (curr != NULL) {
        node = curr->node;

        fprintf(f, "%sSTACK[%d] ", prefix, idx++); fpath_visualize(f, node->path); fprintf(f,"\n");

        /* get next StackNode */
        curr = curr->next;
    }
}


// Stack * stack_pop_shallowest(Stack *stack, int factor) {
//     StackNode *bsf = NULL;
//     StackNode *bsf_parent = NULL;
//     StackNode *parent = NULL;
//     StackNode *curr = stack->top;
//     Node *node;
//     int count = 0;

//     /* create a new stack to return. */
//     Stack * ret;
//     DYNALLOCSTACK(ret, "stack_pop_shallowest");
//     ret->size = 0;
//     ret->top = NULL;

//     while (curr != NULL) {
//         node = curr->node;
//         if (!bsf || node->level < bsf->node->level) {
//             bsf = curr;
//             bsf_parent = parent;
//             count = 0;
//         }
//         ++count;
//         parent = curr;
//         curr = curr->next;
//     }

//     if (bsf == NULL) return ret;

//     int sz = bsf->node->path->sz * factor;  /* calculate how many nodes to return */
//     // int sz = (count+1)*percent_to_send;  /* calculate how many nodes to return */
//     if (sz > count) sz = count; /* make sure it's not more than we have at the best level */
//     if (sz > (stack->size+1)/2) sz = (stack->size+1)/2; /* make sure we don't give away more than half our work (rounded up) - This one is a guess*/

//     /* get the "next" stacknode after the last one we will return */
//     StackNode *end = bsf;
//     for (int i = 0; i < sz-1; ++i) end = end->next;
//     parent = end;  /* resing parent here to not create a new pointer */
//     end = end->next; /* setting end to the next node, which stays on this stack (if not NULL)*/
//     parent->next = NULL; /* set next on the last node to be returned to NULL */
    
//     /* remove the list of stacknodes from the stack */
//     if (bsf_parent) {
//         /* bsf is NOT the top of the stack */
//         bsf_parent->next = end;
//     } else {
//         /* bsf IS the top of the stack */
//         stack->top = end;
//     }
    
//     stack->size -= sz; /* adjust the stack size for the removed nodes*/

//     /* Set values for return stack.  giving up ownership of all the nodes and the stack to the calling function */

//     ret->top = bsf;
//     ret->size = sz;

//     return ret;  
// }