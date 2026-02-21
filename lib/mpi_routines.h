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

#ifndef _MPI_ROUTINES_H_
#define _MPI_ROUTINES_H_

#include "nsb_stack.h"
#include "mpi.h"
#include <stdlib.h>
#include <time.h>

#define __DEBUG_MPI__ FALSE
#define __DEBUG_MPI_PEROIDIC_STATUS__ FALSE
#define __DEBUG_MPI_PERIODIC_STATUS_UPDATE_PERIOD__ 2.0
#define __DEBUG_MPI_PRUNE__ FALSE
#define __DEBUG_MPI_SENDWORK__ FALSE 
#define __DEBUG_MPI_WORK_RETURN__ FALSE 
#define __DEBUG_MPI_NEWAUTO__ FALSE
#define __DEBUG_MPI_NEWCL__ FALSE

#define MPI_CONST_SPLIT_RATIO 2
#define MPI_CONST_MAX_WORK_SIZE_TO_SEND 10
#define MPI_CONST_NODES_BETWEEN_COMM_POLLS 50    /* How many nodes should we process between polling for new messages? */

#define MPI_STATE_WORK_END -1
#define MPI_STATE_WORKING 0
#define MPI_STATE_ASKING_FOR_WORK 10
#define MPI_STATE_WORK_RECEIVED 11
#define MPI_STATE_WAITING_FOR_START 99
#define MPI_STATE_REJECTED 101
#define MPI_STATE_WAIT_FOR_WORK_END 102
#define MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES 200
#define MPI_STATE_WAIT_FOR_WORK_RETURN 201

typedef struct MPIState MPIState;
struct MPIState{
    int my_rank;                    /* my MPI rank */
    int num_processes;              /* number of processes running */
    int state;                      /* current state machine state */
    // int partner_rank;               /* partner_rank, if needed, else leave -1 */
    // int workstop_detection_state;   /* used for Dijkstra's modified detection algorithm, use MPI_TOKEN_STATE_CLEAN/DIRTY */
	float mpi_last_peroidic_update; /* used to determine if periodic update debug messages should be sent */
};


typedef struct CanonState CanonState;
struct CanonState {
    Stack *stack;
	Path **canonpath;
	CodePath **canoncodes;
	int *canonlab;
	int *canonlevel;
	int *samerows;
	set *worktop;
	set *workspace;
	set **fmptr;
	int *queue_sizes;
	int *process_waiting_for_work;
	void (*mpi_received_active_prune)(MPIState *, Path *, int, Stack *);
	int max_work_size_to_send;
	int work_split;
	unsigned long numnodes;
	int n;
	int m;
};

#define MPI_QUEUE_SIZE_EMPTY 0
#define MPI_QUEUE_SIZE_REQUESTED_WORK -1
#define MPI_QUEUE_SIZE_UNDEFINED -2

#define MPI_MSG_ABORT -10
#define MPI_MSG_TIMEOUT -1
#define MPI_MSG_FOUND_FIRST_LEAF 10
#define MPI_MSG_NEED_WORK 100
#define MPI_MSG_TAKE_WORK 1100
#define MPI_MSG_TAKE_WORK_UL 1101
#define MPI_MSG_STACK_SIZE 1150
#define MPI_MSG_RETURN_WORK 1155
#define MPI_MSG_TAKE_RETURN 1160
#define MPI_MSG_TAKE_RETURN_UL 1161
#define MPI_MSG_REJECT_RETURN_WORK 1169
#define MPI_MSG_REJECT_NEED_WORK 1199
#define MPI_MSG_NEW_CL 200
#define MPI_MSG_NEW_AUTO 300
#define MPI_MSG_ACTIVE_PRUNE 400




void mpi_start_messages(MPIState *mpi_state, Path **firstpath, int *firstlab, CodePath **firstcodes, int *firsttc, int firstleaflevel, int *canonlevel, int *allsamelevel, Path **canonpath, int *canonlab, CodePath **canoncodes, int n);
// void mpi_send_start_message(MPIState *mpi_state, Path *firstpath, int *firstlab, CodePath *firstcodes, int *firsttc, int firstleaflevel, int canonlevel, int allsamelevel, int n);
// void mpi_send_start_abort_message(MPIState *mpi_state);
// int mpi_wait_for_start (MPIState *mpi_state, Path **firstpath, int *firstlab, CodePath **firstcodes, Path **canonpath, int *canonlab, CodePath **canoncodes, int *firsttc, int *canonlevel, int *allsamelevel, int *samerows, int n);
void mpi_poll_for_messages  (MPIState *mpi_state, CanonState *canonstate);
void mpi_report_stack_size(MPIState *mpi_state, CanonState *canonstate);
void mpi_ask_for_work       (MPIState *mpi_state, CanonState *canonstate);
void mpi_process_work_request_with_empty_queue(MPIState *mpi_state, CanonState *canonstate, int requesting_rank);
void mpi_broadcast_work_request_reject(MPIState *mpi_state, CanonState *canonstate);
void mpi_ask_for_work_return(MPIState *mpi_state, CanonState *canonstate);
void mpi_return_work(MPIState *mpi_state, CanonState *canonstate);
// void mpi_query_work_end     (MPIState *mpi_state, CanonState *canonstate);
void mpi_send_new_best_cl   (MPIState *, Path *, CodePath *, int *, int, int, int, CanonState*);
void mpi_send_new_automorphism(MPIState *mpi_state, set *fix, set *mcr, int m, CanonState*);
void mpi_send_active_prune_message(MPIState *mpi_state, Path *path, int to_level, CanonState*);
void _mpi_send_work_to(MPIState *mpi_state, CanonState *canonstate, int target_process);
void _mpi_resend_pseudo_broadcast(MPIState *mpi_state, const void *msg, int msg_sz, MPI_Datatype datatype, int tag, MPI_Comm comm, int origin, CanonState *canonstate);
void _mpi_pseudo_send(MPIState *mpi_state, const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, CanonState *canonstate);
void _mpi_pseudo_send_multiple(MPIState *mpi_state, int count, void **buf, int *buff_sz, MPI_Datatype *datatype, int dest, int *tags, MPI_Comm comm, CanonState *canonstate);
void _mpi_pseudo_send_to_many(MPIState *mpi_state, const void *buf, int buff_sz, MPI_Datatype datatype, int *dest, int dest_sz, int tag, MPI_Comm comm, CanonState *canonstate);

#endif /* _MPI_ROUTINES_H_ */
