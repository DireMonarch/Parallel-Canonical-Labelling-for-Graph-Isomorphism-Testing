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

#include "mpi_routines.h"
#include "path.h"
#include "code_path.h"
#include "util.h"


//TODO: remove all referencs to this, it is for debug only
/**
 * This function simply prints out the current queue_sizes array for debugging.
 */
void _debug_queue_size(MPIState *mpi_state, CanonState *canonstate) {
    printf("MPI: Process %d: queue_sizes:", mpi_state->my_rank);
    for (int i = 0; i < mpi_state->num_processes; ++i) {
        printf("  %d", canonstate->queue_sizes[i]);
    }
    printf("\n");
}

/**
 * Need to send the following details:
 * 
 *      firstpath -> firstpath & canonpath
 *      firstlab -> firstlab & canonlab
 *      firstcodes[:canoncodes_sz+1] -> firstcodes & canoncodes[:sz-1]
 *      firsttc
 *      stglb_canonlevel
 *      stglb_allsamelevel
 *      
 */


void mpi_start_messages(MPIState *mpi_state, Path **firstpath, int *firstlab, CodePath **firstcodes, int *firsttc, int firstleaflevel, 
    int *canonlevel, int *allsamelevel, Path **canonpath, int *canonlab, CodePath **canoncodes, int n) {

    if (__DEBUG_MPI__) printf("MPI: Process %d: in mpi_start_messages\n",mpi_state->my_rank);
    int singles[4];
    if (mpi_state->my_rank == 0) {
        singles[0] = (*firstpath)->sz;
        singles[1] = firstleaflevel+2;
        singles[2] = *canonlevel;
        singles[3] = *allsamelevel;
    }

    /* Originally tried this wit MPI_Send, changed to MPI_Bcast.  Should be more efficient, and synchronize all processes at the start.
        This is the only true synchronied section */

    /* singles contins the lengths and single integer elements */
    MPI_Bcast(&singles, 4, MPI_INT, 0, MPI_COMM_WORLD);
    
    /* broadcast firstpath */
    if (mpi_state->my_rank != 0) {
        FREEPATH((*firstpath));
        DYNALLOCPATH((*firstpath), (singles[0]), "mpi_start_messages alloc new firstpath");
    }
    MPI_Bcast((*firstpath)->data, singles[0], MPI_INT, 0, MPI_COMM_WORLD);

    /* broadcast firstlab */
    MPI_Bcast(firstlab, n, MPI_INT, 0, MPI_COMM_WORLD);

    /* broadcast firstcodes */
    if (mpi_state->my_rank != 0) {
        FREECODEPATH((*firstcodes));
        DYNALLOCCODEPATH((*firstcodes), ((*firstpath)->sz+2), "mpi_start_messages alloc new firstcodes")// firstcodes is 2 larger than firstpath
    }
    MPI_Bcast((*firstcodes)->data, (*firstpath)->sz+2, MPI_SHORT, 0, MPI_COMM_WORLD);

    /* broadcast firsttc */
    MPI_Bcast(firsttc, singles[1], MPI_INT, 0, MPI_COMM_WORLD);

    if (mpi_state->my_rank != 0) {
        *canonlevel = singles[2];
        *allsamelevel = singles[3];


        /* replicate necessary data */
        for (int i = 0; i < n; ++i) {
            canonlab[i] = firstlab[i];
        }

        FREEPATH((*canonpath));
        DYNALLOCPATH((*canonpath), ((*firstpath)->sz), "mpi_start_messages alloc new firstpath");
        for (int i = 0; i < (*firstpath)->sz; ++i) {
            (*canonpath)->data[i] = (*firstpath)->data[i];
        }
        
        FREECODEPATH((*canoncodes));
        DYNALLOCCODEPATH((*canoncodes), ((*firstcodes)->sz-1), "mpi_start_messages alloc new firstcodes")// canoncodes is one smaller than firstcodes
        for (int i = 0; i < (*canoncodes)->sz; ++i) {
            (*canoncodes)->data[i] = (*firstcodes)->data[i];
        }        
    }
    if (__DEBUG_MPI__) printf("MPI: Process %d: exit mpi_start_messages\n",mpi_state->my_rank);

}



/**
 * This function polls for general messages coming from other processes
 * 
 * It does not handle specific messages, like sending work between processes, or
 * work end token state communications (other than responding to other processes tokens)
 * 
 * It is important that this function doesn't change the mpi_state->state variable.  That is used
 * by the calling funciton to keep track of what state the current process is actually in.  The
 * messages handled by this function should not chagne that state.
 * 
 * One potential exception to this is the work end process, will likely send a broadcast to 
 * stop all work, that might change to the final work end state.
 */
void mpi_poll_for_messages (MPIState *mpi_state, CanonState *canonstate) {

    if (__DEBUG_MPI_PEROIDIC_STATUS__ && mpi_state->my_rank == 0){
        /**
         * This is a periodic update for debugging purposes
         */
        float wtime = MPI_Wtime();
        if( wtime > mpi_state->mpi_last_peroidic_update + __DEBUG_MPI_PERIODIC_STATUS_UPDATE_PERIOD__) {
            printf("%011f MPI: Process %d: PSU: [0] %d, -1", wtime, mpi_state->my_rank, stack_size(canonstate->stack));
            for(int i = 1; i < mpi_state->num_processes; ++i) printf("   [%d] %d, %d", i, canonstate->queue_sizes[i], canonstate->process_waiting_for_work[i]);
            printf("\n");
            mpi_state->mpi_last_peroidic_update = wtime;
        }
    }  

    MPI_Status recv_status;  /* MPI_Recv status variable */
    int flag = 0;           /* flag used for MPI_Iprobe to report if there are messages */
    int nomsg;          /* dummy variable used when we don't care about the incoming message */
    int msg_count = 0;
    MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD , &flag , &recv_status);

    while (flag) { /* we have a message waiting */
        ++msg_count;
        switch (recv_status.MPI_TAG) /* what we do depends on what type of message this is, and our state */
        {
            
        /**
         * Rewriting this to send every canonstate->work_split item, instead of the top of the stack.
         */
        case MPI_MSG_NEED_WORK: {
            MPI_Recv(&nomsg, 1, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);
            
            /* if only rank 0 can give out work, no other rank should be here! */
            if (mpi_state->my_rank != 0) {
                printf("MPI: Process %d: Received Need Work message, but only Process Zero should!\n", mpi_state->my_rank);
                exit(MPI_MSG_NEED_WORK);
            }
            
            if (canonstate->queue_sizes[recv_status.MPI_SOURCE] == MPI_QUEUE_SIZE_REQUESTED_WORK) {
                printf("MPI: Process %d: Received Need Work message from %d, but that process is already marked as waiting!\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
                exit(MPI_MSG_NEED_WORK);
            }
            
                if (stack_size(canonstate->stack) > 0) { 
                    _mpi_send_work_to(mpi_state, canonstate, recv_status.MPI_SOURCE);

                    /**
                     * commented this chunk out, and moved it to _mpi_send_work_to so it can be reused.
                    // int start_stack_sz = stack_size(canonstate->stack);

                    // /* we have enough work to send */
                    // int send_sz = (stack_size(canonstate->stack) + (canonstate->work_split - 1)) / canonstate->work_split;  /* send amount of work, rounded up (hence the + (canonstate->work_split -1)) */
                    // if (send_sz > canonstate->max_work_size_to_send) send_sz = canonstate->max_work_size_to_send;  /* limit amount of work to send in one chunk */
                    // Node *curr;
                    // int buff_sz = 1;  /* start with 1 for the record count */
                    // int ul_buff_sz = 0; /* used for unsigned long serialization */

                    // Stack *temp_stack;
                    // DYNALLOCSTACK(temp_stack, "temp_stack in MPI_MSG_NEED_WORK");

                    // /**
                    //  * Walk through the stack, from the bottom up, popping out every canonstate->work_split value.
                    //  * Push each of these onto temp_stack, using enqueue, to preserve FIFO order (bottom of original stack items go on top of temp_stack, so they are sent first, and end up on bottom of destination stack
                    //  * While doing this loop, we can also calculate the buffer size needed for sending these nodes.
                    //  */
                    // for (int i = send_sz - 1; i >= 0; --i) {
                    //     curr = stack_pop_at(canonstate->stack, i * canonstate->work_split); /* pop every canonstate->work_split node from the stack, in revers order, up to send_sz */
                    //     stack_enqueue(temp_stack, curr);  /* using enqueue here, so temp_stack works like a FIFO queue, this preserves the order */
                    //     node_get_serialization_size(curr, &buff_sz, &ul_buff_sz, canonstate->n, canonstate->m); /* get the size of the node at stack index i */
                    // }

                    // int *msg = (int*)malloc(sizeof(int)*buff_sz);   /* allocate message buffer */
                    // unsigned long *ul_msg = (unsigned long*)malloc(sizeof(unsigned long)*ul_buff_sz); /* allocate unsigned long message buffer */

                    
                    // /* fill the message buffer */
                    // int idx = 0;                                          /* variable to be used as the message index */
                    // int ul_idx = 0;                                       /* variable to be used as the unsigned long message index */
                    // msg[idx++] = send_sz;                                 /* first word of the message is the number of nodes sent */
                    
                    // /* remove nodes from temp_stack and add to messages */
                    // for (int i = 0; i < send_sz; ++i) {
                    //     curr = stack_pop(temp_stack);                        
                    //     node_serialize(curr, canonstate->n, canonstate->m, &idx, msg, &ul_idx, ul_msg);  /* serialize the node into the message buffer */
                    //     FREENODE(curr);                         /* free the node, as we don't need it anymore */
                    // }

                    // // set up arrays to pass to send_multiple
                    // void *buffers[2] = {msg, ul_msg};
                    // int sizes[2] = {buff_sz, ul_buff_sz};
                    // MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};
                    // int tags[2] = {MPI_MSG_TAKE_WORK, MPI_MSG_TAKE_WORK_UL};
                    
                    // if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: Sending %d nodes of %d (%d) on stack to %d in NEED_WORK [%zu]\n", mpi_state->my_rank, send_sz,start_stack_sz, stack_size(canonstate->stack), recv_status.MPI_SOURCE, canonstate->numnodes);
                    // _mpi_pseudo_send_multiple(mpi_state, 2, buffers, sizes, types, recv_status.MPI_SOURCE, tags, MPI_COMM_WORLD, canonstate);

                    // /* free message buffer */
                    // FREES(msg);
                    // FREES(ul_msg);


                } else {  //if (stack_size(canonstate->stack) > 0)
                    /* we are here because we were asked for work, but don't have enough to give */

                    /* Mark process recv_status.MPI_SOURCE as requested work */
                    canonstate->process_waiting_for_work[recv_status.MPI_SOURCE] = TRUE;
                    canonstate->queue_sizes[recv_status.MPI_SOURCE] = 0;
                    
                }  //if (stack_size(canonstate->stack) > 0)
                break;
            }

        case MPI_MSG_NEW_CL: {
            /* we receieved a new best canonical label */

            int msg_sz;
            MPI_Get_count(&recv_status, MPI_INT, &msg_sz);
            int *msg = (int*)malloc(sizeof(int)*msg_sz);    /* allocate space for message buffer */

            MPI_Recv(msg, msg_sz, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);

            int pos = 0; /* variable used to walk through the message */

            Path *tmp_canonpath;
            DYNALLOCPATH(tmp_canonpath, msg[pos], "Path MPI_MSG_NEW_CL")    /* Allocate space for path */
            ++pos;    /* need to pull this out ofhte DYNALLOCPATH statement, as it would increment more than once */

            /* extract path from message */
            for (int i = 0; i < tmp_canonpath->sz; ++i) {
                tmp_canonpath->data[i] = msg[pos++];
            }

            CodePath *tmp_canoncodes;
            DYNALLOCCODEPATH(tmp_canoncodes, msg[pos], "CodePath MPI_MSG_NEW_CL")    /* Allocate space for code path */
            ++pos;    /* need to pull this out ofhte DYNALLOCPATH statement, as it would increment more than once */

            /* extract code path from message */
            for (int i = 0; i < tmp_canoncodes->sz; ++i) {
                tmp_canoncodes->data[i] = msg[pos++];
            }
            int none;
            int compcanon = codepath_find_eqlev_level(tmp_canoncodes, (*(canonstate->canoncodes)), &none); /* find the equivalence level of the new code path */

            if (compcanon < 1) {
                if (__DEBUG_MPI_NEWCL__) {printf("MPI: Process %d: Ignoring New Best CL message from %d   compcanon: %d", mpi_state->my_rank, recv_status.MPI_SOURCE, compcanon); path_visualize(tmp_canonpath); printf("\n");}
                FREES(tmp_canonpath); /* free the temporary canonical path */
                FREES(tmp_canoncodes); /* free the temporary canonical codes */
            } else {

                


                //IFF this is a better canonical label
                FREEPATH((*(canonstate->canonpath)));
                (*(canonstate->canonpath)) = tmp_canonpath; /* update the canonical path */

                FREECODEPATH((*(canonstate->canoncodes)));
                (*(canonstate->canoncodes)) = tmp_canoncodes; /* update the canonical codes */    

                /* extract partition lab form message */
                for (int i = 0; i < canonstate->n; ++i) {
                    canonstate->canonlab[i] = msg[pos++];
                }
                
                /* extract canonlevel and samerows */
                *(canonstate->canonlevel) = msg[pos++];
                *(canonstate->samerows) = msg[pos++];

                if (__DEBUG_MPI_NEWCL__) {printf("MPI: Process %d: Received New Best CL in %d words from %d  compcanon: %d  ", mpi_state->my_rank, msg_sz, recv_status.MPI_SOURCE, compcanon); path_visualize(*(canonstate->canonpath)); printf("  "); printf("\n");}
                
                /**
                 * All processes other than zero only send to zero, process zero forwards to the rest.
                 */
                if (mpi_state->my_rank == 0) {
                    _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_NEW_CL, MPI_COMM_WORLD, recv_status.MPI_SOURCE, canonstate);
                }
 
                // mpi_handle_new_best_cononical_label(status, path, pi);  /* seems we don't need this at all, as we are updating the details directly?? */
            }
            
            /* free memory */
            FREES(msg);

            break;
            }

        case MPI_MSG_NEW_AUTO: {
            /* we receieved a new automorphism */

            int msg_sz;
            MPI_Get_count(&recv_status, MPI_UNSIGNED_LONG, &msg_sz);
            unsigned long *msg = (unsigned long*)malloc(sizeof(unsigned long)*msg_sz);    /* allocate space for message buffer */

            if (msg_sz != 2 * canonstate->m) {
                printf("MPI: Process %d: ERROR: Received message size %d, expected %d in MPI_MSG_NEW_AUTO\n", mpi_state->my_rank, msg_sz, 2 * canonstate->m);
                exit(MPI_MSG_NEW_AUTO);
            }

            MPI_Recv(msg, msg_sz, MPI_UNSIGNED_LONG, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);
            if (__DEBUG_MPI_NEWAUTO__) {
                printf("MPI: Process %d: received %d words from %d in MPI_MSG_NEW_AUTO  workspace avail: %lu\n",mpi_state->my_rank, msg_sz, recv_status.MPI_SOURCE, (canonstate->worktop - *(canonstate->fmptr)) / (2*canonstate->m));
                // for (int i = 0; i < msg_sz; ++i) {
                //     printf("%zu ", msg[i]);
                // }
                // printf("\n");
            }



            int pos = 0; /* variable used to walk through the message */

            if ((*(canonstate->fmptr)) == canonstate->worktop) (*(canonstate->fmptr)) -= 2 * canonstate->m; /* the algorithm replaces the last automorphism, if it is out of space to store more automorphisms */
            /* extract fix and mcr from message */
            for (int i = 0; i < canonstate->m*2; ++i) {
                (*(canonstate->fmptr))[i] = msg[pos++];
            }


            (*(canonstate->fmptr)) += canonstate->m * 2; /* move fmptr to the next location */

            /**
             * All processes other than zero only send to zero, process zero forwards to the rest.
             */
            if (mpi_state->my_rank == 0) {
                    _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_UNSIGNED_LONG, MPI_MSG_NEW_AUTO, MPI_COMM_WORLD, recv_status.MPI_SOURCE,canonstate);
            }

            /* pass ownership of pi (the automorphism) to the main function, don't free it here! */
            // mpi_handle_new_automorphism(status, pi); /* seems we don't need this at all, as we are updating the details directly?? */
            
            /* free memory */
            free(msg);
            break;
            }            

        case MPI_MSG_ACTIVE_PRUNE: {
            /* we received an active prune message, this is a message that tells us to prune the path at a given level */

            int msg_sz;
            MPI_Get_count(&recv_status, MPI_INT, &msg_sz);
            int *msg = (int*)malloc(sizeof(int)*msg_sz);    /* allocate space for message buffer */

            MPI_Status newstat;
            newstat.MPI_ERROR = 0;
            MPI_Recv(msg, msg_sz, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &newstat);

            int to_level;
            Path *path;
            int idx = 0; /* variable used to walk through the message */
            
            if (msg[idx] < 1) {
                printf("MPI: Process %d: received ACTIVE_PRUNE from %d with path len %d. ABORTING!\n", mpi_state->my_rank, recv_status.MPI_SOURCE, msg[idx]);
                exit(MPI_MSG_ACTIVE_PRUNE);
            }
            DYNALLOCPATH(path, msg[idx], "Path MPI_MSG_ACTIVE_PRUNE"); /* Allocate space for path */
            ++idx; /* need to pull this out of the DYNALLOCPATH statement, as it would increment more than once */

            /* extract path from message */
            for (int i = 0; i < path->sz; ++i) {
                path->data[i] = msg[idx++];
            } 
            /* extract to_level from message */
            to_level = msg[idx++];
            int stack_start = stack_size(canonstate->stack);
            (*(canonstate->mpi_received_active_prune))(mpi_state, path, to_level, canonstate->stack); /* call the active prune function with the path and level */

            if (__DEBUG_MPI_PRUNE__ && stack_size(canonstate->stack) != stack_start) {
                printf("MPI: Process %d: received active prune from %d to level %d  %d -> %d\n", mpi_state->my_rank, recv_status.MPI_SOURCE, to_level, stack_start, stack_size(canonstate->stack));
            }

            /**
             * All processes other than zero only send to zero, process zero forwards to the rest.
             */
            if (mpi_state->my_rank == 0) {
                _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_ACTIVE_PRUNE, MPI_COMM_WORLD, recv_status.MPI_SOURCE,canonstate);
            }
            /* free memory */
            FREES(msg);
            FREEPATH(path); 

            break;
        }

// TODO:  This needs to change for new process    Does it???  what was I thinking here?      
        case MPI_MSG_REJECT_NEED_WORK: {
            if (mpi_state->state == MPI_STATE_ASKING_FOR_WORK) {
                return;
            }
            /* probably can ignore this, as it does seem maybe we've moved on??? */
            printf("MPI: Process %d: received REJECT WORK message from %d at an unexpected time\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
            exit(MPI_MSG_REJECT_NEED_WORK);
            break;
            }

        case MPI_MSG_STACK_SIZE: {
            /**
             * This message is a status update to process 0 with the current stack size of the sending process
             */
            if(mpi_state->my_rank != 0) {
                printf("MPI: Process %d: Received MPI_MSG_STACK_SIZE from %d, but only process 0 should receive these messages!\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
                exit(MPI_MSG_STACK_SIZE);
            }
            int msg;
            MPI_Recv(&msg, 1, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);
            canonstate->queue_sizes[recv_status.MPI_SOURCE] = msg;

            // if (mpi_state->state == MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES) {
            //     /* Only process the queue size report if we are still waiting for queue sizes, otherwise, it's from a process that already ran out of work and
            //             asked for more, and this request can be dropped on the floor */
                        
            //     if (canonstate->queue_sizes[recv_status.MPI_SOURCE] != MPI_QUEUE_SIZE_REQUESTED_WORK) {
            //         canonstate->queue_sizes[recv_status.MPI_SOURCE] = msg;
            //     }

            //     /* check to see if all processes have reported back in */
            //     int largest = MPI_QUEUE_SIZE_UNDEFINED;
            //     int largest_idx = MPI_QUEUE_SIZE_UNDEFINED;
            //     for (int i = 0; i < mpi_state->num_processes; ++i){
            //         if (canonstate->queue_sizes[i] == MPI_QUEUE_SIZE_UNDEFINED) {
            //             largest = MPI_QUEUE_SIZE_UNDEFINED;
            //             largest_idx = MPI_QUEUE_SIZE_UNDEFINED;
            //             break;
            //         }
            //         if (canonstate->queue_sizes[i] > largest) {
            //             largest = canonstate->queue_sizes[i];
            //             largest_idx = i;
            //         }
            //     }
            //     // printf("MPI: Process %d: received MPI_MSG_CURRENT_QUEUE_SIZE message from %d: %d   ", mpi_state->my_rank, recv_status.MPI_SOURCE, msg);
            //     // printf("state: %d  queue_sizes:", mpi_state->state); for (int i = 0; i < mpi_state->num_processes; ++i) printf("  %d", canonstate->queue_sizes[i]); printf("\n");
            //     if (largest != MPI_QUEUE_SIZE_UNDEFINED) {
            //         if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: All processes reported back, largest is %d with %d work available\n", mpi_state->my_rank, largest_idx, largest);
            //         if (__DEBUG_MPI_SENDWORK__) _debug_queue_size(mpi_state, canonstate);
            //         mpi_ask_for_work_return(mpi_state, canonstate, largest_idx);
            //     }
            // }
            break;
        }

        case MPI_MSG_RETURN_WORK: {
            int msg;
            MPI_Recv(&msg, 1, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);

            if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: received a work return request from %d\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
            mpi_return_work(mpi_state, canonstate);
            break;
            }

        case MPI_MSG_TAKE_RETURN: {
            if (mpi_state->my_rank != 0){
                printf("MPI: Process %d: Received MPI_MSG_TAKE_RETURN from %d, but only process 0 should!\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
                exit(MPI_MSG_TAKE_RETURN);
            }
            /* if we get this we should return to let the proper function handle it */
            return;
        }
        case MPI_MSG_REJECT_RETURN_WORK: {
            /* this function shouldn't process this message, just return and let the proper function handle it */
            return;
        }

        default: {
            if (mpi_state->state == MPI_STATE_ASKING_FOR_WORK && (recv_status.MPI_TAG == MPI_MSG_REJECT_NEED_WORK || recv_status.MPI_TAG ==MPI_MSG_TAKE_WORK)) {
                return;
            }
            /* if we get an unknown message type, print a message an quit */
            printf("Unkown message type received.  My Rank: %d   Sender Rank: %d    Tag: %d\n", mpi_state->my_rank, recv_status.MPI_SOURCE, recv_status.MPI_TAG);
            exit(99);
            break;
            }
        }

        /* check for another message, we want to handle all of them */
        MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD , &flag , &recv_status);
    }
}


void mpi_report_stack_size(MPIState *mpi_state, CanonState *canonstate) {
    /**
     * This function sends a report of the current stack size to process 0
     */
    int msg[1];
    msg[0] = stack_size(canonstate->stack);
    MPI_Request request;
    MPI_Isend(&msg, 1, MPI_INT, 0, MPI_MSG_STACK_SIZE, MPI_COMM_WORLD, &request);
}


void mpi_ask_for_work(MPIState *mpi_state, CanonState *canonstate) {
    /** In Ask for Work state */

    MPI_Status recv_status;  /* MPI status used for probe and receive functions */


    mpi_state->state = MPI_STATE_ASKING_FOR_WORK;   /* set our state machine to asking for work */
    
    int nomsg = MPI_MSG_NEED_WORK;
    if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: Asking for more work [%zu]\n", mpi_state->my_rank, canonstate->numnodes);
    _mpi_pseudo_send(mpi_state, &nomsg, 1, MPI_INT, 0, MPI_MSG_NEED_WORK, MPI_COMM_WORLD, canonstate); /* send request for work */
    int flag;

    /* this while loop allows us to keep probing for status to our work request, even if we get another message while waiting */
    while(mpi_state->state == MPI_STATE_ASKING_FOR_WORK) {
        /** wait for a response */
        flag = 0;

        MPI_Iprobe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &recv_status);
        if (flag) {
            /* depending on the tag received take action */
            switch (recv_status.MPI_TAG)
            {
            case MPI_MSG_REJECT_NEED_WORK:
                /* we got a reject, need to clear the message, and try again */
                MPI_Recv(&nomsg, 1, MPI_INT, 0, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);
                mpi_state->state = MPI_STATE_REJECTED;
                break;
                
            case MPI_MSG_TAKE_WORK: 
                ; /* fake statement to allow for a declaration to come next */
                /* we got some work */
                int start_stack_sz = stack_size(canonstate->stack); /* used for debugging, to see how much work we had before we received more */
                mpi_state->state = MPI_STATE_WORK_RECEIVED; /* set state to work received */

                int msg_sz;
                MPI_Get_count(&recv_status, MPI_INT, &msg_sz);
                int *msg = (int*)malloc(sizeof(int)*msg_sz);    /* allocate space for message buffer */

                MPI_Recv(msg, msg_sz, MPI_INT, recv_status.MPI_SOURCE, MPI_MSG_TAKE_WORK, MPI_COMM_WORLD, &recv_status);

                int ul_msg_sz = msg[0] * 3 * canonstate->m;
                unsigned long *ul_msg = (unsigned long*)malloc(sizeof(unsigned long)*ul_msg_sz);    /* allocate space for message buffer */                    
                MPI_Recv(ul_msg, ul_msg_sz, MPI_UNSIGNED_LONG, recv_status.MPI_SOURCE, MPI_MSG_TAKE_WORK_UL, MPI_COMM_WORLD, &recv_status);

                int pos = 0; /* used to track where in the msg a given object is */
                int ul_pos = 0; /* used to track where in the ul_msg a given object is */

                Node *curr; /* used to hold the current node being deserialized */
                int recv_sz = msg[pos++]; /* first word of the message is the number of nodes sent */
                /* deserialize messages and push to stack */
                for(int i = 0; i < recv_sz; ++i) {
                    curr = node_deserialize(canonstate->n, canonstate->m, &pos, msg, &ul_pos, ul_msg); /* deserialize the node from the message */
                    stack_push(canonstate->stack, curr);    /* push current node to stack, stack now owns it, we don't free it here */
// TODO: DELETE THIS
// fprintf(canonstate->tracefile, "MPI: Process %d: %.10f Recv Node [%d/%d] from %d: ",mpi_state->my_rank, MPI_Wtime(), i, recv_sz, recv_status.MPI_SOURCE); fpath_visualize(canonstate->tracefile, curr->path); fprintf(canonstate->tracefile, "\n");
// TO HERE                          
                }
                if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: Received %d nodes stack %d(%d) in TAKE_WORK [%zu]\n", mpi_state->my_rank, recv_sz, stack_size(canonstate->stack), start_stack_sz, canonstate->numnodes);
                /* free message buffer */
                FREES(msg); 
                FREES(ul_msg);  
                break;
        
            default:
                /* if we get here, then we got a message from process 0, but it wasn't related to our work request */
                mpi_poll_for_messages(mpi_state, canonstate);

// TODO:  Review this logic it might not be needed.            
                /* if we return from poll messages in a work end state, then we should return out of this funciton */
                if (mpi_state->state == MPI_STATE_WORK_END) {
                    return;
                }
                break;
            } /* switch (recv_status.MPI_TAG) */
        } /* if (flag)  meaning if the iprobe saw a message from process 0 */
        else {
            /* we get here if there was no message from process 0 waiting,  poll for messages. */
            mpi_poll_for_messages(mpi_state, canonstate);
                
// TODO:  Review this logic it might not be needed.            
            /* if we return form poll messages in a work end state, then we should return out of this funciton */
            if (mpi_state->state == MPI_STATE_WORK_END) {
                return;
            }
        } /* else for if(flag) */
    }

    /* if state is work received, return */
    if (mpi_state->state == MPI_STATE_WORK_RECEIVED) {
        mpi_state->state = MPI_STATE_WORKING;
        return;
    }

    if (mpi_state->state = MPI_STATE_REJECTED) {
        mpi_state->state = MPI_STATE_WAIT_FOR_WORK_END;
        return;
    }

    /* if we are here, who knows how or why, but throw an error */
    printf("MPI: Process %d: Droped out of ask for work with an invaid state %d\n", mpi_state->my_rank, mpi_state->state);
    exit(mpi_state->state);
}





// void mpi_process_work_request_with_empty_queue(MPIState *mpi_state, CanonState *canonstate, int requesting_rank) {
//     /**
//      * This function if for process 0 only, it is designed to handle a work request from another process, if we don't have any work to give.
//      */
//     if (mpi_state->my_rank != 0){
//         /* if we are here, who knows how or why, but throw an error */
//         printf("MPI: Process %d: is requesting a queue size report, but only process 0 should!\n", mpi_state->my_rank);
//         exit(1);
//     }

//     /* set the process that is requesting work, to the requested state */
//     canonstate->queue_sizes[requesting_rank] = MPI_QUEUE_SIZE_REQUESTED_WORK;

//     if (mpi_state->state != MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES && mpi_state->state != MPI_STATE_WAIT_FOR_WORK_RETURN) {
//         /**
//          * If we are not already in the MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES state, then we need to send a queue size request, and enter this state
//          */
//         int dest[mpi_state->num_processes];
//         if (__DEBUG_MPI_SENDWORK__) _debug_queue_size(mpi_state, canonstate);
//         int pos = 0;
//         for (int i = 1; i < mpi_state->num_processes; ++i){
//             /* we only want to send queue size requests to processes that aren't already waiting for work */
//             if (canonstate->queue_sizes[i] != MPI_QUEUE_SIZE_REQUESTED_WORK) {
//                 dest[pos] = i;  // set up the destination array for the send function
//                 ++pos;
//                 canonstate->queue_sizes[i] = MPI_QUEUE_SIZE_UNDEFINED; // set status for all queue sizes to requested
//             }
//         }
//         canonstate->queue_sizes[0] = stack_size(canonstate->stack);
//         int msg[1];
//         msg[0] = MPI_MSG_REPORT_QUEUE_SIZE;
//         // for (int i = 0; i < mpi_state->num_processes-1; ++i) printf("DEBUG: dest[%d] = %d\n", i, dest[i]);
//         if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: state %d, sending queue size request to %d processes.\n", mpi_state->my_rank, mpi_state->state, pos);
//         if (__DEBUG_MPI_SENDWORK__) {
//             printf("\tprocesses:");
//             for (int i = 0; i < pos; ++i) {
//                 printf("  %d", dest[i]);
//             }
//             printf("\n");
//         }
//         if (__DEBUG_MPI_SENDWORK__) _debug_queue_size(mpi_state, canonstate);

//         _mpi_pseudo_send_to_many(mpi_state, &msg, 1, MPI_INT, dest, pos, MPI_MSG_REPORT_QUEUE_SIZE, MPI_COMM_WORLD, canonstate);
//         mpi_state->state = MPI_STATE_WAIT_FOR_QUEUE_SIZE_RESPONSES;
//     }


//     /* check if all other processes have already requested work */
//     for (int i = 1; i < mpi_state->num_processes; ++i) {
//         if (canonstate->queue_sizes[i] != MPI_QUEUE_SIZE_REQUESTED_WORK) return;  // return if there is still a process running
//     }

//     /* if we are here, then all other processes have requested more work, and hence, are done */
//     mpi_state->state = MPI_STATE_WAIT_FOR_WORK_END;
//     int msg = MPI_MSG_REJECT_NEED_WORK;
//     _mpi_resend_pseudo_broadcast(mpi_state, &msg, 1, MPI_INT, MPI_MSG_REJECT_NEED_WORK, MPI_COMM_WORLD, 0, canonstate);
// }

void mpi_broadcast_work_request_reject(MPIState *mpi_state, CanonState *canonstate) {
    int msg = MPI_MSG_REJECT_NEED_WORK;
    _mpi_resend_pseudo_broadcast(mpi_state, &msg, 1, MPI_INT, MPI_MSG_REJECT_NEED_WORK, MPI_COMM_WORLD, 0, canonstate);
}

void mpi_ask_for_work_return(MPIState *mpi_state, CanonState *canonstate) {
    /**
     * Function should calculate if target_rank has enough work to request a return, and if so, send that request.  If not, return to working state.
     */
    if (mpi_state->my_rank != 0){
        /* if we are here, who knows how or why, but throw an error */
        printf("MPI: Process %d: is trying to request work return, only process 0 is permitted!\n", mpi_state->my_rank);
        exit(1);
    }

    /* figure out which process to request from (really just the one with the largest reported stack size) */
    int target_rank = -1;
    int largest_queue = 0;
    for (int i = 1; i < mpi_state->num_processes; ++i) {
        if (canonstate->queue_sizes[i] > largest_queue) {
            target_rank = i;
            largest_queue = canonstate->queue_sizes[i];
        }
    }

    if (target_rank < 1) {
        printf("MPI: Process %d: is trying to request work return, but no valid target was found!\n", mpi_state->my_rank);
        return;  /* if we ended up not finding a valid targe, no need to continue, just return */
    }

    /**
     * First stab at this:
     * 
     *      if the target's queue size is > max_work_size_to_send, we ask it to return work.  We're expecting it to return 50% of it's queue
     */

    if (canonstate->queue_sizes[target_rank] > canonstate->max_work_size_to_send) {
        

        MPI_Status recv_status;  /* MPI status used for probe and receive functions */

        mpi_state->state = MPI_STATE_WAIT_FOR_WORK_RETURN;   /* set our state machine to asking for work */
        
        int nomsg = MPI_MSG_RETURN_WORK;
        if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: sent work return request to %d\n", mpi_state->my_rank, target_rank);
        _mpi_pseudo_send(mpi_state, &nomsg, 1, MPI_INT, target_rank, MPI_MSG_RETURN_WORK, MPI_COMM_WORLD, canonstate); /* send request for work */
        int flag;

        /* this while loop allows us to keep probing for status to our work request, even if we get another message while waiting */
        while(mpi_state->state == MPI_STATE_WAIT_FOR_WORK_RETURN) {
            /** wait for a response */
            flag = 0;
            MPI_Iprobe(target_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &recv_status);
            if (flag) {
                /* depending on the tag received take action */
                switch (recv_status.MPI_TAG)
                {
                case MPI_MSG_REJECT_RETURN_WORK:
                    /* we got a reject, need to clear the message, and try again */
                    MPI_Recv(&nomsg, 1, MPI_INT, recv_status.MPI_SOURCE, recv_status.MPI_TAG, MPI_COMM_WORLD, &recv_status);
                    mpi_state->state = MPI_STATE_WORKING;
                    if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: Received a reject for work return from  %d\n", mpi_state->my_rank, recv_status.MPI_SOURCE);
                    return;
                    break;
                    
                case MPI_MSG_TAKE_RETURN: 
                    ; /* fake statement to allow for a declaration to come next */
                    /* we got some work */
                    int start_stack_sz = stack_size(canonstate->stack); /* used for debugging, to see how much work we had before we received more */
                    mpi_state->state = MPI_STATE_WORK_RECEIVED; /* set state to work received */

                    int msg_sz;
                    MPI_Get_count(&recv_status, MPI_INT, &msg_sz);
                    int *msg = (int*)malloc(sizeof(int)*msg_sz);    /* allocate space for message buffer */

                    MPI_Recv(msg, msg_sz, MPI_INT, recv_status.MPI_SOURCE, MPI_MSG_TAKE_RETURN, MPI_COMM_WORLD, &recv_status);

                    int ul_msg_sz = msg[0] * 3 * canonstate->m;
                    unsigned long *ul_msg = (unsigned long*)malloc(sizeof(unsigned long)*ul_msg_sz);    /* allocate space for message buffer */                    
                    MPI_Recv(ul_msg, ul_msg_sz, MPI_UNSIGNED_LONG, recv_status.MPI_SOURCE, MPI_MSG_TAKE_RETURN_UL, MPI_COMM_WORLD, &recv_status);

                    int pos = 0; /* used to track where in the msg a given object is */
                    int ul_pos = 0; /* used to track where in the ul_msg a given object is */

                    Node *curr; /* used to hold the current node being deserialized */
                    int recv_sz = msg[pos++]; /* first word of the message is the number of nodes sent */
                    /* deserialize messages and push to stack */
                    for(int i = 0; i < recv_sz; ++i) {
                        curr = node_deserialize(canonstate->n, canonstate->m, &pos, msg, &ul_pos, ul_msg); /* deserialize the node from the message */
                        stack_push(canonstate->stack, curr);    /* push current node to stack, stack now owns it, we don't free it here */
                        
                    }
                    if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: Received %d nodes from %d stack %d(%d) in TAKE_RETURN\n", mpi_state->my_rank, recv_sz, recv_status.MPI_SOURCE, stack_size(canonstate->stack), start_stack_sz);

                    if (__DEBUG_MPI_PEROIDIC_STATUS__ && mpi_state->my_rank == 0){
                        /**
                         * This is a periodic update for debugging purposes
                         */
                        float wtime = MPI_Wtime();
                        printf("%011f MPI: Process %d: PSU(TR): [0] %d, -1", wtime, mpi_state->my_rank, stack_size(canonstate->stack));
                        for(int i = 1; i < mpi_state->num_processes; ++i) printf("   [%d] %d, %d", i, canonstate->queue_sizes[i], canonstate->process_waiting_for_work[i]);
                        printf("\n");
                        mpi_state->mpi_last_peroidic_update = wtime;
                    }  


                    /* free message buffer */
                    FREES(msg); 
                    FREES(ul_msg);  
                    break;
            
                default:
                    /* if we get here, then we got a message from process target_rank, but it wasn't related to our work return request */
                    mpi_poll_for_messages(mpi_state, canonstate);

    // TODO:  Review this logic it might not be needed.            
                    /* if we return from poll messages in a work end state, then we should return out of this funciton */
                    if (mpi_state->state == MPI_STATE_WORK_END) {
                        return;
                    }
                    break;
                } /* switch (recv_status.MPI_TAG) */
            } /* if (flag)  meaning if the iprobe saw a message from process target_rank */
            else {
                /* we get here if there was no message from process target_rank waiting,  poll for messages. */
                mpi_poll_for_messages(mpi_state, canonstate);
                    
    // TODO:  Review this logic it might not be needed.            
                /* if we return form poll messages in a work end state, then we should return out of this funciton */
                if (mpi_state->state == MPI_STATE_WORK_END) {
                    return;
                }
            } /* else for if(flag) */
        }  //while(mpi_state->state == MPI_STATE_WAIT_FOR_WORK_RETURN)

        /* if state is work received, return */
        if (mpi_state->state == MPI_STATE_WORK_RECEIVED) {
//TODO:  HERE   =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
            mpi_state->state = MPI_STATE_WORKING;
            for (int i = 1; i < mpi_state->num_processes; ++i) {
                /**
                 * the only reason we'd be here is if we ran out of work,requested work return, and got some work returned, therefore
                 *      we should send work to any processes that have requsted it.
                 */

                if (canonstate->process_waiting_for_work[i] == TRUE) {
                    if (stack_size(canonstate->stack) < 1) {
                        /* We are here because a process (i) has request work, but we are once again out, simply print and move on */
                        // if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: has had work returned, but no longer has any to distribute to %d\n", mpi_state->my_rank, i);
                        // continue;  // Trying this, so the loop can finish, and print out the processes that we can't service.  in prod, we can return here.
                        break;  // Just breaking so we can send debug message at end.  in prod, we can return here.
// TODO: Can change this continue to a return.                        
                        // return;   // if our queue is empty again, exit loop so we can return to normal
                    }                    
                    /**
                     * if here, the current stack/queue is not empty, and there is is a process, i, that has previously requested work and is waiting,
                     *      send it some work!
                     */
                    // if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: has had work returned, and is sending work to %d\n", mpi_state->my_rank, i);
                    canonstate->process_waiting_for_work[i] = FALSE;
                    _mpi_send_work_to(mpi_state, canonstate, i);
                }
            }


            if (__DEBUG_MPI_PEROIDIC_STATUS__ && mpi_state->my_rank == 0){
                /**
                 * This is a periodic update for debugging purposes
                 */
                float wtime = MPI_Wtime();
                printf("%011f MPI: Process %d: PSU(ATR): [0] %d, -1", wtime, mpi_state->my_rank, stack_size(canonstate->stack));
                for(int i = 1; i < mpi_state->num_processes; ++i) printf("   [%d] %d, %d", i, canonstate->queue_sizes[i], canonstate->process_waiting_for_work[i]);
                printf("\n");
                mpi_state->mpi_last_peroidic_update = wtime;
            } 



            return;
        }


        /* if we are here, who knows how or why, but throw an error */
        printf("MPI: Process %d: Droped out of ask for work return with an invaid state %d\n", mpi_state->my_rank, mpi_state->state);
        exit(mpi_state->state);        
    }
}


void mpi_return_work(MPIState *mpi_state, CanonState *canonstate) {
    if (stack_size(canonstate->stack) > 1) { // we want to send at least one back!
        int start_stack_sz = stack_size(canonstate->stack);
        int WORKSPLIT = 2;

        /* we have enough work to send */
        int send_sz = stack_size(canonstate->stack) / WORKSPLIT;  
        // if (send_sz > canonstate->max_work_size_to_send) send_sz = canonstate->max_work_size_to_send;  /* limit amount of work to send in one chunk */
        Node *curr;
        int buff_sz = 1;  /* start with 1 for the record count */
        int ul_buff_sz = 0; /* used for unsigned long serialization */

        Stack *temp_stack;
        DYNALLOCSTACK(temp_stack, "temp_stack in mpi_return_work");

        /**
         * Walk through the stack, from the bottom up, popping out every canonstate->work_split value.
         * Push each of these onto temp_stack, using enqueue, to preserve FIFO order (bottom of original stack items go on top of temp_stack, so they are sent first, and end up on bottom of destination stack
         * While doing this loop, we can also calculate the buffer size needed for sending these nodes.
         */
        for (int i = send_sz - 1; i >= 0; --i) {
            curr = stack_pop_at(canonstate->stack, i * WORKSPLIT); /* pop every canonstate->work_split node from the stack, in revers order, up to send_sz */
            stack_enqueue(temp_stack, curr);  /* using enqueue here, so temp_stack works like a FIFO queue, this preserves the order */
            node_get_serialization_size(curr, &buff_sz, &ul_buff_sz, canonstate->n, canonstate->m); /* get the size of the node at stack index i */
        }

        int *msg = (int*)malloc(sizeof(int)*buff_sz);   /* allocate message buffer */
        unsigned long *ul_msg = (unsigned long*)malloc(sizeof(unsigned long)*ul_buff_sz); /* allocate unsigned long message buffer */

        
        /* fill the message buffer */
        int idx = 0;                                          /* variable to be used as the message index */
        int ul_idx = 0;                                       /* variable to be used as the unsigned long message index */
        msg[idx++] = send_sz;                                 /* first word of the message is the number of nodes sent */
        
        /* remove nodes from temp_stack and add to messages */
        for (int i = 0; i < send_sz; ++i) {
            curr = stack_pop(temp_stack);                        
            node_serialize(curr, canonstate->n, canonstate->m, &idx, msg, &ul_idx, ul_msg);  /* serialize the node into the message buffer */
            FREENODE(curr);                         /* free the node, as we don't need it anymore */
        }

        // set up arrays to pass to send_multiple
        void *buffers[2] = {msg, ul_msg};
        int sizes[2] = {buff_sz, ul_buff_sz};
        MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};
        int tags[2] = {MPI_MSG_TAKE_RETURN, MPI_MSG_TAKE_RETURN_UL};
        
        if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: Returning %d nodes of %d (%d) on stack to 0 in RETURN_WORK [%zu]\n", mpi_state->my_rank, send_sz,start_stack_sz, stack_size(canonstate->stack), canonstate->numnodes);
        _mpi_pseudo_send_multiple(mpi_state, 2, buffers, sizes, types, 0, tags, MPI_COMM_WORLD, canonstate);

        /* free message buffer */
        FREES(msg);
        FREES(ul_msg);
    } else {
        int msg = MPI_MSG_REJECT_RETURN_WORK;
        MPI_Request request;
        MPI_Isend(&msg, 1, MPI_INT, 0, MPI_MSG_REJECT_RETURN_WORK, MPI_COMM_WORLD, &request);
        if (__DEBUG_MPI_WORK_RETURN__) printf("MPI: Process %d: Rejecting work request, because stack is empty\n", mpi_state->my_rank);
    }
}

// TODO:  Proboably need to change these so only process 0 does broadcasts

void mpi_send_new_best_cl(MPIState *mpi_state, Path *canonpath, CodePath *canoncodes, int *canonlab, int canonlevel, int samerows, int n, CanonState *canonstate) {
    int msg_sz = canonpath->sz + canoncodes->sz + n + 4;  /* Size of message to send */

    int *msg = (int*)malloc(sizeof(int)*msg_sz);   /* allocate message buffer */

    int pos = 0;  /* variable to be used as the message index */
    
    /** Add the path to the message */
    msg[pos++] = canonpath->sz;                      
    /* loop through the path and put the path words into the message */
    for (int i = 0; i < canonpath->sz; ++i) {
        msg[pos++] = canonpath->data[i];             
    }
    /** */

    /** Add the code path to the message */
    msg[pos++] = (int)(canoncodes->sz);
    /* loop through the code path and put the code path words into the message */
    for (int i = 0; i < canoncodes->sz; ++i) {
        msg[pos++] = canoncodes->data[i];
    }
    /** */

    /** Add the parition (pi) to the message */
    /* loop through the partition and put the lab words into the message */
    for (int i = 0; i < n; ++i) {
        msg[pos++] = canonlab[i];
    }
    /** */
    
    /** Add the canonlevel and samerows to the message */
    msg[pos++] = canonlevel;
    msg[pos++] = samerows;
    /** */

    if (__DEBUG_MPI_NEWCL__) {printf("MPI: Process %d: Start Broadcast New Best CL in %d words  state %d\n", mpi_state->my_rank, msg_sz, mpi_state->state); }

    // _mpi_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_NEW_CL, MPI_COMM_WORLD, canonstate);
    if (mpi_state->my_rank == 0) {
        /* if process zero is sending the message, just do the pseudobroadcast */
        _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_NEW_CL, MPI_COMM_WORLD, 0, canonstate); 
    } else {
        /* Otherwise, send only to zero */
        _mpi_pseudo_send(mpi_state, msg, msg_sz, MPI_INT, 0, MPI_MSG_NEW_CL, MPI_COMM_WORLD, canonstate); 
    }

    if (__DEBUG_MPI_NEWCL__) {printf("MPI: Process %d: Finished Broadcast New Best CL\n", mpi_state->my_rank);}

    FREES(msg);
}

void mpi_send_new_automorphism(MPIState *mpi_state, set *fix, set *mcr, int m, CanonState *canonstate) {

    int msg_sz = 2 * m;  /* Set message size */

    unsigned long *msg = (unsigned long*)malloc(sizeof(unsigned long)*msg_sz);   /* allocate message buffer */

    int pos = 0;  /* variable to be used as the message index */
    
    /* Add the fix partition */
    for (int i = 0; i < m; ++i) {
        msg[pos++] = fix[i];
    }
    /* Add the mcr partition */
    for (int i = 0; i < m; ++i) {
        msg[pos++] = mcr[i];
    }

    if (__DEBUG_MPI_NEWAUTO__) {
        printf("MPI: Process %d: Broadcast New Automorphism in %d words  workspace avail: %ld\n", mpi_state->my_rank, msg_sz, (canonstate->worktop - *(canonstate->fmptr)) / (2*canonstate->m));
    }

    // _mpi_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_UNSIGNED_LONG, MPI_MSG_NEW_AUTO, MPI_COMM_WORLD, canonstate);
    if (mpi_state->my_rank == 0) {
        /* if process zero is sending the message, just do the pseudobroadcast */

        _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_UNSIGNED_LONG, MPI_MSG_NEW_AUTO, MPI_COMM_WORLD, 0, canonstate); 
    } else {    
        /* Otherwise, send only to zero */
        _mpi_pseudo_send(mpi_state, msg, msg_sz, MPI_UNSIGNED_LONG, 0, MPI_MSG_NEW_AUTO, MPI_COMM_WORLD, canonstate); 
    }
    FREES(msg); /* free the message buffer */
}


void mpi_send_active_prune_message(MPIState *mpi_state, Path *path, int to_level, CanonState *canonstate) {
    int msg_sz = path->sz + 2; /* message size is path size + 2 for the level and path->sz */
    int msg[msg_sz]; /* message size is path size + 2 for the level and samerows */
    int idx = 0; /* index for the message */

    if (path->sz < 1) {
        printf("MPI: Process %d: request to send prune with path len zero!\n", mpi_state->my_rank);
        path_visualize(path);
        exit(MPI_MSG_ACTIVE_PRUNE);
    }
    msg[idx++] = path->sz; /* first word is the size of the path */

    for (int i = 0; i < path->sz; ++i) {
        msg[idx++] = path->data[i]; /* copy the path data */
    }
    msg[idx++] = to_level; /* add the level */

    if (__DEBUG_MPI_PRUNE__) {printf("MPI: Process %d: sending active prune to level %d\n", mpi_state->my_rank, to_level);} 

    // _mpi_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_ACTIVE_PRUNE, MPI_COMM_WORLD, canonstate);
    if (mpi_state->my_rank == 0) {
        /* if process zero is sending the message, just do the pseudobroadcast */
        _mpi_resend_pseudo_broadcast(mpi_state, msg, msg_sz, MPI_INT, MPI_MSG_ACTIVE_PRUNE, MPI_COMM_WORLD, 0, canonstate); 
    } else {
        /* Otherwise, send only to zero */    
        _mpi_pseudo_send(mpi_state, msg, msg_sz, MPI_INT, 0, MPI_MSG_ACTIVE_PRUNE, MPI_COMM_WORLD, canonstate);
    }
}


void _mpi_send_work_to(MPIState *mpi_state, CanonState *canonstate, int target_process){
    int start_stack_sz = stack_size(canonstate->stack);

    /* assuming, if we are here, then we have enough work to send */
    int send_sz = (stack_size(canonstate->stack) + (canonstate->work_split - 1)) / canonstate->work_split;  /* send amount of work, rounded up (hence the + (canonstate->work_split -1)) */
    if (send_sz > canonstate->max_work_size_to_send) send_sz = canonstate->max_work_size_to_send;  /* limit amount of work to send in one chunk */
    Node *curr;
    int buff_sz = 1;  /* start with 1 for the record count */
    int ul_buff_sz = 0; /* used for unsigned long serialization */

    Stack *temp_stack;
    DYNALLOCSTACK(temp_stack, "temp_stack in MPI_MSG_NEED_WORK");

    /**
     * Walk through the stack, from the bottom up, popping out every canonstate->work_split value.
     * Push each of these onto temp_stack, using enqueue, to preserve FIFO order (bottom of original stack items go on top of temp_stack, so they are sent first, and end up on bottom of destination stack
     * While doing this loop, we can also calculate the buffer size needed for sending these nodes.
     */
    for (int i = send_sz - 1; i >= 0; --i) {
        curr = stack_pop_at(canonstate->stack, i * canonstate->work_split); /* pop every canonstate->work_split node from the stack, in revers order, up to send_sz */
        stack_enqueue(temp_stack, curr);  /* using enqueue here, so temp_stack works like a FIFO queue, this preserves the order */
        node_get_serialization_size(curr, &buff_sz, &ul_buff_sz, canonstate->n, canonstate->m); /* get the size of the node at stack index i */
    }

    int *msg = (int*)malloc(sizeof(int)*buff_sz);   /* allocate message buffer */
    unsigned long *ul_msg = (unsigned long*)malloc(sizeof(unsigned long)*ul_buff_sz); /* allocate unsigned long message buffer */

    
    /* fill the message buffer */
    int idx = 0;                                          /* variable to be used as the message index */
    int ul_idx = 0;                                       /* variable to be used as the unsigned long message index */
    msg[idx++] = send_sz;                                 /* first word of the message is the number of nodes sent */
    
    /* remove nodes from temp_stack and add to messages */
    for (int i = 0; i < send_sz; ++i) {
        curr = stack_pop(temp_stack);                        
        node_serialize(curr, canonstate->n, canonstate->m, &idx, msg, &ul_idx, ul_msg);  /* serialize the node into the message buffer */
        FREENODE(curr);                         /* free the node, as we don't need it anymore */
    }

    // set up arrays to pass to send_multiple
    void *buffers[2] = {msg, ul_msg};
    int sizes[2] = {buff_sz, ul_buff_sz};
    MPI_Datatype types[2] = {MPI_INT, MPI_UNSIGNED_LONG};
    int tags[2] = {MPI_MSG_TAKE_WORK, MPI_MSG_TAKE_WORK_UL};
    
    if (__DEBUG_MPI_SENDWORK__) printf("MPI: Process %d: Sending %d nodes of %d (%d) on stack to %d in NEED_WORK [%zu]\n", mpi_state->my_rank, send_sz,start_stack_sz, stack_size(canonstate->stack), target_process, canonstate->numnodes);
    canonstate->queue_sizes[target_process] = send_sz; /* we are sending send_sz nodes to target_process, so we can assume that process was out of work and will now have only the work we are sending it */
    _mpi_pseudo_send_multiple(mpi_state, 2, buffers, sizes, types, target_process, tags, MPI_COMM_WORLD, canonstate);

    /* free message buffer */
    FREES(msg);
    FREES(ul_msg);

}

void _mpi_resend_pseudo_broadcast(MPIState *mpi_state, const void *msg, int msg_sz, MPI_Datatype datatype, int tag, MPI_Comm comm, int origin, CanonState *canonstate) {
    if (mpi_state->my_rank != 0) {
        printf ("MPI: Process %d: Called resend _mpi_resend_pseudo_broadcast, but not process 0\n", mpi_state->my_rank);
        exit(-1);
    }


    int dest[mpi_state->num_processes];
    int pos = 1, dest_sz = 0;
    while (pos < mpi_state->num_processes) {
        // if (pos != origin && canonstate->processes_running[pos]) dest[dest_sz++] = pos;
        if (pos != origin) dest[dest_sz++] = pos;
        ++pos;
    }

    if (dest_sz > 0) {
        _mpi_pseudo_send_to_many(mpi_state, msg, msg_sz, datatype, dest, dest_sz, tag, comm, canonstate);
    }
}

void _mpi_pseudo_send(MPIState *mpi_state, const void *buf, int buff_sz, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, CanonState *canonstate){
    MPI_Request request;  /* used to track request and ensure it completes */
    
    /* send message */
    MPI_Isend(buf, buff_sz, datatype, dest, tag, comm, &request);

    /* wait for communication to finish, while waiting, keep processing in coming message*/
    /* this is needed to prevent deadlocks */
    int flag = FALSE;
    while (!flag) {
        MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
        if (!flag) {
            mpi_poll_for_messages(mpi_state, canonstate);
        }
    }
}

void _mpi_pseudo_send_multiple(MPIState *mpi_state, int count, void **buf, int *buff_sz, MPI_Datatype *datatype, int dest, int *tags, MPI_Comm comm, CanonState *canonstate){
    MPI_Request requests[count]; /* used for tracking the requests */
    
    for (int i = 0; i < count; ++i){
        /* send message(s) */
        MPI_Isend(buf[i], buff_sz[i], datatype[i], dest, tags[i], comm, &requests[i]);
    }
    /* wait for communication to finish, while waiting, keep processing in coming message*/
    /* this is needed to prevent deadlocks */
    int flag = FALSE;
    while (!flag) {
        MPI_Testall(count, requests, &flag, MPI_STATUS_IGNORE);
        if (!flag) {
            mpi_poll_for_messages(mpi_state, canonstate);
        }
    }
}


void _mpi_pseudo_send_to_many(MPIState *mpi_state, const void *buf, int buff_sz, MPI_Datatype datatype, int *dest, int dest_sz, int tag, MPI_Comm comm, CanonState *canonstate) {
    
    MPI_Request requests[dest_sz]; /* used for tracking the requests */
    /* initiate the non-blocking sends */
    for (int i = 0; i < dest_sz; ++i){
        if (dest[i] == mpi_state->my_rank) {
            requests[i] = MPI_REQUEST_NULL;
        } else {
            MPI_Isend(buf, buff_sz, datatype, dest[i], tag, comm, &requests[i]);
        }
    }

    /* wait for communication to finish, while waiting, keep processing incoming messages */
    /* this is needed to prevent deadlocks */
    int flag = FALSE;
    while (!flag) {
        MPI_Testall(dest_sz, requests, &flag, MPI_STATUSES_IGNORE);
        if (!flag) {
            mpi_poll_for_messages(mpi_state, canonstate);
        }
    }
}
