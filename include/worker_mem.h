/*
   worker_mem.h

   Created by Petter Taule on 06.08.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include "tables.h"

#ifndef WORKER_MEM_H
#define WORKER_MEM_H

typedef struct {
    table_ptrs_t* data;                      /* Array of table (pointers), one for each process */
    size_t size;                             /* Size of dyn allocated array                     */
    const short int (*sum_table)[N_CONFIGS]; /* Pointer to sum table (workers use the same)     */
    const double* eta;                       /* Pointer to eta array (workers use the same)     */
} worker_mem_t;


void init_worker_mem(
        worker_mem_t* wm,
        size_t init_size,
        const short int sum_table[][N_CONFIGS],
        const double* eta
        );

void resize_worker_mem(worker_mem_t* wm, size_t new_size);
void worker_mem_gc(worker_mem_t* wm);

#endif /* ifndef WORKER_MEM_H */
