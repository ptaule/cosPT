/*
   worker_mem.c

   Created by Petter Taule on 06.08.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include "../include/worker_mem.h"


void init_worker_mem(
        worker_mem_t* wm,
        size_t init_size,
        const short int sum_table[][N_CONFIGS],
        const double* eta
        )
{
    wm->data = (table_ptrs_t*)malloc(init_size * sizeof(table_ptrs_t));
    wm->size = init_size;
    wm->sum_table = sum_table;
    wm->eta = eta;
}


void resize_worker_mem(worker_mem_t* wm, size_t new_size) {
    // Only resize if new_size is larger
    if (new_size <= wm->size) return;
    table_ptrs_t* new = (table_ptrs_t*)realloc(wm->data, new_size * sizeof(table_ptrs_t));
    if (!new) {
        warning("Realloc was unable to allocate memory.");
        return;
    }
    wm->data = new;
    wm->size = new_size;
    printf("New size = %ld\n", new_size);
}

void worker_mem_gc(worker_mem_t* wm) {
    free(wm->data);
    wm->data = NULL;
    wm->size = 0;
    wm->sum_table = NULL;
    wm->eta = NULL;
}
