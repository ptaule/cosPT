/*
   utilities.c

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <math.h>

#include "utilities.h"

void number_to_base3(
        short int number, // in, number in base 10
        char coefficients[], // out, array of base 3 digits
        size_t size          // in, size of coefficients
        )
{
    for (size_t i = 0; i < size; ++i) {
        coefficients[i] = number % 3;
        number /= 3;
    }
}

void base3_to_number(
        const char coefficients[], // in, array of base 3 digits
        size_t size,               // in, size of coefficients
        short int* number          // out, number in base 10
        )
{
    *number = 0;
    for (size_t i = 0; i < size; ++i) {
        *number += coefficients[i] * pow(3,i);
    }
}

