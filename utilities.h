/*
   utilities.h

   Created by Petter Taule on 24.01.2019
   Copyright (c) 2019 Petter Taule. All rights reserved.
*/

#include <stddef.h>

// Constants:
const double PI = 3.14159265359;

// Macros:
//


void number_to_base3(short int number, char coefficients[], size_t size);
void base3_to_number(const char coefficients[], size_t size, short int* number);
