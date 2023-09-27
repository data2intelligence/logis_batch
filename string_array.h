/*
 * string_array.h
 *
 *  Created on: Apr 14, 2014
 *      Author: Peng Jiang
 */

#ifndef STRING_ARRAY_H_
#define STRING_ARRAY_H_

#include <string.h>

// memory compact string array
typedef struct string_array
{
	char **array, *buffer;
	size_t n, length;

}string_array;

// initialize string array from a string line with separator
string_array *init_string_array(const char buffer[], const char sep);

// free space allocated for string array
void free_string_array(string_array *arr);


// output string array to buffer. length is the size of buffer
void write_string_array(const string_array *array,
		char buffer[], const size_t length, const char sep);


// some tools for sanity check

// make sure current buffer has '\0' to end
void check_string_buffer_end(const char buffer[], const size_t length);

#endif /* STRING_ARRAY_H_ */
