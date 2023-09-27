#include "string_array.h"
#include <stdlib.h>
#include <stdio.h>

string_array *init_string_array(const char buffer[], const char sep)
{
	char *ptr, **array;
	size_t i, j, n, length = strlen(buffer);

	string_array *result = (string_array*)malloc(sizeof(string_array));

	// with extra \0 at very end
	ptr = result->buffer = (char*) malloc((length+1)*sizeof(char));

	// copy the parse buffer to string buffer
	memcpy(ptr, buffer, (length+1)*sizeof(char));
	result->length = length;

	// count how many strings are there
	for(n=1,i=0; i<length; i++)
	{
		if(ptr[i] == sep)
		{
			ptr[i] = '\0';
			n++;
		}
	}

	result->n = n;

	// hook up the array to begin of each string
	array = result->array = (char **) malloc(n*sizeof(char*));
	array[0] = ptr;

	for(j=1,i=0 ; i<length ; i++,ptr++)
	{
		if(*ptr == '\0' && j < n)
		{
			array[j] = ptr + 1;
			j++;
		}
	}

	return result;
}

void free_string_array(string_array *arr)
{
	free(arr->buffer);
	free(arr->array);
	free(arr);
}


void write_string_array(const string_array *array,
		char buffer[], const size_t length, const char sep)
{
	size_t i, arr_len = array->length;

	if(arr_len + 1 > length)
	{
		fprintf(stderr, "Insufficient space for string array output\n");
		exit(1);
	}

	memcpy(buffer, array->buffer, (arr_len + 1) * sizeof(char));

	for(i=0;i<arr_len;i++)
	{
		if(buffer[i] == '\0') buffer[i] = sep;
	}
}


void check_string_buffer_end(const char buffer[], const size_t length)
{
	size_t i;
	for(i=0;i<length;i++)
	{
		if(buffer[i] == '\0') break;
	}

	if(i >= length-1)
	{
		fprintf(stderr, "buffer size %lu is not enough to parse string", length);
		exit(1);
	}
}
