/*
 * util.h
 *
 *  Created on: Apr 14, 2014
 *      Author: peng
 */

#ifndef UTIL_H_
#define UTIL_H_

extern "C" {
#include "string_array.h"
}


#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

#define MAXLINE 8000000 // ~20000 genes x ~400 char length

// Utils for parameter reading, parsing and checking
#define DISPLAY_BOOL(x) ((x)?("1 (yes)"):("0 (no)"))

template <class T>
T load_positive_value(const string s, const string fieldname, const T minvalue, const T maxvalue)
{
	istringstream iss(s);
	T n;

	if(s[0] == '-') {
		cerr << "Please input a positive value for field \"" << fieldname << "\"." << endl;
		cerr << s << " is your input." << endl;
		exit(1);
	}

	iss >> n;

	if(iss.fail())
	{
		cerr << "Please input a valid value for field \"" << fieldname << "\"." << endl;
		cerr << s << " is your input." << endl;
		exit(1);
	}

	if(n < minvalue) {
		cerr << "The minimum value for field \"" << fieldname << "\" is " << minvalue << '.' << '\n';
		cerr << s << " is your input." << endl;
		exit(1);
	}

	if(n > maxvalue) {
		cerr << "The maximum value for field \"" << fieldname << "\" is " << maxvalue << '.' << '\n';
		cerr << s << " is your input." << endl;
		exit(1);
	}

	return n;
}



void common_names(
		const string_array *X, const string_array *Y, const string_array *B,
		vector<string> &result);


// read matrix with column names and row names.
// fname: input file name
// data: space of allocated data space.
// Rformat: true if there is no start column name for row names column
// offset: offset all read values
// left_pad, right_pad: number of extra elements to pad in the begin and end

template <class T>
void read_matrix(
		const string &fname, T *&data, string_array *&rownames, string_array *&colnames,
		const bool Rformat = false,
		const T offset=0,
		const vector<string> *left_pad=NULL,
		const vector<string> *right_pad=NULL,
		const T padv = 0)
{
	int streampos;
	size_t i, j, inx, n, p, total;
	T value;

	char buffer[MAXLINE], *buffer_ptr = buffer;

	vector<string>::const_iterator iter;

	ifstream fin(fname.c_str());

	if(fin.fail())
	{
		cerr << "Cannot open " << fname << endl;
		exit(1);
	}

	// read in column names
	inx = 0;		// counter of current read position

	if(left_pad!=NULL)
	{
		for(iter=left_pad->begin(); iter!=left_pad->end(); iter++)
		{
			total = iter->size();
			memcpy(buffer_ptr, iter->data(), total * sizeof(char));
			buffer_ptr[total] = '\t';
			buffer_ptr += (total + 1);
			inx += (total + 1);
		}
	}

	if(inx >= MAXLINE-1)
	{
		cerr << "Exceed maximum parse buffer limit. Please change MAXLINE and compile again." << endl;
		exit(1);
	}

	fin.getline(buffer_ptr, MAXLINE - inx, '\n');
	buffer_ptr += strlen(buffer_ptr);
	inx += strlen(buffer_ptr);

	if(fin.fail() && !fin.eof())
	{
		cerr << "Exceed maximum parse buffer limit. Please change MAXLINE and compile again." << endl;
		exit(1);
	}

	if(right_pad!=NULL)
	{
		*buffer_ptr = '\t';
		buffer_ptr++;

		for(iter=right_pad->begin(); iter!=right_pad->end(); iter++)
		{
			total = iter->size();
			memcpy(buffer_ptr, iter->data(), total * sizeof(char));
			buffer_ptr[total] = '\t';
			buffer_ptr += (total + 1);
			inx += (total + 1);
		}

		if(inx >= MAXLINE)
		{
			cerr << "Exceed maximum parse buffer limit. Please change MAXLINE and compile again." << endl;
			exit(1);
		}

		*(buffer_ptr-1) = '\0';
	}

	check_string_buffer_end(buffer, MAXLINE);
	colnames = init_string_array(buffer, '\t');

	// read in row names
	if(!Rformat)
	{
		fin.getline(buffer, MAXLINE, '\n');

		if(fin.fail() && !fin.eof())
		{
			cerr << "Exceed maximum parse buffer limit. Please change MAXLINE and compile again." << endl;
			exit(1);
		}

	}else{
		streampos = fin.tellg();
		buffer_ptr = buffer;

		// count current read position
		inx = 0;

		while(fin.getline(buffer_ptr, MAXLINE - inx, '\n'))
		{
			for(; *buffer_ptr != '\t'; buffer_ptr++, inx++);

			buffer_ptr++;	// for next read in
			inx++;
		}

		if(fin.fail() && !fin.eof())
		{
			cerr << "Exceed maximum parse buffer limit. Please change MAXLINE and compile again." << endl;
			exit(1);
		}

		*(buffer_ptr-1) = '\0';

		// rewind to previous position
		fin.clear();
		fin.seekg(streampos);
	}

	check_string_buffer_end(buffer, MAXLINE);
	rownames = init_string_array(buffer, '\t');

	p = colnames->n;
	n = rownames->n;

	data = (T*)malloc(n*p*sizeof(T));

	// get the numbers to read from data matrix
	total = p;
	if(left_pad!=NULL) total -= left_pad->size();
	if(right_pad!=NULL) total -= right_pad->size();

	for (inx=i=0; i<n; i++)
	{
		// jump row name for R format
		if(Rformat){fin.getline(buffer, MAXLINE, '\t');}

		if(left_pad != NULL)
		{
			for(j=0; j<left_pad->size(); j++,inx++)
			{
				data[inx] = padv;
			}
		}

		for(j=0; j<total; j++,inx++)
		{
			fin >> value;

			if(fin.fail())
			{
				cerr << "Failure on reading numbers" << endl;
				exit(1);
			}

			data[inx] = value - offset;
		}

		if(right_pad != NULL)
		{
			for(j=0; j<right_pad->size(); j++,inx++)
			{
				data[inx] = padv;
			}
		}
	}

	fin.close();
}


template <class T>
void write_matrix(const string &fname,
		const T *values, const string_array *rownames, const string_array *colnames,
		const bool Rformat=false)
{
	size_t i, j, inx, n = rownames->n, p = colnames->n;

	ofstream fout(fname.c_str());

	for(i=0;i<p;i++){fout << colnames->array[i] << (i==(p-1)?'\n':'\t');}

	if(!Rformat)
	{
		for(i=0;i<n;i++){fout << rownames->array[i] << (i==(n-1)?'\n':'\t');}
	}

	//fout.precision(20); << std::scientific

	for(inx=i=0;i<n;i++)
	{
		if(Rformat)
		{
			fout << rownames->array[i] << '\t';
		}

		for(j=0;j<p;j++,inx++)
		{
			fout << values[inx] << (j==(p-1)?'\n':'\t');
		}
	}
	fout.close();
}


template <class T>
void intersection(const set<T> &v1, const set<T> &v2, set<T> &result)
{
	result.clear();

	typename std::set<T>::const_iterator
		first1 = v1.begin(), last1 = v1.end(),
		first2 = v2.begin(), last2 = v2.end();

	while (first1 != last1 && first2 != last2)
	{
		if (*first1 < *first2)
			first1++;
		else if (*first2 < *first1)
			first2++;
		else {
			result.insert(*first1);
			first1++;
			first2++;
		}
	}
}

template <class T>
void compact_matrix(
		string_array *names, T data[], const size_t stride,
		const vector<string> &name_order, const bool sequential)
{
	char *sptr_in, *sptr_out;
	size_t i, pos;
	map<string, size_t> inxmap;
	map<string, size_t>::iterator miter;

	vector<string>::const_iterator viter;

	T *buffer_in = NULL, *buffer_out = NULL;

	if(!sequential)
	{
		buffer_in = new T[stride];
		buffer_out = new T[stride];
	}

	for(viter=name_order.begin(), i=0; viter!=name_order.end(); viter++, i++)
	{
		inxmap[*viter] = i;
	}

	for(i= 0 ; i < names->n ; i++)
	{
		// not included
		miter = inxmap.find(names->array[i]);
		if(miter == inxmap.end()) continue;

		pos = miter->second;

		// no need to change
		if(i == pos) continue;

		if(sequential)
		{
			memcpy(data + pos*stride, data + i*stride, sizeof(T)*stride);
			names->array[pos] = names->array[i];

		}else{
			// cyclic permutation

			memcpy(buffer_in, data + i*stride, sizeof(T)*stride);
			sptr_in = names->array[i];
			inxmap.erase(miter);		// break this in the chain

			while(pos !=i)
			{
				// move out the pos target to buffer
				memcpy(buffer_out, data + pos*stride, sizeof(T)*stride);
				sptr_out = names->array[pos];

				// move current to pos target
				memcpy(data + pos*stride, buffer_in, sizeof(T)*stride);
				names->array[pos] = sptr_in;

				// next place
				miter = inxmap.find(sptr_out);
				if(miter == inxmap.end()) break;

				swap(buffer_in, buffer_out);
				swap(sptr_in, sptr_out);

				pos = miter->second;
			}

			if(pos == i)
			{	// close cycle
				memcpy(data + i*stride, buffer_in, sizeof(T)*stride);
				names->array[i] = sptr_in;
			}
		}
	}

	names->n = name_order.size();

	if(!sequential)
	{
		delete[] buffer_in;
		delete[] buffer_out;
	}
}


#endif /* UTIL_H_ */
