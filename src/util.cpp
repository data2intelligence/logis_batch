/*
 * util.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: peng
 */

#include "util.h"

void common_names(
		const string_array *X, const string_array *Y, const string_array *B,
		vector<string> &result)
{
	size_t i;
	string s;
	set<string> sX, sY, sM, sB;

	for(i=0;i<X->n;i++) sX.insert(X->array[i]);
	for(i=0;i<Y->n;i++) sY.insert(Y->array[i]);
	for(i=0;i<B->n;i++) sB.insert(B->array[i]);

	intersection(sX, sY, sM);
	intersection(sM, sB, sY);

	// set gene order by Y, assumed to be the largest matrix
	result.clear();

	for(i=0;i<Y->n;i++)
	{
		s = Y->array[i];
		if(sY.find(s) != sY.end()) result.push_back(s);
	}
}

