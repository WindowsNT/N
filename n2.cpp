// n2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#define N_MT
#define N_SMALLVECTOR
#include "n.h"


/*
void* operator new(size_t size)
{
	void* p = malloc(size);
	return (void*)p;
}

void operator delete(void* p)
{
	if (!p)
		return; // OK
	free(p);
}

*/

int main()
{
	N<> n1 = 123LL;
	N<> n2 = "109000654700000673000000";
	auto n3 = n2;


	// Test
	std::cout << n3.s() << std::endl;

	// Num Digits
	std::cout << "Digits of n3: " << n3.NumDigits() << std::endl;

	// bits
	n3 <<= 4;
	std::cout << n3.s() << std::endl;

	// Add
	n3 += n2;
	std::cout << n3.s() << std::endl;

	// mul
	n3 = n2 * n1;
	std::cout << n3.s() << std::endl;

	// parts
	n3 = n2.upperpart(4);
	std::cout << n3.s() << std::endl;
	
	// Clear
	n3.clear();
	std::cout << n3.s() << std::endl;

	// Rand
	n3 = N<>::rand();
	std::cout << n3.s() << std::endl;


}

