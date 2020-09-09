#ifndef __DNA_H__
#define __DNA_H__

#include "stdio.h"
#include <cstring>
#include<iostream>
#include <stdexcept>

class InvalidNucleotide:public std::invalid_argument
{
public:
	InvalidNucleotide();
	/*virtual*/ const char* what() const throw();
};
inline InvalidNucleotide::InvalidNucleotide():std::invalid_argument(""){}
inline const char* InvalidNucleotide::what() const throw()
{
	return "invalid dna";
}
#endif //__DNA_H__
