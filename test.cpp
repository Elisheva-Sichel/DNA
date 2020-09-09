#include <gtest/gtest.h>
#include "dna_sequence.h"

TEST(Dna_constructor,constructor)
{
	ASSERT_THROW(DnaSequence dna1(NULL),std::runtime_error);
}

TEST(Dna_constructor2,Invalid_Nucleotide)
{
	ASSERT_THROW(DnaSequence dna2("ASD"),InvalidNucleotide);
}

TEST(Dna_constructor3,Invalid_Nucleotide)
{
	DnaSequence dna3("GCA");
	ASSERT_STREQ(dna3.getDna().c_str(), "GCA");
}

TEST(Dna_constructor4,Invalid_Nucleotide)
{
	DnaSequence dna4("GCA");
	DnaSequence dna5=dna4;
	ASSERT_STREQ(dna5.getDna().c_str(), "GCA");
}

TEST(length,check_length)
{
	DnaSequence dna4("GCATTC");
	ASSERT_EQ(dna4.length(),6);
}

TEST(compare,comp)
{
	DnaSequence dna4("GCATTC");
	DnaSequence dna5("GCATTC");
	DnaSequence dna6("CCATTC");
	ASSERT_TRUE(dna4==dna5);
	ASSERT_TRUE(dna4!=dna6);
}

TEST(index,com)
{	
	DnaSequence dna4("GCATTC");
	dna4[3]='C';
	ASSERT_TRUE(dna4[2]=='A');
	ASSERT_TRUE(dna4[3]=='C');
}

TEST(slice,check_slice)
{
	DnaSequence dna4("GCA");
	DnaSequence dna5("GC");
	ASSERT_TRUE(dna4.slice(0,2)==dna5);
}
