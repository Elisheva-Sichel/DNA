#ifndef __MYDnaSequence_H__
#define __MYDnaSequence_H__

#include <map>
#include <string>
#include <list>
#include <fstream>

#include "invalid_nucleotide.h"

class DnaSequence 
{
private:
	class Nucleotides
	{
	public:
		Nucleotides(){}
		Nucleotides(const char& dna);
		Nucleotides& operator = (const char& dna);
		operator char();
		char get_char()const;
		
	private:
		char m_char;
	};
	Nucleotides* m_dna;
	static bool check(const char str);
	//const static std::map<char, char> nuc_map;
public:
	
	DnaSequence(const std::string& dna);
	DnaSequence(const char *dna);
	DnaSequence(std::fstream& myfile);
	DnaSequence(const DnaSequence& dna);
	~DnaSequence();
	std::string getDna()const;
	size_t length()const;
	DnaSequence& operator = (const DnaSequence& dna);
	DnaSequence::Nucleotides& operator [](const size_t ind)const;
	void writeToFile(const std::string& fileName);
	DnaSequence slice(size_t head,size_t tail);
	DnaSequence pairs();
	size_t find(const DnaSequence& sub);
	size_t count(const DnaSequence& sub);
	std::list<size_t> find_all(const DnaSequence& sub);
	std::list<DnaSequence> find_all_consensus();
	friend std::ostream& operator << (std::ostream& cout,const DnaSequence::Nucleotides& dna);
};

inline DnaSequence::DnaSequence(const std::string& dna)
{
	m_dna=new Nucleotides[strlen(dna.c_str())];
	for(size_t i=0;i<=strlen(dna.c_str());i++)
	{
		m_dna[i]=dna.c_str()[i];
	}
}

inline DnaSequence::DnaSequence(const char* dna)
{	
	if(dna==NULL)
		throw std::runtime_error("can not get null");	
	m_dna=new Nucleotides[strlen(dna)];
	try
	{
		for(size_t i=0;i<=strlen(dna);i++)
		{
			m_dna[i]=dna[i];
		}
	}
	catch(InvalidNucleotide& exp)
	{
		delete[] m_dna;
		throw;
	}

}

inline DnaSequence::DnaSequence(std::fstream& myfile):m_dna(NULL)
{
	std::string line;
	std::string sequenceFromFile;
	if (myfile.is_open())
	{
  		while (getline(myfile,line))
  		{
      			sequenceFromFile += line;
    		}
   		 myfile.close();
  	}
	else
	{
		throw std::ifstream::failure("opening file failed");
	}
	m_dna =  new Nucleotides[sequenceFromFile.length()];
	try
	{
		for(size_t i=0; i < sequenceFromFile.length(); i++)
			m_dna[i] = sequenceFromFile[i];
	}
	catch(InvalidNucleotide& e)
	{
		delete[] m_dna;
		m_dna = NULL;
		throw;
	}
}

inline DnaSequence::DnaSequence(const DnaSequence& dna):m_dna(new Nucleotides[dna.length()])
{
	try
	{
		for(size_t i=0;i<=dna.length();i++)
		{
			m_dna[i]=dna[i].get_char();
		}
	}
	catch(InvalidNucleotide& exp)
	{
		delete[] m_dna;
		throw;
	}
	
}

inline DnaSequence::~DnaSequence()
{
	delete[] m_dna;
	m_dna = NULL;
}

inline DnaSequence& DnaSequence::operator = (const DnaSequence& dna)
{
	if(this!=&dna)
	{
		delete[] m_dna;
		m_dna=new Nucleotides[dna.length()];
		for(size_t i=0;i<=dna.length();i++)
		{
			m_dna[i]=dna[i];
		}
	}
	return *this;
}

inline DnaSequence::Nucleotides& DnaSequence::operator [] (const size_t ind)const
{
	if(ind>length())
	{			
		throw std::out_of_range("index not in range");
	}		
	return m_dna[ind];
}

inline bool operator == (const DnaSequence& str1,const DnaSequence& dna)
{
	size_t min_len=str1.length()<dna.length()?str1.length():dna.length();
	for(size_t i=0;i<=min_len;i++)
	{
		if(str1[i].get_char()!=dna[i].get_char())
		{
			return false;
		}
	}
	return true;
}

inline bool operator != (const DnaSequence& str1,const DnaSequence& dna)
{
	return !(str1==dna);
}

inline std::ostream& operator << (std::ostream& cout,const DnaSequence& dna)
{
	for(size_t i=0;i<dna.length();i++)
	{
		std::cout<<dna[i].get_char();
	}
	return cout;
}

inline void DnaSequence::writeToFile(const std::string &fileName)
{
	std::ofstream myfile;
 	myfile.open(fileName.c_str());
	if (myfile.is_open())
	{
		for(size_t i=0;i<length();i++)
		{
			myfile<<m_dna[i].get_char();
		}
		myfile.close();
		return;
	}
	else
	{
		throw std::ifstream::failure("opening file failed");
	}
}

inline std::string DnaSequence::getDna()const
{
	std::string sequence;
	for(size_t i=0; i< length(); i++)
	sequence += (*this)[i].get_char();
	return sequence;
}

inline size_t DnaSequence::length()const
{
	size_t i=0;
	while(m_dna[i].get_char()!='\0')
		i++;
	return i;
}

inline DnaSequence DnaSequence::slice(size_t head,size_t tail)
{
	char* dna=new char[tail-head-1];
	
	for(size_t i=head;i<tail;i++)
	{
		
		dna[i-head]=m_dna[i].get_char();
	}
	DnaSequence d(dna);
	delete [] dna;
	return d;
	
}
inline DnaSequence DnaSequence::pairs()
{
	static std::map<char, char>nuc_ma;
	nuc_ma['A']='T';
	nuc_ma['T']='A';
	nuc_ma['G']='C';
	nuc_ma['C']='G';
	char* dna=new char[length()];
	for(size_t i=length();i>0;i--)
	{
		dna[length()-i]=nuc_ma[(m_dna[i-1].get_char())];
	}
	DnaSequence d(dna);
	delete [] dna;
	return d;
}

inline size_t DnaSequence::count(const DnaSequence& sub)
{	
	return find_all(sub).size();
}

inline DnaSequence::Nucleotides::Nucleotides(const char& dna):m_char(dna)
{
	try
	{
		if(check(dna)==false)
		{
			throw InvalidNucleotide();
		}
	}
	catch(InvalidNucleotide& exp)
	{
		throw;
	}
}

inline DnaSequence::Nucleotides& DnaSequence::Nucleotides::operator = (const char& dna)
{
	try
	{
		if(check(dna)==false)
		{
			throw InvalidNucleotide();
		}
		m_char=dna;
	}
	catch(InvalidNucleotide& exp)
	{
		throw;
	}
	return *this;
}

inline std::ostream& operator << (std::ostream& cout,const DnaSequence::Nucleotides& dna)
{
	std::cout<<dna.get_char();
	return cout;
}

inline DnaSequence::Nucleotides::operator char()
{
	return m_char;
}
inline char DnaSequence::Nucleotides::get_char()const
{
	return m_char;
}

#endif //__MYDnaSequence_H__
