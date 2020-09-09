#include "dna_sequence.h"
bool DnaSequence::check(const char str)
{
	if ((str!='A')&(str!='C')&(str!='G')&(str!='T')&(str!='\0'))
		return false;
	return true;
}

size_t DnaSequence::find(const DnaSequence& sub)
{
	size_t j=0,i=0;
	for(i=0;i<length()-sub.length()+1;i++)
	{
		if(m_dna[i].get_char()==sub[0].get_char())
		{
			size_t tmp=i;
			for(j=1;j<sub.length();j++)
			{
				i++;
				if(m_dna[i].get_char()!=sub[j].get_char())
				{
					i=tmp;
					break;
				}
			}
			if(j==sub.length())
			{
				
				return tmp;
			}
		}
	}
	if(i==length()-sub.length()+1)
	{
		throw std::out_of_range("there is not such a subsequence");		
	}
	return 0;
}
std::list<size_t> DnaSequence::find_all(const DnaSequence& sub)
{
	std::list<size_t> vec;
	size_t tmp=0;
	while(tmp<length()-sub.length())
	{
		try
		{
			tmp=tmp+slice(tmp,length()+1).find(sub)+1;
			vec.push_back(tmp-1);
		}
		catch(std::out_of_range& e)
		{
			return vec;
		}
	}
	return vec;
}
std::list<DnaSequence> DnaSequence::find_all_consensus()
{
	std::list<DnaSequence> sol;
	std::list<size_t> start=find_all("ATG");
	std::list<size_t> end1=find_all("TAG");
	std::list<size_t> end2=find_all("TAA");
	std::list<size_t> end3=find_all("TGA");
	std::list<size_t>::iterator it_s;
	std::list<size_t>::iterator it_e;
	for(it_s = start.begin(); it_s != start.end(); ++it_s)
	{
		for(it_e = end1.begin(); it_e != end1.end(); ++it_e)
		{
			
			if((((*it_e-*it_s)%3)==0)&(*it_e>*it_s))
			{
				sol.push_back(slice(*it_s,*it_e));
			}
		}
		for(it_e = end2.begin(); it_e != end2.end(); ++it_e)
		{
			if((((*it_e-*it_s)%3)==0)&(*it_e>*it_s))
			{
				sol.push_back(slice(*it_s,*it_e));
			}
		}
		for(it_e = end3.begin(); it_e != end3.end(); ++it_e)
		{
			if((((*it_e-*it_s)%3)==0)&(*it_e>*it_s))
			{

				sol.push_back(slice(*it_s,*it_e));
			}
		}
	}
	return sol;
}
