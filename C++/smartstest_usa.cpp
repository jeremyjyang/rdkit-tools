#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

int main(int argc, char *argv[])
{
  std::vector<std::string> smis;
  smis.push_back("C1CC=CCC1");
  smis.push_back("C1CCCCC1");
  smis.push_back("CC1CCC3C4CCCCC4CCC3C1");

  std::vector<std::string> smas;
  smas.push_back("*~1~*~*~*~*~*1");
  smas.push_back("*~1~*~*~*~*~*~1");
  smas.push_back("*~1~*~*~*~*~*~*~*~*~*~*~*~*~*1");

  for (int i=0;i<smis.size();++i)
  {
    std::cerr << "smi: " << smis[i] << std::endl;
    for (int j=0;j<smas.size();++j)
    {
      std::cerr << "\tsma: " << smas[j] << std::endl;

      ROMol *mol=SmilesToMol(smis[i]);
      ROMol *pattern=SmartsToMol(smas[j]);

      // a MatchVect is a vector of std::pairs with (patIdx,molIdx):
      std::vector<MatchVectType> matches;
      unsigned int nMatches;
      // 4th arg is "uniquify".
      nMatches=SubstructMatch(*mol,*pattern,matches,false,true,false,false);
      std::cerr << "\t\tmatches.size() = " << matches.size() << std::endl;
      if (nMatches==0) continue;
      nMatches=SubstructMatch(*mol,*pattern,matches,true,true,false,false);
      std::cerr << "\t\tunique matches.size() = " << matches.size() << std::endl;

      for (int k=0;k<matches.size();++k)
      {
        std::cerr << "\t\t( " ;
        for (int kk=0;kk<matches[k].size();++kk)
        {
          std::cerr << matches[k][kk].second ;
          if (kk<matches[k].size()-1) std::cerr << "," ;
        }
        std::cerr << " )" << std::endl;
      }

      delete pattern;
      delete mol;
    }
    std::cerr << std::endl;
  }
}

