#include "ParticleSetExplicit.hpp"
#include <RandomSeqGenerator.h>

namespace qmcplusplus
{

template<class RANDOM>
void ParticleSet::randomizeFromSourceWithEngine(const ParticleSet& src, RANDOM& rng)
{
  const SpeciesSet& srcSpSet(src.getSpeciesSet());
  SpeciesSet& spSet(getSpeciesSet());
  int srcChargeIndx = srcSpSet.getAttribute("charge");
  int srcMemberIndx = srcSpSet.getAttribute("membersize");
  int ChargeIndex   = spSet.addAttribute("charge");
  int MemberIndx    = spSet.addAttribute("membersize");
  int Nsrc          = src.getTotalNum();
  int Nptcl         = getTotalNum();
  int NumSpecies    = spSet.TotalNum;
  int NumSrcSpecies = srcSpSet.TotalNum;
  //Store information about charges and number of each species
  std::vector<int> Zat, Zspec, NofSpecies, NofSrcSpecies, CurElec;
  Zat.resize(Nsrc);
  Zspec.resize(NumSrcSpecies);
  NofSpecies.resize(NumSpecies);
  CurElec.resize(NumSpecies);
  NofSrcSpecies.resize(NumSrcSpecies);
  for (int spec = 0; spec < NumSrcSpecies; spec++)
  {
    Zspec[spec]         = (int)round(srcSpSet(srcChargeIndx, spec));
    NofSrcSpecies[spec] = (int)round(srcSpSet(srcMemberIndx, spec));
  }
  for (int spec = 0; spec < NumSpecies; spec++)
  {
    NofSpecies[spec] = (int)round(spSet(MemberIndx, spec));
    CurElec[spec]    = first(spec);
  }
  int totQ = 0;
  for (int iat = 0; iat < Nsrc; iat++)
    totQ += Zat[iat] = Zspec[src.GroupID[iat]];
  app_log() << "  Total ion charge    = " << totQ << std::endl;
  totQ -= Nptcl;
  app_log() << "  Total system charge = " << totQ << std::endl;
  // Now, loop over ions, attaching electrons to them to neutralize
  // charge
  int spToken = 0;
  // This is decremented when we run out of electrons in each species
  int spLeft = NumSpecies;
  std::vector<PosType> gaussRand(Nptcl);
  makeGaussRandomWithEngine(gaussRand, rng);
  for (int iat = 0; iat < Nsrc; iat++)
  {
    // Loop over electrons to add, selecting round-robin from the
    // electron species
    int z = Zat[iat];
    while (z > 0 && spLeft)
    {
      int sp = spToken++ % NumSpecies;
      if (NofSpecies[sp])
      {
        NofSpecies[sp]--;
        z--;
        int elec = CurElec[sp]++;
        app_log() << "  Assigning " << (sp ? "down" : "up  ") << " electron " << elec << " to ion " << iat
                  << " with charge " << z << std::endl;
        double radius = 0.5 * std::sqrt((double)Zat[iat]);
        R[elec]       = src.R[iat] + radius * gaussRand[elec];
      }
      else
        spLeft--;
    }
  }
  // Assign remaining electrons
  int ion = 0;
  for (int sp = 0; sp < NumSpecies; sp++)
  {
    for (int ie = 0; ie < NofSpecies[sp]; ie++)
    {
      int iat       = ion++ % Nsrc;
      double radius = std::sqrt((double)Zat[iat]);
      int elec      = CurElec[sp]++;
      R[elec]       = src.R[iat] + radius * gaussRand[elec];
    }
  }
}

template void ParticleSet::randomizeFromSourceWithEngine<RandomBase<double>>(const ParticleSet& pset_source,
                                                                             RandomBase<double>& rng);

template void ParticleSet::randomizeFromSourceWithEngine<RNGThreadSafe<StdRandom<double>>>(
    const ParticleSet& pset_source,
    RNGThreadSafe<StdRandom<double>>& rng);

} // namespace qmcplusplus
