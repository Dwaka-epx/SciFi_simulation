#ifndef _NEUTFSIPART_H_
#define _NEUTFSIPART_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

class NeutFsiPart : public TObject {

 public:
  NeutFsiPart();

  ~NeutFsiPart(){};

  // Variables 

  Int_t          fPID;     // Particle ID                (PDG)
  TLorentzVector fDir;     // Direction of particle
  Float_t        fMomLab;  // Momentum of particle in lab frame (MeV/c)
  Float_t        fMomNuc;  // Momentum of particle in nucleon rest frame (MeV/c)
  Int_t          fVertStart; // Starting vertex index in NeutFsiVert array
  Int_t          fVertEnd;  // Ending vertex index in NeutFsiVert array
  
  
  
  ClassDef(NeutFsiPart, 1)
};

#endif // _NEUTFSIPART_H_
