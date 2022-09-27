#ifndef _NEUTNUCFSISTEP_H_
#define _NEUTNUCFSISTEP_H_

#include <TObject.h>

class NeutNucFsiStep : public TObject {

 public:
  NeutNucFsiStep();

  ~NeutNucFsiStep(){};

  // Variables 

  Double_t       fECMS2; // CMS energy squared
  Double_t       fProb;  // Probability
  
  ClassDef(NeutNucFsiStep, 1)
};

#endif // _NEUTNUCFSISTEP_H_
