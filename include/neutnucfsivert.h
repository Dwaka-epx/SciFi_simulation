#ifndef _NEUTNUCFSIVERT_H_
#define _NEUTNUCFSIVERT_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

class NeutNucFsiVert : public TObject {

 public:
  NeutNucFsiVert();

  ~NeutNucFsiVert(){};

  // Variables 

  Int_t          fVertFlag;       // Vertex flag
  Int_t          fVertFirstStep;  // ID of the First step
  TLorentzVector fPos;            // Vertex position 
  TLorentzVector fMom;            // 4 momentum
  
  ClassDef(NeutNucFsiVert, 1)
};

#endif // _NEUTNUCFSIVERT_H_
