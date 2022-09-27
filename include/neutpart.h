#ifndef _NEUTPART_H_
#define _NEUTPART_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

class NeutPart : public TObject {

 public:
  NeutPart();

  ~NeutPart(){};

  // Variables 

  Int_t          fPID;     // Particle ID                (PDG)
  Double_t       fMass;    // Mass of the particle       (MeV/c2)
  TLorentzVector fP;       // 4 momentum of the particle (MeV/c, MeV)
  
  Bool_t         fIsAlive; // Particle should be tracked or not
                           //       ( in the detector simulator )

  Int_t          fStatus;  // Status flag of this particle
						   /*   -2: Non existing particle
						        -1: Initial state particle
								 0: Normal
								 1: Decayed to the other particle
								 2: Escaped from the detector
								 3: Absorped
								 4: Charge exchanged
								 5: Pauli blocked
								 6: N/A
								 7: Produced child particles
								 8: Inelastically scattered
						   */

  TLorentzVector fPosIni;  // Initilal position in the nucleus
  TLorentzVector fPosFin;  // Final(current) position in the nucleus

  ClassDef(NeutPart, 1)
};

#endif // _NEUTPART_H_
