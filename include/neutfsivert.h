#ifndef _NEUTFSIVERT_H_
#define _NEUTFSIVERT_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

class NeutFsiVert : public TObject {

 public:
  NeutFsiVert();

  ~NeutFsiVert(){};

  // Variables 

  TLorentzVector fPos;     // Position of vertex
  Int_t          fVertID;  // Interaction type at vertex (*10 FOR HI-NRG)
//                  -1 : ESCAPE
//                   0 : INITIAL (or unmatched parent vertex if I>1)
//                   1 :
//                   2 : 
//                   3 : ABSORPTION
//                   4 : CHARGE EXCHANGE
//                   5 : 
//                   6 : 
//                   7 : HADRON PRODUCTION (hi-nrg only, i.e. 70)
//                   8 : QUASI-ELASTIC SCATTER
//                   9 : FORWARD (ELASTIC-LIKE) SCATTER


  
  
  
  ClassDef(NeutFsiVert, 1)
};

#endif // _NEUTFSIVERT_H_
