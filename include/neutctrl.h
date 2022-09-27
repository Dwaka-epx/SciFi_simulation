#ifndef _NEUTCTRL_H_
#define _NEUTCTRL_H_

#include <TObject.h>

class NeutCtrl : public TObject {

 public:

  NeutCtrl() {};
  ~NeutCtrl() {};

  /* Control flags */
  Bool_t IsFerm;      //          Consider Fermi momentum or not
                      //          True  : consider
                      //          False : ignore

  Bool_t IsPauli;     //          Consider Pauli blocking effect or not
                      //          True  : consider
                      //          False : ignore

  Int_t ModeSelect;   //          Primary neutrino interaction mode selection
                      //           0 : default
                      //          -1 : Use cross-section weight (CrsWgt)
                      //           N : default (N>0) : select one mode

  Float_t CrsWgt[48]; //          Cross-section weight table
                      //           From  0 to 23 : for neutrino
                      //           From 24 to 47 : for anti-neutrino
  
  /* Nuclear effect related */
  Bool_t IsNucEff;    //          Consider pi/K/p/n/eta rescattering or not
                      //          True  : consider
                      //          False : ignore

  Bool_t  IsFZone;    //          Formation zone flag
                      //          True  : on
                      //          False : off

  /* deep inelastic scattering related */
  Int_t PDFID;        //          Type of the parton distribution function.
                      //           7 : GRV94
                      //          12 : GRV98

  Bool_t IsBdkYn;     //          Use Bodek-Yang correction or not
                      //          True  : Use
                      //          False : Not use

  ClassDef(NeutCtrl, 1)

};
#endif
