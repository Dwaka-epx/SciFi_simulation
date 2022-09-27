#ifndef _NEUTVTX_H_
#define _NEUTVTX_H_

#include <TObject.h>
#include <TObjArray.h>
#include <TLorentzVector.h>

class NeutVtx : public TObject {

 public:
  NeutVtx(Int_t np = 0);

  ~NeutVtx();

  Int_t EventNo;
  void  SetNvtx(Int_t nvtx);
  Int_t Nvtx(void)        const {return fNvtx; }

  TLorentzVector *Pos(Int_t vtxid) const
	{ return (TLorentzVector *)((fNvtx>vtxid) ? (fPos->At(vtxid)) : NULL);};
  
  void  SetPos(Int_t vtxid, TLorentzVector pos);
  void  SetPos(Int_t nvtx , TLorentzVector *pos);
  void  Dump();

 private:
  /* Generated vector information */
  Int_t fNvtx;         //          Number of particles

  TObjArray *fPos;     // ->       To store vertexes

  ClassDef(NeutVtx, 1)
};

#endif
