
#ifndef _N_NEUTROOTTREESINGLETON_H_
#define _N_NEUTROOTTREESINGLETON_H_

#include "TTree.h"
#include "TChain.h"
#include "neutvect.h"

 class NeutrootTreeSingleton { 

 public:
   static NeutrootTreeSingleton* Instance (std::string filename);
   static NeutrootTreeSingleton* Instance (TTree *a_intree);

   NeutrootTreeSingleton();
   ~NeutrootTreeSingleton();

   NeutVect* GetNeutVectAddress();
   Int_t GetEntry(Long64_t entry = 0, Int_t getall = 0);
   void LoadTree(std::string filename);
   void LoadTree(TTree *a_intree);

   TChain *tree_neutroot;
   TBranch *br_neutvect;
   NeutVect *nvect;

   int f_nEvents;
   int f_nFiles;

   static NeutrootTreeSingleton * fInstance;

   struct Cleaner {
     void DummyMethodAndSilentCompiler() { }
     ~Cleaner() {
       if (NeutrootTreeSingleton::fInstance !=0) {
	 delete NeutrootTreeSingleton::fInstance;
	 NeutrootTreeSingleton::fInstance = 0;
       }
     }
   };
   friend struct Cleaner;

   ClassDef(NeutrootTreeSingleton, 1)    
};

#endif

