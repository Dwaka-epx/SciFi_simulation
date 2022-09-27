//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WLSRunAction.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/WLSRunAction.hh
/// \brief Definition of the WLSRunAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSRunAction_h
#define WLSRunAction_h 1

#include "globals.hh"
#include "g4root.hh"

#include "G4String.hh"
#include "G4UserRunAction.hh"

//#include "TFile.h"
//#include "TTree.h"

#include <fstream>

//#define NTHLAYER 201 
//#define NTHLAYER 9 
//#define NTHLAYER 5 
//#define NTHLAYER 3 
#define NTHLAYER 1 

class G4Run;
class WLSRunActionMessenger;
class WLSDetectorConstruction;

#include "TFile.h"
#include "TTree.h"
class TFile;
class TTree;

class WLSRunAction : public G4UserRunAction
{
  public:

    WLSRunAction(WLSDetectorConstruction*,G4String);
    //WLSRunAction(G4String);
    virtual ~WLSRunAction();


    void Set0ChargedSecondary (){ fEnergyCharged = 0; }
    double GetChargedSecondary (){ return fEnergyCharged; }
    void AddChargedSecondary (G4double ekin)
                 {fEnergyCharged += ekin; fNbCharged++;
                  if (ekin<fEmin[0]) fEmin[0] = ekin;
                  if (ekin>fEmax[0]) fEmax[0] = ekin;
                 };

    void Set0NeutralSecondary(){ fEnergyNeutral = 0; }
    double GetNeutralSecondary(){ return fEnergyNeutral; }
    void AddNeutralSecondary (G4double ekin)
                 {fEnergyNeutral += ekin; fNbNeutral++;
                  if (ekin<fEmin[1]) fEmin[1] = ekin;
                  if (ekin>fEmax[1]) fEmax[1] = ekin;
                 };

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void  SetRndmFreq(G4int val) { fSaveRndm = val; }
    G4int GetRndmFreq()          { return fSaveRndm; }

    inline void SetAutoSeed (const G4bool val) { fAutoSeed = val; }

	// to call file pointer
	std::ofstream outputf;

	//
	// event action
	//
	inline void   DefineTreeEvtAct1() { fTreeEvtAct1 = new TTree("treeEvtAct1", "tree evt act 1"); }
	TTree        *GetTreeEvtAct1()    { return fTreeEvtAct1; }

   inline void   DefineTreeEvtAct2() { fTreeEvtAct2 = new TTree("treeEvtAct2", "tree evt act 2"); }
   TTree        *GetTreeEvtAct2()    { return fTreeEvtAct2; }

	// stacking action
	inline void   DefineTreeStkAct() { fTreeStkAct = new TTree("treeStkAct", "tree stacking action"); }
	TTree        *GetTreeStkAct()    { return fTreeStkAct; }

	// stacking action
	inline void   DefineTreeStpAct() { fTreeStpAct = new TTree("treeStpAct", "tree stepping action"); }
	TTree        *GetTreeStpAct()    { return fTreeStpAct; }

	private:
 
   WLSRunActionMessenger* fRunMessenger;
   WLSDetectorConstruction *fDetConstruction;

   G4double fEnergyDeposit;
   G4double fTrackLength;
   G4double fEnergyCharged, fEnergyNeutral;
   G4double fEmin[2], fEmax[2];

   G4long   fNbSteps;
   G4int    fNbCharged, fNbNeutral;

   G4int fSaveRndm;
   G4bool fAutoSeed;
   G4String fName;

	TFile *fFile;
   TTree *fTreeEvtAct1;
   TTree *fTreeEvtAct2;
   TTree *fTreeStkAct;
   TTree *fTreeStpAct;
};

#endif
