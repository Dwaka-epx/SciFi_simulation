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
// $Id: WLSStackingAction.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/WLSStackingAction.cc
/// \brief Implementation of the WLSStackingAction class
//
//
#include "WLSRunAction.hh"
#include "WLSEventAction.hh"
#include "WLSStackingAction.hh"

#include "G4RunManager.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


WLSStackingAction::WLSStackingAction(WLSRunAction* RA, WLSEventAction* EA) 
: fRunaction(RA), fEventaction(EA),
  fPhotonCounter(0) 
{ 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSStackingAction::~WLSStackingAction() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
WLSStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{

	int q = 0;
	double cene = 0;
	double nene = 0; 
	#if 0 
   if ( !fRunaction->GetTreeStkAct() ) {
      //std::cerr << "tree is defined in evt = " << nevt << std::endl;
      fRunaction->DefineTreeStkAct();// = new TTree("tree", "");
      fRunaction->GetTreeStkAct()->Branch("q",&q,"q/I");
      fRunaction->GetTreeStkAct()->Branch("cene",&cene,"cene/D");
      fRunaction->GetTreeStkAct()->Branch("nene",&nene,"nene/D");
	}
	#endif
  	G4ParticleDefinition* particleType = aTrack->GetDefinition();

   // keep primary particle
   if ( aTrack->GetParentID() == 0 ) return fUrgent;

   //
   //energy spectrum of secondaries
   //
   G4double energy = aTrack->GetKineticEnergy();
   bool charged = (aTrack->GetDefinition()->GetPDGCharge() != 0.);
   #if 0 
   G4cerr << "\n[StackingAction]: " 
          <<" charged = " << charged
          //<< " edep = " << std::setw(13) << edep
          << " kin_ene = " << std::setw(13) << energy 
          << "\n" << G4endl;
   #endif
	#if 1
   //fEventaction->AddSecondary(energy);  
   if (charged) {
		q = 1;
		cene    = energy;
      //fRunaction->AddChargedSecondary(energy);
      //analysisManager->FillH1(4,energy);
   } else {
		q = 0;
		nene    = energy;
      //fRunaction->AddNeutralSecondary(energy);
      //analysisManager->FillH1(5,energy);
   }
   //fRunaction->GetTreeStkAct()->Fill();	
	#endif

  	if ( particleType == G4OpticalPhoton::OpticalPhotonDefinition() ) {
   	// keep optical photon
     	fPhotonCounter++;
     	return fUrgent;
  	} else {
     	// discard all other secondaries
     	// return fKill;
  	}
  	return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSStackingAction::NewStage() 
{
   G4cout << "\n\n##### Number of optical photons produces in this event : "
          << fPhotonCounter << " #####\n\n" << G4endl;
}

int WLSStackingAction::GetOpticalNPhotons() 
{
	return fPhotonCounter;	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSStackingAction::PrepareNewEvent() { fPhotonCounter = 0; }
