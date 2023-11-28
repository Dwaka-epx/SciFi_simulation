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
// $Id: WLSPrimaryGeneratorAction.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/WLSPrimaryGeneratorAction.hh
/// \brief Definition of the WLSPrimaryGeneratorAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSPrimaryGeneratorAction_h
#define WLSPrimaryGeneratorAction_h 1

#include "sizeOfFiberArray.hh"
#include "globals.hh"
//#include "WLSEventAction.hh"
#include "WLSRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4AffineTransform.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

class G4GeneralParticleSource;
class G4ParticleGun;

class G4Event;
class G4PhysicsTable;

class WLSDetectorConstruction;
class WLSPrimaryGeneratorMessenger;

// ROOT
#include "TFile.h"
#include "TTree.h"
class TFile;
class TTree;

// from NEUT
#include "neutvect.h"
#include "neutpart.h"
#include "neutvtx.h"


class WLSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:

		//WLSPrimaryGeneratorAction(WLSDetectorConstruction*);
		//WLSPrimaryGeneratorAction(WLSDetectorConstruction*, WLSEventAction*);
	//
	// TString has some linker trouble: 
	 //  undefined : TString::TString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const&)
	 //  => use G4String
		//WLSPrimaryGeneratorAction(WLSDetectorConstruction*, G4String, int, int);
		WLSPrimaryGeneratorAction(WLSDetectorConstruction*, G4String, int, int, WLSRunAction*);
		virtual ~WLSPrimaryGeneratorAction();

	public:

		virtual void GeneratePrimaries(G4Event*);
		void BuildEmissionSpectrum();
		void SetOptPhotonPolar(G4double);
		void SetDecayTimeConstant(G4double);

		//const G4GeneralParticleSource* GetSouce();
		//const G4ParticleGun* GetSouce();
		G4ParticleGun* GetSouce();
		int GetNeutrinoMode()         { return _neutMode;    }
		int GetNeutrinoNSeed()        { return _neutNSeeds;  }
		float GetNeutrinoEnergy()     { return _neutEnergy;  }
		float GetNeutrinoMomentum()   { return _neutMomentum;}
		
		int GetNParticles()            { return _nParticles; } //add 23/11/27 note:_nGenerate < _nParticles
		int GetNNeutrons()            { return _nNeutrons; }
		int GetNMuons()            { return _nMuons; }
		int GetNProtons()             { return _nProtons;  }
		int GetNPions()               { return _nPions;    }
		void SetNeutrinoMode(int a)   { _neutMode   = a; }
		void SetNeutrinoNSeed(int a)  { _neutNSeeds = a; }
		int fPID[nmaxParticles];


	protected:

		G4PhysicsTable* fIntegralTable;

	private:

		void SetOptPhotonPolar();
		void SetOptPhotonTime();
 
		WLSRunAction* fRunAction;
		WLSDetectorConstruction*   fDetector;
		//G4GeneralParticleSource*   fParticleGun;
		G4ParticleGun*              fParticleGun;
		WLSPrimaryGeneratorMessenger* fGunMessenger;

		static G4bool fFirst;
		G4double fTimeConstant;

		//WLSEventAction* fEventAction; // add
	 TTree *_NeutTree;
	 NeutVect *_NeutVec;
	 G4String _NeutName;
	 int _NeutSttNum;
	 int _NeutGenNum;
	 int _neutMode;
	 int _neutNSeeds;
	 float _neutEnergy;
	 float _neutMomentum;

	 int _nParticles; // all particles //add 23/11/27
	 int _nNeutrons; // neutron
	 int _nProtons; // proton
	int _nPions; // pion
	int _nMuons;
};

#endif
