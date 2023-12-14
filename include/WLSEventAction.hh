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
// $Id: WLSEventAction.hh 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/include/WLSEventAction.hh
/// \brief Definition of the WLSEventAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef WLSEventAction_h
#define WLSEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"

#include "WLSPrimaryGeneratorAction.hh"
#include "WLSStackingAction.hh"
//#include "WLSRunAction.hh"

//define nlayers_dummy, nfibers_dummy
#include "MyConst.hh"

class WLSRunAction;
class WLSEventActionMessenger;
class WLSPrimaryGeneratorAction;
class WLSStackingAction;

//const int nlayers_dummy =1200;
//const int nfibers_dummy =1200;

class WLSEventAction : public G4UserEventAction
{
  public:

    //WLSEventAction(WLSRunAction*,);
    WLSEventAction(WLSRunAction*,WLSPrimaryGeneratorAction*);
    //WLSEventAction(WLSRunAction*,WLSPrimaryGeneratorAction*,WLSStackingAction*);
    virtual ~WLSEventAction();

  public:

    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);

    G4int GetEventNo();
    void SetEventVerbose(G4int);

    void SetDrawFlag(G4String val)  { fDrawFlag = val; };
    void SetPrintModulo(G4int val)  { fPrintModulo = val; };

    void SetForceDrawPhotons(G4bool b){fForceDrawPhotons=b;}
    void SetForceDrawNoPhotons(G4bool b){fForceNoPhotons=b;}

    void SetBeamPrimaryX(G4int a) { fPrimaryX = a; } // add
    void SetBeamPrimaryY(G4int a) { fPrimaryY = a; } // add
    void SetBeamPrimaryZ(G4int a) { fPrimaryZ = a; } // add

   inline void AddPhotCountX(G4int a, G4int ch) { fPhotCountX[ch] += a; } // add
   inline void AddPhotCountY(G4int a, G4int ch) { fPhotCountY[ch] += a; } // add
   inline void AddPhotCountZ(G4int a, G4int ch) { fPhotCountZ[ch] += a; } // add

   //inline void AddEnergyDepositLayerScintiX(int m, int a, float b) { fEdepoLayerScintiX[m][a] += b; } 
   //inline void AddEnergyDepositLayerScintiY(int m, int a, float b) { fEdepoLayerScintiY[m][a] += b; } 
   //inline void AddEnergyDepositLayerScinti     (int a, float b) { fEdepoLayerScinti     [a] += b; } 
   inline void AddEnergyDepositLayerScinti  (int m, int a, float b) { 
		fEdepoLayerScinti[m][a] += b; 
	} 

   inline void AddEnergyDepositAndStepLength(float a, float b) { 
		fEdep += a; 
		fStepLength += b; 
	}
 
 
	inline void SetFastestDetectorTiming(G4float time) {
   	if ( fFastestDetResponce > time ) fFastestDetResponce = time;   
	}

   inline void SetFastestTimeX(G4float time, G4int ch) { 
		if ( fFastestTimeX[ch] > time ) fFastestTimeX[ch] = time; 
	} // add
   inline void SetFastestTimeY(G4float time, G4int ch) { 
		if ( fFastestTimeY[ch] > time ) fFastestTimeY[ch] = time; 
	} // add
   inline void SetFastestTimeZ(G4float time, G4int ch) {
		if ( fFastestTimeZ[ch] > time ) fFastestTimeZ[ch] = time; 
	} // add

	void GiveParticleInitialPosi(G4ThreeVector a);

	int GetEventID() { return fEventID; }
  	private:

    WLSRunAction* fRunAction;
    WLSEventActionMessenger* fEventMessenger;
	 WLSPrimaryGeneratorAction *fPrimarysource;
	 WLSStackingAction* fStacking;

    G4int fVerboseLevel;
    G4int fPrintModulo;
 
    G4int fMPPCCollID;

	 int fEventID; 
    G4String fDrawFlag;

    G4bool fForceDrawPhotons;
    G4bool fForceNoPhotons;
	
	int fPrimaryX; // add
	int fPrimaryY; // add
	int fPrimaryZ; // add
	int fPhotCountX[NTHLAYER*NTHLAYER]; // add
	int fPhotCountY[NTHLAYER*NTHLAYER]; // add
	int fPhotCountZ[NTHLAYER*NTHLAYER]; // add
   
   //float fEdepoLayerScintiX[999][299]; // add
   //float fEdepoLayerScintiY[999][299]; // add
   // const int nlayers_dummy = 999;
   // const int nfibers_dummy = 999;
   float fEdepoLayerScinti [nlayers_dummy][nfibers_dummy]; // add

   float fEdep;
	float fStepLength;
   float fFastestDetResponce;

	float fFastestTimeX[NTHLAYER*NTHLAYER]; // add
	float fFastestTimeY[NTHLAYER*NTHLAYER]; // add
	float fFastestTimeZ[NTHLAYER*NTHLAYER]; // add
};

#endif
