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
// $Id: WLSSteppingAction.cc 75292 2013-10-30 09:25:15Z gcosmo $
//
/// \file optical/wls/src/WLSSteppingAction.cc
/// \brief Implementation of the WLSSteppingAction class
//
//
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4TouchableHandle.hh"

#include "WLSSteppingAction.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSSteppingActionMessenger.hh"
#include "WLSPhotonDetSD.hh"
#include "WLSStackingAction.hh"

#include "G4ParticleTypes.hh"

#include "WLSUserTrackInformation.hh"

#include "G4ProcessManager.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

// Purpose: Save relevant information into User Track Information

static const G4ThreeVector ZHat = G4ThreeVector(0.0,0.0,1.0);

G4int WLSSteppingAction::fMaxRndmSave = 10000;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//WLSSteppingAction::WLSSteppingAction(WLSDetectorConstruction* detector)
WLSSteppingAction::WLSSteppingAction(
   WLSDetectorConstruction* detector
   , WLSRunAction* RA
   , WLSEventAction* eventAction) // add
  : fDetector(detector), 
    fRunaction(RA), 
    fEventAction(eventAction) // add
{
  fSteppingMessenger = new WLSSteppingActionMessenger(this);

  fCounterEnd = 0;
  fCounterMid = 0;
  //fBounceLimit = 100000;
  fBounceLimit = 10000000;

  fOpProcess = NULL;
  ResetCounters();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingAction::~WLSSteppingAction()
{
  	delete fSteppingMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  WLSSteppingAction::SetBounceLimit(G4int i)   {fBounceLimit = i;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfBounces()      {return fCounterBounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfClad1Bounces() {return fCounterClad1Bounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfClad2Bounces() {return fCounterClad2Bounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::GetNumberOfWLSBounces()   {return fCounterWLSBounce;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSSteppingAction::ResetSuccessCounter()     {
      G4int temp = fCounterEnd; fCounterEnd = 0; return temp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void WLSSteppingAction::saveRandomStatus(G4String subDir) 
// save the random status into a sub-directory
// Pre: subDir must be empty or ended with "/"
{
 
    // don't save if the maximum amount has been reached
    if (WLSSteppingAction::fMaxRndmSave == 0) return;

    G4RunManager* theRunManager = G4RunManager::GetRunManager();
    G4String randomNumberStatusDir = theRunManager->GetRandomNumberStoreDir();
 
    G4String fileIn  = randomNumberStatusDir + "currentEvent.rndm";

    std::ostringstream os;

    os << "run" << theRunManager->GetCurrentRun()->GetRunID() << "evt"
       << theRunManager->GetCurrentEvent()->GetEventID() << ".rndm" << '\0';

    G4String fileOut = randomNumberStatusDir + subDir + os.str();

    G4String copCmd = "/control/shell cp "+fileIn+" "+fileOut;
    G4UImanager::GetUIpointer()->ApplyCommand(copCmd);

    WLSSteppingAction::fMaxRndmSave--;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void WLSSteppingAction::UserSteppingAction(const G4Step* theStep)
{
	//G4cout << "CALLED: WLSSteppingAction::UserSteppingAction" << G4endl;
  	G4Track* theTrack = theStep->GetTrack();
  	WLSUserTrackInformation* trackInformation
      = (WLSUserTrackInformation*) theTrack->GetUserInformation();

   G4String pname = theTrack->GetParticleDefinition()->GetParticleName();
	float  track_x = theTrack->GetPosition().x();
	float  track_y = theTrack->GetPosition().y();
	float  track_z = theTrack->GetPosition().z();
	//float  track_x = theTrack->GetVertexPosition().x();
	//float  track_y = theTrack->GetVertexPosition().y();
	//float  track_z = theTrack->GetVertexPosition().z();

   //
   // here, check the placement tree
   // 0 is the last placement object
   //
   //G4TouchableHandle aTouch = theTrack->GetTouchableHandle(); 
   G4TouchableHandle aTouch = theStep->GetTrack()->GetTouchableHandle(); 
   G4String vname = aTouch->GetVolume(0)->GetName();
   G4String vnumb = aTouch->GetVolume(0)->GetCopyNo();
   G4ThreeVector vpos1 = aTouch->GetVolume(0)->GetFrameTranslation();
   G4ThreeVector vpos2 = aTouch->GetVolume(0)->GetTranslation();
#if 0 
   if (vname!="World")
      G4cerr << "step volume name = " << vname  << " and copyNo = " << vnumb 
             << " pos : " << vpos1.x() << " " << vpos1.y() << " " << vpos1.z() 
             << " pos : " << vpos2.x() << " " << vpos2.y() << " " << vpos2.z() 
             << G4endl;
#endif

#if 0 //23/11/27
	if ( vname=="FbrScin")
      G4cerr << "this is FiberScin. the placement tree: upper one = "
             << aTouch->GetVolume(1)->GetName() << " copyNo = "
             << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
	// this is FiberScin. the placement tree: upper one = FbrCoat copyNo = 76

   if ( vname=="FbrCoat" ) {
      // this is coating, 1 must be Mother volume 
      // =>      this is FiberCoat. the placement tree: upper one = FiberLayerMother.1
      //G4cerr << "this is FiberCoat. the placement tree: upper one = " 
      //       << aTouch->GetVolume(1)->GetName() << " copyNo = "
      //       << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
      std::string fn = aTouch->GetVolume(1)->GetName();
      std::string::size_type idx = fn.rfind('.');
      std::string num = fn.substr(idx+1); // number
      int xxxx = atoi(num.c_str());
      G4cerr << "this is FiberCoat. the placement tree: upper one = " << aTouch->GetVolume(1)->GetName()
             << " copyNo = " << aTouch->GetVolume(1)->GetCopyNo()
             << " volume Num = " << xxxx << G4endl;
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.42 copyNo = 0 volume Num = 42
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.43 copyNo = 0 volume Num = 43
   }
#endif

#if 0 
   if ( vname=="Sct1" )
      // this is scintilator, 1 must be Coating of the fiber 
      // =>      this is FiberScin. the placement tree: upper one = Joint1 copyNo = 36
      G4cerr << "this is FiberScin. the placement tree: upper one = "
             << aTouch->GetVolume(1)->GetName() << " copyNo = "
             << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 0
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 12
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 24

   if ( vname=="Sct2" )
      // this is scintilator, 1 must be Coating of the fiber 
      // =>      this is FiberScin. the placement tree: upper one = Joint2 copyNo = 36
      G4cerr << "this is FiberScin. the placement tree: upper one = " 
             << aTouch->GetVolume(1)->GetName() << " copyNo = "
             << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 0
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 12
	// this is FiberScin. the placement tree: upper one = Joint2 copyNo = 24

   if ( vname=="Joint1" ) {
      // this is coating, 1 must be Mother volume 
      // =>      this is FiberCoat. the placement tree: upper one = FiberLayerMother.1
      //G4cerr << "this is FiberCoat. the placement tree: upper one = " 
      //       << aTouch->GetVolume(1)->GetName() << " copyNo = "
      //       << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
      std::string fn = aTouch->GetVolume(1)->GetName();
      std::string::size_type idx = fn.rfind('.');
      std::string num = fn.substr(idx+1); // number
      int xxxx = atoi(num.c_str());
      G4cerr << "this is FiberCoat. the placement tree: upper one = " << aTouch->GetVolume(1)->GetName() 
             << " copyNo = " << aTouch->GetVolume(1)->GetCopyNo()
             << " volume Num = " << xxxx << G4endl;
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.0 copyNo = 0 volume Num = 0
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.1 copyNo = 0 volume Num = 1
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.2 copyNo = 0 volume Num = 2
   }
   if ( vname=="Joint2" ) { 
      // this is coating, 1 must be Mother volume 
      // =>      this is FiberCoat. the placement tree: upper one = FiberLayerMother.1
      //G4cerr << "this is FiberCoat. the placement tree: upper one = " 
      //       << aTouch->GetVolume(1)->GetName() << " copyNo = "
      //       << aTouch->GetVolume(1)->GetCopyNo() << G4endl;
      std::string fn = aTouch->GetVolume(1)->GetName();
      std::string::size_type idx = fn.rfind('.');
      std::string num = fn.substr(idx+1); // number
      int xxxx = atoi(num.c_str());
      G4cerr << "this is FiberCoat. the placement tree: upper one = " << aTouch->GetVolume(1)->GetName() 
             << " copyNo = " << aTouch->GetVolume(1)->GetCopyNo() 
             << " volume Num = " << xxxx << G4endl;
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.0 copyNo = 0 volume Num = 0
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.1 copyNo = 0 volume Num = 1
	// this is FiberCoat. the placement tree: upper one = FiberLayerMother.2 copyNo = 0 volume Num = 2
	}
#endif

	//
   // take energy deposit and step length
	//
   float edep  = theStep->GetTotalEnergyDeposit();
   float stepl = theStep->GetStepLength();

  	G4StepPoint* prePoint = theStep->GetPreStepPoint (); // pre step
  	G4StepPoint* pstPoint = theStep->GetPostStepPoint(); // pst step
   float prePointTime = prePoint->GetGlobalTime();
   float pstPointTime = pstPoint->GetGlobalTime();

	// get GetTouchableHandle for each
	G4TouchableHandle preTouch = theStep->GetPreStepPoint ()->GetTouchableHandle();
	G4TouchableHandle pstTouch = theStep->GetPreStepPoint ()->GetTouchableHandle();

  	//G4VPhysicalVolume* prePhysVol = prePoint->GetPhysicalVolume();
  	//G4VPhysicalVolume* pstPhysVol = pstPoint->GetPhysicalVolume();
  	G4VPhysicalVolume* prePhysVol = preTouch->GetVolume();
  	G4VPhysicalVolume* pstPhysVol = pstTouch->GetVolume();

  	G4String prePhysVolname;// = prePhysVol->GetName();
  	G4String pstPhysVolname;// = pstPhysVol->GetName();
	int Preid=-1;// = prePhysVol->GetCopyNo();
	int Pstid=-1;// = prePhysVol->GetCopyNo();
	int Motherid=-1;

   if ( prePhysVol ) {
      prePhysVolname = prePhysVol->GetName();
      //Preid = prePhysVol->GetCopyNo();
      #if 1
      if ( prePhysVolname=="Sct1" || prePhysVolname=="Sct2" || 
           prePhysVolname=="FbrScin" ) { // scintillator
			// take mother id and copyNum
         std::string fn  = preTouch->GetVolume(2)->GetName(); // scinti - coating - mother - world
         std::string::size_type idx = fn.rfind('.');
         std::string num = fn.substr(idx+1);
         Motherid = atoi(num.c_str());
         Preid = preTouch->GetVolume(1)->GetCopyNo(); 
         //Preid = preTouch->GetVolume(1)->GetCopyNo() + nfibers_in_layer * Preid;
         //G4cerr << ">>>>> FiberSci Motherid = " << Motherid << " Preid = " << Preid << " edep = " << edep << " fn = " << fn << G4endl;
			// >>>>> FiberSci Motherid = 16 Preid = 76 edep = 0.0225326 fn = FiberLayerMother.16
      }
      #endif
   }
   if ( pstPhysVol ) {
      pstPhysVolname = pstPhysVol->GetName();
      //Pstid = pstPhysVol->GetCopyNo();
      #if 1
      if ( pstPhysVolname=="Sct1" || pstPhysVolname=="Sct2" || 
           pstPhysVolname=="FbrScin" ) { // scintillator
         std::string fn = pstTouch->GetVolume(2)->GetName(); // scinti - coating - Ymother - world
         std::string::size_type idx = fn.rfind('.');
         std::string num = fn.substr(idx+1);
         Motherid = atoi(num.c_str());
         Pstid = pstTouch->GetVolume(1)->GetCopyNo();
         //Pstid = pstTouch->GetVolume(1)->GetCopyNo() + nfibers_in_layer * Pstid;
         //G4cerr << ">>>>> FiberSci Motherid = " << Motherid << " Pstid = " << Pstid << " edep = " << edep << " fn = " << fn << G4endl;
			// >>>>> FiberSci Motherid = 16 Pstid = 76 edep = 0.0225326 fn = FiberLayerMother.16
      }
      #endif
   }

	//
   // just output
	//
	if ( 0 ) {
	//if ( edep!=0 && (prePhysVolname=="Sct1"     || prePhysVolname=="Sct2" ) ) {
   	G4cerr << "### fEventAction->GetEventID() = " << fEventAction->GetEventID() << G4endl; 
      G4cerr << "### name = " << pname << ", edep = " << edep << G4endl;
   	G4cerr << "x = " << track_x << G4endl;
   	G4cerr << "y = " << track_y << G4endl;
   	G4cerr << "z = " << track_z << G4endl;
		G4cerr << "prePointTime = " << prePointTime << " pstPointTime = " << pstPointTime << G4endl;	
   	G4cerr << "### prePhysVolname : " << prePhysVolname << "  its Preid = " << Preid << G4endl;
   	G4cerr << "### pstPhysVolname : " << pstPhysVolname << "  its Pstid = " << Pstid << G4endl;
	}

	//
  	// Recording data on a primary track : GetParentID()==0
	//
  	if ( theTrack->GetParentID()==0 ) {
      //G4cerr << "Name: " << theTrack->GetParticleDefinition()->GetParticleName() << G4endl;
      //G4cerr << "GetParentID()==0: PhysicalVolume : " << prePhysVolname << G4endl;
      //G4cerr << "GetParentID()==0: PhysicalVolume : " << pstPhysVolname << G4endl;
      if ( prePhysVolname == "ScCube" || 
           pstPhysVolname == "ScCube" ) {
      	fEventAction->AddEnergyDepositAndStepLength(edep, stepl);
		}
     	if ( theTrack->GetCurrentStepNumber() == 1 ) {
			//G4double pz = theTrack->GetVertexMomentumDirection().z();
			//G4double fInitTheta = theTrack->GetVertexMomentumDirection().angle(ZHat);
     	}
  	}

	#define __ROOTRUNACTION__

   #ifdef __ROOTRUNACTION__
	int evt, detid(0), pcode, trackid;
	float x,y,z;
   if ( !fRunaction->GetTreeStpAct() ) {
		G4cerr << "##### DefineTreeStpAct #####" << G4endl;
      fRunaction->DefineTreeStpAct();// = new TTree("tree", "");
      fRunaction->GetTreeStpAct()->Branch("evt",  &evt,  "evt/I");
      fRunaction->GetTreeStpAct()->Branch("detid",&detid,"detid/I");
      fRunaction->GetTreeStpAct()->Branch("trackid",&trackid,"trackid/I");
      fRunaction->GetTreeStpAct()->Branch("code",&pcode,"code/I");
      fRunaction->GetTreeStpAct()->Branch("x",   &x,    "x/F");
      fRunaction->GetTreeStpAct()->Branch("y",   &y,    "y/F");
      fRunaction->GetTreeStpAct()->Branch("z",   &z,    "z/F");
      fRunaction->GetTreeStpAct()->Branch("edep",&edep, "edep/F");
   }
	evt = fEventAction->GetEventID();
   pcode = theTrack->GetParticleDefinition()->GetPDGEncoding();
	trackid = theTrack->GetTrackID();
   #endif

	//
   // energy is stored in array... mother-id, fiber-id, edep 
	//
   if ( prePhysVolname=="Sct1" || prePhysVolname=="Sct2" ||
        prePhysVolname=="FbrScin" ) { // scintillator
      //G4cerr << "fill edepo" << G4endl;
		//fEventAction->AddEnergyDepositLayerScinti(Preid, edep);
		fEventAction->AddEnergyDepositLayerScinti(Motherid, Preid, edep);
		if ( Motherid%2 == 0 ) detid= 2;
		if ( Motherid%2 == 1 ) detid= 1;
	}

	#ifdef __ROOTRUNACTION__
	if ( pname!="optical photon" ) {
      x = track_x;
      y = track_y;
      z = track_z;
      fRunaction->GetTreeStpAct()->Fill();
	}
   #endif


	//
  	// Retrieve the status of the optical photon
	//
  	G4OpBoundaryProcessStatus theStatus = Undefined;
  	G4ProcessManager* OpManager 
		= G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  	if ( OpManager ) {
     	int MAXofPostStepLoops = OpManager->GetPostStepProcessVector()->entries();
     	G4ProcessVector* fPostStepDoItVector = OpManager->GetPostStepProcessVector(typeDoIt);

     	for (int i=0; i<MAXofPostStepLoops; i++) {
         G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
         fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
         if (fOpProcess) { 
				theStatus = fOpProcess->GetStatus(); 
				break;
			}
     	}
  	}

  	// Find the skewness of the ray at first change of boundary
  	if ( fInitGamma == -1 &&
       (theStatus == TotalInternalReflection ||
        theStatus == FresnelReflection ||
        theStatus == FresnelRefraction)
        && trackInformation->isStatus(InsideOfFiber) ) {

        G4double px = theTrack->GetVertexMomentumDirection().x();
        G4double py = theTrack->GetVertexMomentumDirection().y();
        G4double x_1  = theTrack->GetPosition().x();
        G4double y_1  = theTrack->GetPosition().y();

        fInitGamma = x_1 * px + y_1 * py;
        fInitGamma = fInitGamma / std::sqrt(px*px + py*py) / std::sqrt(x_1*x_1 + y_1*y_1);
        fInitGamma = std::acos(fInitGamma*rad);

        if ( fInitGamma / deg > 90.0)  { fInitGamma = 180 * deg - fInitGamma;}
  	}
  	// Record Photons that missed the photon detector but escaped from readout
  	if ( !pstPhysVol && trackInformation->isStatus(EscapedFromReadOut) ) {
		//UpdateHistogramSuccess(pstPoint,theTrack);
     	ResetCounters();
     	return;
  	}

  	// Assumed photons are originated at the fiber OR
  	// the fiber is the first material the photon hits
  	switch ( theStatus ) {
     	// Exiting the fiber
     	case FresnelRefraction:
     	case SameMaterial:

     	G4bool isFiber;
      isFiber = pstPhysVolname == "WLSFiber"
             || pstPhysVolname == "Clad1"
             || pstPhysVolname == "Clad2";
 
     	if ( isFiber ) {
      	if (trackInformation->isStatus(OutsideOfFiber)) {
         	trackInformation->AddStatusFlag(InsideOfFiber);
			}
     	// Set the Exit flag when the photon refracted out of the fiber
      } else if (trackInformation->isStatus(InsideOfFiber)) {

      	// EscapedFromReadOut if the z position is the same as fiber's end
       	if (theTrack->GetPosition().z() == fDetector->GetWLSFiberEnd()) {
           	trackInformation->AddStatusFlag(EscapedFromReadOut);
            fCounterEnd++;
        	} else { // Escaped from side
            trackInformation->AddStatusFlag(EscapedFromSide);
            trackInformation->SetExitPosition(theTrack->GetPosition());
            //UpdateHistogramEscape(pstPoint,theTrack);

            fCounterMid++;
            ResetCounters();
         }

         trackInformation->AddStatusFlag(OutsideOfFiber);
      	trackInformation->SetExitPosition(theTrack->GetPosition());
      }

      return;
 
     	// Internal Reflections
     	case TotalInternalReflection:
 
      // Kill the track if it's number of bounces exceeded the limit
      if ( fBounceLimit > 0 && fCounterBounce >= fBounceLimit ) {
          theTrack->SetTrackStatus(fStopAndKill);
          trackInformation->AddStatusFlag(murderee);
          ResetCounters();
          G4cout << "\n Bounce Limit Exceeded" << G4endl;
          return;
      }
 
     	case FresnelReflection:

     	fCounterBounce++;
      if      ( prePhysVolname == "WLSFiber") fCounterWLSBounce++;
      else if ( prePhysVolname == "Clad1"   ) fCounterClad1Bounce++;
      else if ( prePhysVolname == "Clad2"   ) fCounterClad2Bounce++;
 
      // Determine if the photon has reflected off the read-out end
      if (theTrack->GetPosition().z() == fDetector->GetWLSFiberEnd()) {
      	if (!trackInformation->isStatus(ReflectedAtReadOut) &&
              trackInformation->isStatus(InsideOfFiber)) {
         	trackInformation->AddStatusFlag(ReflectedAtReadOut);

            if (fDetector->IsPerfectFiber() &&
                theStatus==TotalInternalReflection) {
            	theTrack->SetTrackStatus(fStopAndKill);
               trackInformation->AddStatusFlag(murderee);
 					//UpdateHistogramReflect(pstPoint,theTrack);
               ResetCounters();
            	return;
         	}
      	}
      }
      return;

     	// Reflection of the mirror
     	case LambertianReflection:
     	case LobeReflection:

     	case SpikeReflection:
       	// Check if it hits the mirror
       	if ( pstPhysVolname == "Mirror" ) {
         	trackInformation->AddStatusFlag(ReflectedAtMirror);
 			}
       	return;

     	case Detection: // Detected by a detector
			//G4cout << "\npstPhysVolname = " << pstPhysVolname << G4endl;
			//G4cerr << "pstPhysVolname=" << pstPhysVolname<< " Pstid=" <<Pstid << G4endl;
		   //G4cerr << " pstPhysVolname = " << pstPhysVolname << " Pstid = " << Pstid 
         //       << " prePointTime = " << prePointTime << " pstPointTime = " << pstPointTime << G4endl;
			if (pstPhysVolname=="PhotonDetX"||
             pstPhysVolname=="PhotonDetY"||
             pstPhysVolname=="PhotonDetZ") {
				ResetCounters();
            theTrack->SetTrackStatus(fStopAndKill);

            // take the fastest timing in the detector
            //  compare if (theFastDetectorResponce > pstPointTime ) 
            fEventAction->SetFastestDetectorTiming(pstPointTime);
 
            if (pstPhysVolname=="PhotonDetX") {
					fEventAction->AddPhotCountX(1,Pstid);
					fEventAction->SetFastestTimeX(pstPointTime, Pstid);
				}  
            if (pstPhysVolname=="PhotonDetY") {
					fEventAction->AddPhotCountY(1,Pstid);
					fEventAction->SetFastestTimeY(pstPointTime, Pstid);
				}
            if (pstPhysVolname=="PhotonDetZ") {
					fEventAction->AddPhotCountZ(1,Pstid);
					fEventAction->SetFastestTimeZ(pstPointTime, Pstid);
            }
				return;
			}
   		//G4cerr << "prePhysVolname =" << prePhysVolname << " Preid =" << Preid << G4endl;
   		//G4cerr << "pstPhysVolname=" << pstPhysVolname<< " Pstid=" <<Pstid << G4endl;
   		//G4cerr << "PhysName0=" << PhysName0 << " PhysName1=" << PhysName1 << G4endl;
   		//G4cerr << "aTouch->GetReplicaNumber(depth=0)=" << aTouch->GetReplicaNumber(0) << G4endl;
   		//G4cerr << "aTouch->GetReplicaNumber(depth=1)=" << aTouch->GetReplicaNumber(1) << G4endl;


       	// Check if the photon hits the detector and process the hit if it does
      	//if (pstPhysVolname=="PhotonDet") {
      	if (pstPhysVolname=="SensitiveDetector") {
          	G4SDManager* SDman = G4SDManager::GetSDMpointer();
          	G4String SDname="WLS/PhotonDet";
          	WLSPhotonDetSD* mppcSD = (WLSPhotonDetSD*)SDman->FindSensitiveDetector(SDname);

          	if (mppcSD) mppcSD->ProcessHits_constStep(theStep,NULL);
          	// Record Photons that escaped at the end
            //if (trackInformation->isStatus(EscapedFromReadOut))
            //    UpdateHistogramSuccess(pstPoint,theTrack);
          	// ----- Stop Tracking when it hits the detector's surface
          	ResetCounters();
          	theTrack->SetTrackStatus(fStopAndKill);
         	return;
      	}
      	break;

   	default: 
			break;
  	}
 
  	// Check for absorbed photons
  	if (theTrack->GetTrackStatus() != fAlive  &&
      trackInformation->isStatus(InsideOfFiber)) {
		//UpdateHistogramAbsorb(pstPoint,theTrack);
     	ResetCounters();
     	return;
	}	
 
}
