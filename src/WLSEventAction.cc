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
// $Id: WLSEventAction.cc 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/src/WLSEventAction.cc
/// \brief Implementation of the WLSEventAction class


#include "WLSEventAction.hh"
#include "WLSRunAction.hh"

#include "WLSEventActionMessenger.hh"

#include "WLSPhotonDetHit.hh"
#include "WLSTrajectory.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"

#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"

#include "Randomize.hh"

#include "TMath.h"

// Purpose: Invoke visualization at the end
//          Also can accumulate statistics regarding hits
//          in the PhotonDet detector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventAction::WLSEventAction(WLSRunAction* runaction, 
                               WLSPrimaryGeneratorAction* primarysource)
 : /* initialize with different name */ fRunAction(runaction), 
   /* initialize with different name */ fPrimarysource(primarysource), 
   /* initialize with different name */ 
   fVerboseLevel(0),
   fPrintModulo(100), 
   fDrawFlag("all")
{
  	fEventMessenger = new WLSEventActionMessenger(this);
  	fForceDrawPhotons = false;
  	fForceNoPhotons   = false;
  	fMPPCCollID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSEventAction::~WLSEventAction()
{
  	delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventAction::BeginOfEventAction(const G4Event* evt)
{
 	int evtNb = evt->GetEventID();
   fEventID  = evt->GetEventID();

 	if (evtNb%fPrintModulo==0) {
    	G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
	}
 	if (fVerboseLevel>0) {
    	G4cout << "<<< Event  " << evtNb << " started." << G4endl;
		G4cout << "<<< NTHLAYER  " << NTHLAYER << G4endl;
	}

	fPrimaryX=0;
	fPrimaryY=0;
	fPrimaryZ=0;
   fEdep = 0;
	fStepLength = 0;
   fFastestDetResponce = 99999;
	for (int i = 0; i < NTHLAYER * NTHLAYER; i++) {
		fPhotCountX[i]=0;
		fPhotCountY[i]=0;
		fPhotCountZ[i]=0;
      fFastestTimeX[i]=99999;
      fFastestTimeY[i]=99999;
      fFastestTimeZ[i]=99999;
	}

	G4cerr << ">>>>> clear data array" << G4endl;
   for (int m = 0; m < nlayers_dummy; m++) {
   	for (int i = 0; i < nfibers_dummy; i++) {
			//fEdepoLayerScintiX[m][i]=0;
			//fEdepoLayerScintiY[m][i]=0;
			fEdepoLayerScinti[m][i]=0;
		}
	}

	fRunAction->Set0ChargedSecondary();
	fRunAction->Set0NeutralSecondary();
}

/*
void WLSEventAction::GiveParticleInitialPosi(G4ThreeVector a)
{
	G4cerr << "\n\nCALLED: WLSEventAction::GiveParticleInitialPosi" << G4endl;
   G4cerr << "a.getX()=" << a.getX() << " a.getY()=" << a.getY() << G4endl;
	fPrimaryX=a.getX();
	fPrimaryY=a.getY();
	fPrimaryZ=a.getZ();
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void WLSEventAction::EndOfEventAction(const G4Event* evt)
{
   int ievt  = evt->GetEventID(); // #of event
  	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   
  	//Visualization of Trajectory
  	if ( pVVisManager ) {
   	G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   	int nTtrajects = 0;
   	if ( trajectoryContainer ) nTtrajects = trajectoryContainer->entries();
		if ( fDrawFlag == "all"    ) G4cout << "draw all trajectories" << G4endl;
   	if ( fDrawFlag == "charged") G4cout << "draw only charged" << G4endl;
   	G4cout << " #of trajects: " << nTtrajects << G4endl;
   	for (int i = 0; i < nTtrajects; i++) { 
			WLSTrajectory* trj = (WLSTrajectory *)((*(evt->GetTrajectoryContainer()))[i]);
         //G4cout << "call DrawTrajectory" << " Particle Name: " << trj->GetParticleName() << G4endl;
	      trj->SetForceDrawTrajectory(true);		
         trj->SetForceNoDrawTrajectory(false);
         #if 0
         
			G4cout << "IF VISULARIZATION IS ON :" << G4endl;
   		G4cout << "<<< fPrimary=Name       " << trj->GetParticleName()           << G4endl; // add
   		G4cout << "<<< fPrimary=momentum.X " << trj->GetInitialMomentum().getX() << G4endl; // add
   		G4cout << "<<< fPrimary=momentum.Y " << trj->GetInitialMomentum().getY() << G4endl; // add
   		G4cout << "<<< fPrimary=momentum.Z " << trj->GetInitialMomentum().getZ() << G4endl; // add
   		G4cout << "<<< fPrimary=momentum.M " << trj->GetInitialMomentum().mag()  << G4endl; // add
			#endif
        	if      ( fDrawFlag == "all" ) {
				trj->DrawTrajectory();
        	} else if ( fDrawFlag == "charged" && trj->GetCharge()!= 0.) {
         	trj->DrawTrajectory();
			}
    	}
  	}

	#if 0

   //Get Initial position, momentum, PID 23/11/27
   struct ParticleInfo {
    Int_t pid;                 // 粒子の識別子
    Double_t initialMomentumX; // 初期運動量の x 成分
    Double_t initialMomentumY; // 初期運動量の y 成分
    Double_t initialMomentumZ; // 初期運動量の z 成分
    Double_t initialPositionX; // 初期位置の x 座標
    Double_t initialPositionY; // 初期位置の y 座標
    Double_t initialPositionZ; // 初期位置の z 座標
   };
   std::cerr << ">>>>> struct ParticleInfo is defined in evt = " << ievt <<" <<<<<\n\n"<< std::endl;

   ParticleInfo particleInfo;
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   //G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
         std::cerr << "trajectoryContainer =" << trajectoryContainer<< std::endl;

   if( trajectoryContainer ){
         std::cerr << ">>>>> trajectoryContainer is defined in evt = " << ievt << std::endl;
      if(!fRunAction->GetTreeMLInfo()){
         //TTree* TreeMLInfo = new TTree("treeMLInfo", "Information for ML such as initialPosition, initialMomentum, PID and so on.");
         std::cerr << ">>>>> TreeMLInfo is defined in evt = " << ievt << std::endl;
         fRunAction->DefineTreeMLInfo();
         fRunAction->GetTreeMLInfo()->Branch("pid", &particleInfo.pid, "pid/I");
         /*should do refunction with TVector3
				 fRunAction->GetTreeMLInfo()->Branch("initialMomentumX", &particleInfo.initialMomentumX, "initialMomentumX/D");
         fRunAction->GetTreeMLInfo()->Branch("initialMomentumY", &particleInfo.initialMomentumY, "initialMomentumY/D");
         fRunAction->GetTreeMLInfo()->Branch("initialMomentumZ", &particleInfo.initialMomentumZ, "initialMomentumZ/D");
         fRunAction->GetTreeMLInfo()->Branch("initialPositionX", &particleInfo.initialPositionX, "initialPositionX/D");
         fRunAction->GetTreeMLInfo()->Branch("initialPositionY", &particleInfo.initialPositionY, "initialPositionY/D");
         fRunAction->GetTreeMLInfo()->Branch("initialPositionZ", &particleInfo.initialPositionZ, "initialPositionZ/D");
					*/
			}else{
         std::cerr << ">>>>> TreeMLInfo was already defined <<<<<\n\n" << std::endl;
      }

      int numTrajectories = trajectoryContainer->entries();
      for(int itrajectory = 0; itrajectory < numTrajectories; ++itrajectory){
	   	WLSTrajectory* trajectory = (WLSTrajectory *)((*(evt->GetTrajectoryContainer()))[itrajectory]);
         G4VTrajectoryPoint* initialPoint = trajectory->GetPoint(0);
         if (trajectory) {
             // PID（粒子の識別子）の取得
             particleInfo.pid = trajectory->GetPDGEncoding();
             // 初期運動量の取得
             G4ThreeVector initialMomentum = trajectory->GetInitialMomentum();
             particleInfo.initialMomentumX = initialMomentum.x();
             particleInfo.initialMomentumY = initialMomentum.y();
             particleInfo.initialMomentumZ = initialMomentum.z();
             // 初期位置の取得
             G4ThreeVector g4initialPosition = initialPoint->GetPosition();
             particleInfo.initialPositionX = g4initialPosition.x();
             particleInfo.initialPositionY = g4initialPosition.y();
             particleInfo.initialPositionZ = g4initialPosition.z();
             // TTreeにデータを詰める
            G4cout << ">>>>> Fill to TreeMLInfo  <<<<<\n\n" << G4endl; 
            fRunAction->GetTreeMLInfo()->Fill();
         }
      }
   }else{
         std::cerr << ">>>>> there is no trajectory." << std::endl;
   }
	#endif

   

   // Save the Random Engine Status at the of event
   if ( fRunAction->GetRndmFreq() == 2 ) {
      G4Random::saveEngineStatus("endOfEvent.rndm");
      G4int evtNb = evt->GetEventID();
      if ( evtNb%fPrintModulo == 0 ) {
         G4cout << "\n---> End of Event: " << evtNb << G4endl;
         G4Random::showEngineStatus();
      }
   }

   // Get Hits from the detector if any
   G4SDManager * SDman = G4SDManager::GetSDMpointer();
   G4String colName = "PhotonDetHitCollection";
   fMPPCCollID = SDman->GetCollectionID(colName);

   G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
   WLSPhotonDetHitsCollection* mppcHC = 0;

   // Get the hit collections
   if ( HCE ) {
      if ( fMPPCCollID>=0 ) mppcHC = (WLSPhotonDetHitsCollection*)( HCE->GetHC(fMPPCCollID) );
   }

   // Get hit information about photons that reached the detector in this event
   if ( mppcHC ) {
	   //G4int n_hit = mppcHC->entries();
   }

	#if 0
   G4AnalysisManager* ana = G4AnalysisManager::Instance();
   ana->FillNtupleDColumn(0, evt->GetEventID());
   ana->FillNtupleDColumn(1, pEne);
   //ana->FillNtupleDColumn(2, a.getX());
   //ana->FillNtupleDColumn(3, a.getY());
   //ana->FillNtupleDColumn(4, a.getZ());
   ana->FillNtupleDColumn(2, fStacking->GetOpticalNPhotons());
   //ana->FillNtupleDColumn(6, fPhotCountX);
   //ana->FillNtupleDColumn(7, fPhotCountY);
   //ana->FillNtupleDColumn(8, fPhotCountZ);
   //ana->AddNtupleRow();
	#endif


   //std::vector <double> Pinitial(ndim,0);
   TVector3 Pinitial;
   TVector3 initialPosition;
	 double mass, kineEnergy, p_abs, pt, edep, stepl, dedx, polar;
   int layer = NTHLAYER;
   int allCh = NTHLAYER * NTHLAYER;
   int PID[nmaxParticles]={0};
	int numOutParticle[NoutCCQE]={0};
   // this is because variable length array can not be initilized by xChPE[allCh]={};
   const int cAllCh = NTHLAYER * NTHLAYER;
   int xChPE[cAllCh]={};
   int yChPE[cAllCh]={};
   int zChPE[cAllCh]={};
   double xChTM[cAllCh]={};
   double yChTM[cAllCh]={};
   double zChTM[cAllCh]={};
	 double fmyMomentumRange = myMomentumRange;

   if ( !fRunAction->GetTreeEvtAct1() ){
      std::cerr << ">>>>> treeEvtAct1 is defined in evt = " << ievt << std::endl;
      fRunAction->DefineTreeEvtAct1();
      fRunAction->GetTreeEvtAct1()->Branch("ievt",  &ievt,  "ievt/I");
      fRunAction->GetTreeEvtAct1()->Branch("PID",  PID,  Form("PID[%d]/I",nmaxParticles));
      fRunAction->GetTreeEvtAct1()->Branch("mass",  &mass,  "mass/D");
      fRunAction->GetTreeEvtAct1()->Branch("kineticEenergy",&kineEnergy,"kineticEnergy/D");
      fRunAction->GetTreeEvtAct1()->Branch("Pinitial",&Pinitial);
      fRunAction->GetTreeEvtAct1()->Branch("initialPosition",&initialPosition);
      fRunAction->GetTreeEvtAct1()->Branch("numOutParticle",numOutParticle,Form("numOutParticle[%d]/I",NoutCCQE));      
      fRunAction->GetTreeEvtAct1()->Branch("myMomentumRange",  &fmyMomentumRange,  "myMomentumRange/D");
      fRunAction->GetTreeEvtAct1()->Branch("edep",  &edep,  "edep/D");
      //fRunAction->GetTreeEvtAct1()->Branch("p_abs", &p_abs, "p_abs/D");
      //fRunAction->GetTreeEvtAct1()->Branch("pt",    &pt,    "pt/D");
      //fRunAction->GetTreeEvtAct1()->Branch("stepl", &stepl, "stepl/D");
      //fRunAction->GetTreeEvtAct1()->Branch("dedx",  &dedx,  "dedx/D");
      //fRunAction->GetTreeEvtAct1()->Branch("layer", &layer, "layer/I");
      //fRunAction->GetTreeEvtAct1()->Branch("allCh", &allCh, "allCh/I"); //  necessary to define below
      //fRunAction->GetTreeEvtAct1()->Branch("xChPE",  xChPE,  "xChPE[allCh]/I");
      //fRunAction->GetTreeEvtAct1()->Branch("yChPE",  yChPE,  "yChPE[allCh]/I");
      //fRunAction->GetTreeEvtAct1()->Branch("zChPE",  zChPE,  "zChPE[allCh]/I");
      //fRunAction->GetTreeEvtAct1()->Branch("xChTM",  xChTM,  "xChTM[allCh]/D");
      //fRunAction->GetTreeEvtAct1()->Branch("yChTM",  yChTM,  "yChTM[allCh]/D");
      //fRunAction->GetTreeEvtAct1()->Branch("zChTM",  zChTM,  "zChTM[allCh]/D");
   } else  {
      std::cerr << ">>>>> treeEvtAct was already defined" << std::endl;
   }

	//for (int i = 0; i < NTHLAYER * NTHLAYER; i++) {
   //	xChPE[i] = fPhotCountX[i];
   //	yChPE[i] = fPhotCountY[i];
   //	zChPE[i] = fPhotCountZ[i];
   //	xChTM[i] = (double) ( fFastestTimeX[i] - fFastestDetResponce );
   //	yChTM[i] = (double) ( fFastestTimeY[i] - fFastestDetResponce );
   //	zChTM[i] = (double) ( fFastestTimeZ[i] - fFastestDetResponce );
   //   if ( xChTM[i] > 999 ) xChTM[i] = 0;
   //   if ( yChTM[i] > 999 ) yChTM[i] = 0;
   //   if ( zChTM[i] > 999 ) zChTM[i] = 0;
	//}


   //if (fVerboseLevel>0)
   // G4GeneralParticleSource: GPS  GetSouce()
   // input in WLSPrimaryGeneratorAction.cc
	
	G4ParticleGun* fParticleGun = fPrimarysource->GetSouce();
	//TTree* treeInitialParticles = fPrimarysource->ftreeInitialParticles;
   G4ParticleDefinition* defPrim = fParticleGun->GetParticleDefinition();   
	
   int Nmuons = fPrimarysource->GetNMuons();numOutParticle[0]=Nmuons;
	int Nproton = fPrimarysource->GetNProtons();numOutParticle[1]=Nproton;
	G4cerr << "Nmuons=" << Nmuons << G4endl;
	G4cerr << "Nproton=" << Nproton << G4endl;
	//G4cerr << "Npion=" << Npion << G4endl;

	int Nparticle = fPrimarysource->GetNParticles();
   std::cout << "Nparticle = " << Nparticle << std::endl;
   for (int iparticle = 0; iparticle < Nparticle; ++iparticle){
      PID[iparticle] = (fPrimarysource->fPID)[iparticle];
   }


   double pMass = defPrim->GetPDGMass();
   double particleKineticEnergy = fParticleGun->GetParticleEnergy();
   //double pMomn = fParticleGun->GetParticleMomentum();
   G4ThreeVector particleMomentumDirection = fParticleGun->GetParticleMomentumDirection ();
   G4ThreeVector g4initialPosition = fParticleGun->GetParticlePosition();
	double polarAngle = particleMomentumDirection.theta();     // get polar angle
	double azimuAngle = particleMomentumDirection.phi();       // get azimuth angle
	double cosTheta   = particleMomentumDirection.cosTheta();  // get cos theta 

   //fStacking->NewStage();

   G4cout << "<<< Event  " << evt->GetEventID() << " ended." << G4endl;
   G4cout << "<<< fPrimary mass             = " << pMass << G4endl; // add
   G4cout << "<<< fPrimary kinetic energy           = " << particleKineticEnergy << G4endl; // add
   G4cout << "<<< fPrimary momentum dir X (unit vector) = " << particleMomentumDirection.getX() << G4endl; // add
   G4cout << "<<< fPrimary momentum dir Y (unit vector) = " << particleMomentumDirection.getY() << G4endl; // add
   G4cout << "<<< fPrimary momentum dir Z (unit vector) = " << particleMomentumDirection.getZ() << G4endl; // add
   G4cout << "<<< fPrimary momentum polar angle   = " << polarAngle << "[rad] " << polarAngle/3.141592*180 << "[deg]" <<G4endl; // add
   G4cout << "<<< fPrimary momentum azimuth angle = " << azimuAngle << "[rad] " << azimuAngle/3.141592*180 << "[deg]" <<G4endl; // add
   G4cout << "<<< fPrimary momentum cos theta     = " << cosTheta   << "[rad] " <<G4endl; // add

   
   double absMom2 = sqrt( pow(particleKineticEnergy+pMass,2) - pow(pMass,2) ); // |p|2 = E2 - m2
  	//Pinitial[0] = absMom2 * particleMomentumDirection.getX(),
		//Pinitial[1] = absMom2 * particleMomentumDirection.getY(),
		//Pinitial[2] = absMom2 * particleMomentumDirection.getZ(); 
		for (int idim = 0; idim < ndim; idim++)
		{
			Pinitial[idim] = absMom2 * particleMomentumDirection(idim);
			initialPosition[idim] = g4initialPosition(idim);
		}
		

   //G4cout << "<<< momentum GetParticleMomentum() = " << pMomn << G4endl;
   G4cout << "<<< fPrimary momentum p   = sqrt( (KinEne + Mass)^2 - Mass^2 ) = " << absMom2 << G4endl;
   G4cout << "<<< fPrimary momentum X   = "<< Pinitial[0] << G4endl;
   G4cout << "<<< fPrimary momentum Y   = "<< Pinitial[1] << G4endl;
   G4cout << "<<< fPrimary momentum Z   = "<< Pinitial[2] << G4endl;

   // transverse P: Pt = P * sin(azimuth)
   //double tanphi  = TMath:: Tan(particleMomentumDirection.getY()/particleMomentumDirection.getX());
   double phi_azim= TMath::ATan(particleMomentumDirection.getY()/particleMomentumDirection.getX());
   //add 23/11/24
   G4cout << "<<< fPrimary Pt      = " << absMom2 * TMath::Sin(phi_azim) << G4endl; //Really...?
   G4cout << "<<< fPrimary posi X  = " << g4initialPosition.getX() << G4endl; // add
   G4cout << "<<< fPrimary posi Y  = " << g4initialPosition.getY() << G4endl; // add
   G4cout << "<<< fPrimary posi Z  = " << g4initialPosition.getZ() << G4endl; // add
   //G4cout << "<<< Ngenerated photon= " << fStacking->GetOpticalNPhotons() << G4endl; // add
   #if 0
   G4cout<< "PhotCountX ";
   for (int i=0; i<NTHLAYER*NTHLAYER; i++) G4cout<< " " << fPhotCountX[i] ;
   G4cout<< G4endl; // 
   G4cout<< "PhotCountY ";
   for (int i=0; i<NTHLAYER*NTHLAYER; i++) G4cout<< " " << fPhotCountY[i] ;
   G4cout<< G4endl; // 
   G4cout<< "PhotCountZ ";
   for (int i=0; i<NTHLAYER*NTHLAYER; i++) G4cout<< " " << fPhotCountZ[i] ;
   G4cout<< G4endl; // 
   #endif

	mass   = pMass;
	kineEnergy = particleKineticEnergy;
	p_abs      = absMom2;
	pt     = absMom2 * TMath::Sin(phi_azim);
   edep   = fEdep;
   stepl  = fStepLength;
   dedx   = fEdep/fStepLength;
   //nGenePEs= 0;//fStacking->GetOpticalNPhotons();

   G4cout << "<<< energy ioize dep  = " << edep << G4endl; // add
   G4cout << "<<< energy delta dep  = " << fRunAction->GetChargedSecondary() << G4endl; // add
   G4cout << "<<< energy brems dep  = " << fRunAction->GetNeutralSecondary() << G4endl; // add
   G4cout << "<<< step legth  = " << stepl<< G4endl; // add
   G4cout << "<<< dedx        = " << dedx << G4endl; // add
 
   fRunAction->GetTreeEvtAct1()->Fill();


	//
	//
	//
	const int chall_dum1 = 999;
	const int chall_dum2 = 999;
   int lyall0 = nlayers_dummy;
   int chall1 = chall_dum1;
   int chall2 = chall_dum2;
   //float chx[chall_dum1]={};
   //float chy[chall_dum1]={};
   float chw[nlayers_dummy][nfibers_dummy];//={};
   //memset(chw, 0, nlayers_dummy*chall_dum2*sizeof(float));

   if ( !fRunAction->GetTreeEvtAct2() ) {
      std::cerr << ">>>>> treeEvtAct2 is defined in evt = " << ievt << std::endl;
      fRunAction->DefineTreeEvtAct2();
      fRunAction->GetTreeEvtAct2()->Branch("ievt",  &ievt,  "ievt/I");
      fRunAction->GetTreeEvtAct2()->Branch("polar", &polarAngle,  "polar/D");
      //fRunAction->GetTreeEvtAct2()->Branch("nlayers_dummy", &nlayers_dummy, "nlayers_dummy/I"); //  necessary to define below
      //fRunAction->GetTreeEvtAct2()->Branch("nfibers_dummy", &nfibers_dummy, "nfibers_dummy/I"); //  necessary to define below
      //fRunAction->GetTreeEvtAct2()->Branch("chall1", &chall1, "chall1/I"); //  necessary to define below
      //fRunAction->GetTreeEvtAct2()->Branch("chall2", &chall2, "chall2/I"); //  necessary to define below
      //fRunAction->GetTreeEvtAct2()->Branch("chx",  chx,  "chx[chall1]/F");
      //fRunAction->GetTreeEvtAct2()->Branch("chy",  chy,  "chy[chall1]/F");
      //fRunAction->GetTreeEvtAct2()->Branch("chw",  chw,  "chw[nlayers_dummy][nfibers_dummy]/F");
      fRunAction->GetTreeEvtAct2()->Branch("chw",  chw,  Form("chw[%d][%d]/F",nlayers_dummy,nfibers_dummy)); // this works
	}

	polar = polarAngle;
   for (int m = 0; m < nlayers_dummy; m++) {
	   for (int i = 0; i < nfibers_dummy; i++) {
			//chx[i] = fEdepoLayerScintiX[m][i];
			//chy[i] = fEdepoLayerScintiY[m][i];
			chw[m][i] = fEdepoLayerScinti[m][i];
         //if ( chw[m][i]>0 ) G4cerr << "fiber layer = " << m << " fiber id = " << i << " edep = " << chw[m][i] << G4endl;
		}
	}

   G4cout << ">>>>> Fill to tree  <<<<<\n\n" << G4endl; 
   fRunAction->GetTreeEvtAct2()->Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int WLSEventAction::GetEventNo()
{
  return fpEventManager->GetConstCurrentEvent()->GetEventID();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSEventAction::SetEventVerbose(G4int level)
{
  fVerboseLevel = level;
}
