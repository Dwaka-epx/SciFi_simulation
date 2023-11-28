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
// $Id: WLSPrimaryGeneratorAction.cc 89119 2015-03-20 09:09:45Z gcosmo $
//
/// \file optical/wls/src/WLSPrimaryGeneratorAction.cc
/// \brief Implementation of the WLSPrimaryGeneratorAction class


#include "G4ios.hh"
#include "G4Event.hh"

#include "G4GeneralParticleSource.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"

#include "Randomize.hh"

#include "WLSPrimaryGeneratorAction.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSPrimaryGeneratorMessenger.hh"

#include "G4SystemOfUnits.hh"

#include "G4AutoLock.hh"


#include <iostream>
#include <fstream>

namespace {
  G4Mutex gen_mutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WLSPrimaryGeneratorAction::fFirst = false;

#if 0
WLSPrimaryGeneratorAction:: WLSPrimaryGeneratorAction(WLSDetectorConstruction* dc)
//WLSPrimaryGeneratorAction:: WLSPrimaryGeneratorAction(WLSDetectorConstruction* dc, WLSEventAction* eventAction)
//  : fEventAction(eventAction) // add
#endif
#if 1
WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction(WLSDetectorConstruction* dc, G4String NeutName, int NeutSttNum, int NeutGenNum,
	WLSRunAction* runaction)
	: _NeutTree(0),
	  _NeutName(NeutName),     // input neut file .root
	  _NeutSttNum(NeutSttNum), // NEUT starting generate event number 
	  _NeutGenNum(NeutGenNum), // Geant generate EventID 
	  _neutMode(0),
	  _neutNSeeds(0),
	  fRunAction(runaction)
#endif
{
  fDetector = dc;
  fIntegralTable = NULL;

  //fParticleGun = new G4GeneralParticleSource();
  fParticleGun = new G4ParticleGun();
  fGunMessenger = new WLSPrimaryGeneratorMessenger(this);

	//fEventAction = new WLSEventAction();

  // G4String particleName;
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  fTimeConstant = 0.;

//  fParticleGun->SetParticleDefinition(particleTable->
//                               FindParticle(particleName="opticalphoton"));

#if 1 
	//
	// preparation NEUT info.
	//
	TString xxx = "[WLSPrimaryGeneratorAction::WLSPrimaryGeneratorAction]"; 
	//TString inFileName = "/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1_build_cos7/run_neuInteraction/neut_5.4.0_600MeV_C.card.vect.root";
	//TString inFileName = _NeutName;
	G4String inFileName = _NeutName;
	std::cout << "\n\n" << xxx << ": input NEUT file\n >>>>> " << inFileName << std::endl;
	TFile *fin = new TFile(inFileName, "READ");
	fin->GetObject("neuttree", _NeutTree);
	_NeutTree->SetBranchAddress("vectorbranch", &_NeutVec);

	std::cout << xxx << ": Reading tree with " << _NeutTree->GetEntries() << " NEUT events." << std::endl;
	std::cout << xxx << ": generation start event number : " << _NeutSttNum
				 << "\n\n" << std::endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPrimaryGeneratorAction::~WLSPrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fGunMessenger;
	if (fIntegralTable) {
	  fIntegralTable->clearAndDestroy();
	  delete fIntegralTable;
	}
}



//const G4GeneralParticleSource* WLSPrimaryGeneratorAction::GetSouce() 
G4ParticleGun* WLSPrimaryGeneratorAction::GetSouce()
{
  	return fParticleGun;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetDecayTimeConstant(G4double time)
{
  fTimeConstant = time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::BuildEmissionSpectrum()
{
	if (fIntegralTable) return;

	const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

	G4int numOfMaterials = G4Material::GetNumberOfMaterials();

	if(!fIntegralTable)fIntegralTable = new G4PhysicsTable(numOfMaterials);

	for (G4int i=0 ; i < numOfMaterials; i++) {

		 G4PhysicsOrderedFreeVector* aPhysicsOrderedFreeVector =
														 new G4PhysicsOrderedFreeVector();

		 G4Material* aMaterial = (*theMaterialTable)[i];

		 G4MaterialPropertiesTable* aMaterialPropertiesTable =
												aMaterial->GetMaterialPropertiesTable();

		 if (aMaterialPropertiesTable) {
			 G4MaterialPropertyVector* theWLSVector =
							 aMaterialPropertiesTable->GetProperty("WLSCOMPONENT");

			 if (theWLSVector) {
				 G4double currentIN = (*theWLSVector)[0];
				 if (currentIN >= 0.0) {
					 G4double currentPM = theWLSVector->Energy(0);
					 G4double currentCII = 0.0;
					 aPhysicsOrderedFreeVector->
												 InsertValues(currentPM , currentCII);
					 G4double prevPM  = currentPM;
					 G4double prevCII = currentCII;
					 G4double prevIN  = currentIN;

					 for (size_t j = 1;
							j < theWLSVector->GetVectorLength();
							j++)
					 {
						currentPM = theWLSVector->Energy(j);
						currentIN = (*theWLSVector)[j];
						currentCII = 0.5 * (prevIN + currentIN);
						currentCII = prevCII + (currentPM - prevPM) * currentCII;
						aPhysicsOrderedFreeVector->
												 InsertValues(currentPM, currentCII);
						prevPM  = currentPM;
						prevCII = currentCII;
						prevIN  = currentIN;
					 }
				 }
			 }
		 }
		 fIntegralTable->insertAt(i,aPhysicsOrderedFreeVector);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
#if 1 

	fevtID = anEvent->GetEventID();
	//TFile* foutFile = fRunAction->GetTFile(); std::cerr << "	TFile* foutFile = fRunAction->GetTFile();" << std::endl;
	//foutFile->cd(); std::cerr << "	foutFile->cd();" << std::endl;
	if ( !fRunAction->GetTreeInitialParticles() ){
		std::cerr << ">>>>> treeInitialParticles is defined in evt = " << fevtID << std::endl;
		fRunAction->DefineTreeInitialParticles();
		std::cerr << fRunAction << std::endl;
     	fRunAction->GetTreeInitialParticles()->Branch("fevtID", &fevtID, "fevtID/I");
	
	}else{
		std::cerr << ">>>>> treeInitialParticles was already defined " << std::endl;
		std::cerr << fRunAction << std::endl;
	}
	fRunAction->GetTreeInitialParticles()->Fill(); std::cerr << "fRunAction->GetTreeInitialParticles()->Fill();" << std::endl;
	//fRunAction->GetTreeInitialParticles()->Scan();
	
	G4cout << "\n\n>>>>> CALLED GeneratePrimaries at the begining of event" << G4endl;
	G4cout << ">>>>> [Geant] Geant generate EventID : " << fevtID << " / " << _NeutGenNum << G4endl;
	G4cout << ">>>>> [NEUT]: NEUT generate event number : " << _NeutSttNum << G4endl;

	//std::ofstream ofs("./sim_output/initial_particles.txt", std::ios::in | std::ios::ate);

	//for (int evt = 0; evt < nevs; ++evt) {
		//_NeutTree->GetEntry( evt );
		//_NeutTree->GetEntry( fevtID );
		_NeutTree->GetEntry( _NeutSttNum );

		int npart = _NeutVec->Npart();
		int mode  = _NeutVec->Mode;
		int nGenerate=0;
		_nParticles=0;
		_nNeutrons=0; // neutron
		_nProtons=0; // proton
		_nPions=0; // pion
		_nMuons=0; // muon

		fPID[nmaxParticles]={0};
		frndx[nmaxParticles][ndim]={0};
		fvertex[nmaxParticles][ndim]={0};
		fInitialMomentum[nmaxParticles][ndim]={0};
		fnumOutParticle[2]={0}; // [Nmuon,Nproton]
		fEnergy[nmaxParticles]={0};


		//if ( evt%1==0 ) 
		std::cout<< "\033[31m>>>>> [NEUT] generate Ev: " << _NeutSttNum
					<< ", Neut mode " << mode << " saved N particles " << npart << " <<<<<\n \033[m" << std::endl;

		for (int pit = 0; pit < npart; ++pit) { // n particles loop
			NeutPart const & part = *_NeutVec->PartInfo(pit);
			if ( part.fIsAlive==1 ) {
				_nParticles++;
				TLorentzVector fP     = part.fP;      // 4 momentum of the particle (MeV/c, MeV)
				TLorentzVector fPosIni= part.fPosIni; // Initilal position in the nucleus
				TLorentzVector fPosFin= part.fPosFin; // Final(current) position in the nucleus
				fEnergy[pit] = fP.E();
				double Pabs = fP.Vect().Mag();
				double Mass = sqrt(fEnergy[pit]*fEnergy[pit]-Pabs*Pabs); // E^2 = M^2+P^2
				if ( part.fPID==12 || part.fPID==14 ) _neutEnergy  = fEnergy[pit]; // take neutrino initial energy
				if ( part.fPID==12 || part.fPID==14 ) _neutMomentum= Pabs; // take neutrino initial momentum
				if ( part.fPID==2112 ) _nNeutrons++ ; // neutron
				if ( part.fPID==2212 ) _nProtons++ ; // proton
				if ( part.fPID==211  ) _nPions++ ; // pion
				if ( part.fPID==13  ) _nMuons++ ; // pion
				fPID[pit] = part.fPID;
				G4ParticleDefinition* particleDefinition
					= G4ParticleTable::GetParticleTable()->FindParticle(part.fPID);
				double G4Mass = particleDefinition->GetPDGMass();
				//G4cout << "set g4 particle: " << particleDefinition->GetParticleName() << G4endl;
				fParticleGun->SetParticleDefinition(particleDefinition);
				// incident/initial position
				//const int ndim=3;
				//double rndx[ndim];
				//double vertex[ndim]={-30.,-30.,60.};
				//double width = 60.;
				//for (int idim = 0; idim < ndim; ++idim){
				//	frndx[pit][idim] = G4UniformRand()*width;
				//	fInitialMomentum[pit][idim] = fP[idim];
				//}
				
				//fvertex[pit][0] = frndx[pit][0] - 30.;
				//fvertex[pit][1] = frndx[pit][1] - 30.;
				//fvertex[pit][2] = frndx[pit][2] + 70.;
				
	 			double rndx = G4UniformRand();
	 			double rndy = G4UniformRand();
	 			double rndz = G4UniformRand();
	 			double vertex = (rndx*60.-30.); //0.;	
	 			double vertey = (rndy*60.-30.); //0.; 	
	 			double vertez = (rndz*60.+70.); //100.; 	

	 			//fParticleGun->SetParticlePosition(G4ThreeVector(fvertex[pit][0]*mm,fvertex[pit][1]*mm,fvertex[pit][2]*mm)); //G4UniformRand()
	 			fParticleGun->SetParticlePosition(G4ThreeVector(vertex*mm,vertey*mm,vertez*mm)); //G4UniformRand()
				fParticleGun->SetParticleMomentum(G4ThreeVector(fP.Px(),fP.Py(),fP.Pz()));

				// Energy you set must be "kinetic energy"
				// E=T+M
				// P=sqrt{T(T+2M)}
				// T=sqrt(P^2+M^2)-M の"T"。
				//fParticleGun->SetParticleEnergy(0.511*MeV); // kinetic 
				//double Ekin = sqrt(Pabs*Pabs+G4Mass*G4Mass) - G4Mass;
				//fParticleGun->SetParticleEnergy(Ekin*MeV); // kinetic

				G4ParticleMomentum fParticleMomDirection = fParticleGun->GetParticleMomentumDirection();
				//TVector3 fParticleMomDirection = fParticleGun->GetParticleMomentumDirection();
				std::cerr << " >>>>> NeutPart info. pid = " << std::right << std::setw(4) << part.fPID 
							 << " IsAlive = " << std::right << std::setw(2) << part.fIsAlive
							 << std::endl;
				std::cerr
					  << " >>>>> NeutPart info."
					  //<< " >>>>> NeutPart info. pid = " << std::right << std::setw(4) << part.fPID
					  //<< " Status = " << std::right << std::setw(2) << part.fStatus
					  << std::setprecision(5)
					  << "  M "  << std::right << std::setw(6) << Mass
					  << "  P "  << std::right << std::setw(6) << Pabs
					  << "  Energy "  << std::right << std::setw(6) << fEnergy
					  << "  Pex direction"  << std::right << std::setw(6) << fP.Px()/Pabs
					  << "  Pey direction"  << std::right << std::setw(6) << fP.Py()/Pabs
					  << "  Pez direction"  << std::right << std::setw(6) << fP.Pz()/Pabs
					  //<< "  Xfx " << right << setw(9) << fPosFin.Px() 
					  //<< "  Xfy " << right << setw(9) << fPosFin.Py() 
					  //<< "  Xfz " << right << setw(9) << fPosFin.Pz()
				<< std::endl;
				double G4MomAbs = fParticleGun->GetParticleMomentum();
				std::cerr
					  << " >>>>> Set G4: " << particleDefinition->GetParticleName()
					  << "  M: " << std::right << std::setw(5) << G4Mass
					  << "  P: " << std::right << std::setw(5) << G4MomAbs
					  << "  E(sqrt(M^2+p^2)): "   << std::right << std::setw(5) << sqrt(G4Mass*G4Mass+G4MomAbs*G4MomAbs)
					  << "  Ekin(GetPartEne()): " << std::right << std::setw(5) << fParticleGun->GetParticleEnergy()
					  << std::setprecision(4)
					  << "  Pex: " << std::right << std::setw(3) << fParticleMomDirection.x()
					  << "  Pey: " << std::right << std::setw(3) << fParticleMomDirection.y()
					  << "  Pez: " << std::right << std::setw(3) << fParticleMomDirection.z()
					  << std::endl;
					std::cerr << "\033[31m >>>>> IS REGISTERED <<<<<\n \033[m" << std::endl;
					fParticleGun->GeneratePrimaryVertex(anEvent);
					nGenerate++;

					//Fill
			} else {
				std::cerr << std::endl;
			}
		}
	//fnumOutParticle[0] = _nMuons;
	//fnumOutParticle[1] = _nProtons;

	SetNeutrinoMode(mode);
	SetNeutrinoNSeed(nGenerate); //???
	_NeutSttNum++;
	return;
#endif

	G4cout << " GeneratePrimaries is called at the begining of event" << G4endl;
	if (!fFirst) {
		fFirst = true;
	  	BuildEmissionSpectrum();	
	}

#ifdef use_sampledEnergy
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4double sampledEnergy = 3*eV;

  for (size_t j=0 ; j<theMaterialTable->size() ; j++) {
		G4Material* fMaterial = (*theMaterialTable)[j];
		if (fMaterial->GetName() == "PMMA" ) {
			G4MaterialPropertiesTable* aMaterialPropertiesTable =
												  fMaterial->GetMaterialPropertiesTable();
			const G4MaterialPropertyVector* WLSIntensity =
						 aMaterialPropertiesTable->GetProperty("WLSCOMPONENT");

			if (WLSIntensity) {
				G4int MaterialIndex = fMaterial->GetIndex();
				G4PhysicsOrderedFreeVector* WLSIntegral =
				  (G4PhysicsOrderedFreeVector*)((*fIntegralTable)(MaterialIndex));

				G4double CIImax = WLSIntegral->GetMaxValue();
				G4double CIIvalue = G4UniformRand()*CIImax;

				sampledEnergy = WLSIntegral->GetEnergy(CIIvalue);
			}
		}
  }

  //fParticleGun->SetParticleEnergy(sampledEnergy);
	
  //particleGun->SetParticleEnergy( 1.*GeV );
  //で設定する値は運動量や全エネルギーではなく"運動エネルギー"となる。
  // E = T + M
  // P = sqrt{T(T+2M)}
  // T = sqrt(P^2+M^2) - M
  // の"T"。
#endif

  //The code behing this line is not thread safe because polarization
  //and time are randomly selected and GPS properties are global
  G4AutoLock l(&gen_mutex);
  if(fParticleGun->GetParticleDefinition()->GetParticleName()=="opticalphoton"){
	 SetOptPhotonPolar();
	 SetOptPhotonTime();
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonPolar()
{
  G4double angle = G4UniformRand() * 360.0*deg;
  SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
  if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
  {
	  G4cout << "-> warning from WLSPrimaryGeneratorAction::SetOptPhotonPolar()"
				<< ":  the ParticleGun is not an opticalphoton" << G4endl;
	  return;
  }

  G4ThreeVector normal (1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  G4double modul2       = product*product;

  G4ThreeVector e_perpend (0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
  G4ThreeVector e_paralle    = e_perpend.cross(kphoton);

  G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPrimaryGeneratorAction::SetOptPhotonTime()
{
	G4double time = -std::log(G4UniformRand())*fTimeConstant;
	fParticleGun->SetParticleTime(time);
}
