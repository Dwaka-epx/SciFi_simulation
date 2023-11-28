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
// $Id: WLSRunAction.cc 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/src/WLSRunAction.cc
/// \brief Implementation of the WLSRunAction class


#include "WLSRunAction.hh"
#include "WLSRunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"

#include "WLSDetectorConstruction.hh"
#include "WLSSteppingAction.hh"
//#include "WLSPrimaryGeneratorAction.hh"

#include <ctime>
#include <sstream>
#include <string.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//WLSRunAction::WLSRunAction(G4String name)
WLSRunAction::WLSRunAction(WLSDetectorConstruction* det,G4String name)
  : fSaveRndm(0), 
	fAutoSeed(false), 
    fDetConstruction(det),
	fName(name),
	fFile(0),// initialize
	fTreeEvtAct1(0), // initialize
	fTreeEvtAct2(0), // initialize
	fTreeInitialParticles(0), // initialize
	fTreeStkAct(0), // initialize
	fTreeStpAct(0) // initialize
{
  fRunMessenger = new WLSRunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSRunAction::~WLSRunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  	G4cout << "mytest: Run " << aRun->GetRunID() << " start." << G4endl;
  	G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  	G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");

	#if 0
  	G4AnalysisManager* ana = G4AnalysisManager::Instance();
  	G4cout << "mytest: G4AnalysisManager : " << ana->GetType() << G4endl;
	// ----- open an output file
	ana->OpenFile(fName);
	// ----- create ntuple
	//ana->SetVerboseLevel(1);
	ana->CreateNtuple("cube", "cube info");
	ana->CreateNtupleDColumn("evt"); 
	ana->CreateNtupleDColumn("pene"); 
	ana->CreateNtupleDColumn("nphotos"); 	
	// it doesn't work someting like 
	//tree->Branch("nhits", &nhits,"nhits/I");
   //tree->Branch("order", order, "order[nhits]/D");
	#endif

	std::string str = fName;//+".root";
	std::cout << "fFile =" << fName << std::endl;
		fFile = new TFile(str.c_str(),"RECREATE"); std::cerr << "fFile = new TFile(str.c_str(),RECREATE);" << std::endl;
		//TTree* fTreeInitialParticles = new TTree("dummy_treeInitialParticle", "");
	#if 0
   tree = new TTree("tree", "");
   int evt; // #of event
   int layer=NTHLAYER;
	int allCh=NTHLAYER*NTHLAYER;
   double xChPE[999];
   double yChPE[999];
   double zChPE[999];
   tree->Branch("evt",   &evt,  "evt/I");
   tree->Branch("layer", &layer,"layer/I");
   tree->Branch("allCh", &allCh,"allCh/I");
   tree->Branch("xChPE",  xChPE,  "xChPE[allCh]/D");
   tree->Branch("yChPE",  yChPE,  "yChPE[allCh]/D");
   tree->Branch("zChPE",  zChPE,  "zChPE[allCh]/D");
	#endif
	
	
#if 0 // 
	for (int i=0; i<NTHLAYER; i++) {
		for (int j=0; j<NTHLAYER; j++) {
			std::stringstream name;
			name << "nPhotXch" << i*NTHLAYER +j << std::ends;
			ana->CreateNtupleDColumn(name.str().data()); 
		}
	}
	for (int i=0; i<NTHLAYER; i++) {
		for (int j=0; j<NTHLAYER; j++) {
			std::stringstream name;
			name << "nPhotiYch" << i*NTHLAYER +j << std::ends;
			ana->CreateNtupleDColumn(name.str().data()); 
		}
	}
   for (int i=0; i<NTHLAYER; i++) {
      for (int j=0; j<NTHLAYER; j++) {
			std::stringstream name;
			name << "nPhotZch" << i*NTHLAYER +j << std::ends;
			ana->CreateNtupleDColumn(name.str().data()); 
		}
	}
#endif
	//ana->FinishNtuple(0);

  	if ( fAutoSeed ) {
     // automatic (time-based) random seeds for each run
     G4cout << "*******************" << G4endl;
     G4cout << "*** AUTOSEED ON ***" << G4endl;
     G4cout << "*******************" << G4endl;
     long seeds[2];
     time_t systime = time(NULL);
     seeds[0] = (long) systime;
     seeds[1] = (long) (systime*G4UniformRand());
     G4Random::setTheSeeds(seeds);
     G4Random::showEngineStatus();
  	} else {
     G4Random::showEngineStatus();
  	}

  	if (fSaveRndm > 0) {
		G4Random::saveEngineStatus("BeginOfRun.rndm");
	}

	// ----- output text file
	//outputf.open(fName+".txt", std::ios::app);
}


void WLSRunAction::EndOfRunAction(const G4Run* )
{
	G4cerr << "\n>>>>> CALLED: WLSRunAction::EndOfRunAction" << G4endl;
	
	//fTreeInitialParticles->Print();

   //int nbEvents = aRun->GetNumberOfEvent();
   //double DepositedBySecondaries = fEnergyCharged/nbEvents + fEnergyNeutral/nbEvents;
   //double TotalDepositByPrimary  = analysisManager->GetH1(1)->mean() + DepositedBySecondaries;

  	if (fSaveRndm == 1) {
     G4Random::showEngineStatus();
     G4Random::saveEngineStatus("endOfRun.rndm");
  	}

   TTree *tree2 = new TTree("tree_yread", "readout position vector");
   int id; TVector3 *v3 = 0;
   tree2->Branch("v3", &v3);
   tree2->Branch("id", &id);

	std::vector<int>    sciID;
   std::vector<TVector3> vec;
	sciID = fDetConstruction->GetReadoutYID();
   vec   = fDetConstruction->GetReadoutY();
   G4cerr << ">>>>> ReadoutYID.size = " << sciID.size() << G4endl;
   G4cerr << ">>>>> ReadoutY .size = " << vec.size()   << G4endl;
   G4cerr << ">>>>> Remmber position is defined whintin the mother volume" << G4endl;
   for (unsigned int i = 0; i<vec.size(); i++) {
      v3 = &vec[i];
		id = sciID[i];	
      //G4cerr << "readout point no." << i << " fiber id = " << id
      //       << " z = " << v3->Z() << " y = " << v3->Y() << " x = " << v3->X()
      //       << G4endl;
      tree2->Fill();
   }
	
   TTree *tree3 = new TTree("tree_fread", "readout position vector");
   tree3->Branch("v3", &v3);
   tree3->Branch("id", &id);

	sciID = fDetConstruction->GetReadoutFID();
   vec   = fDetConstruction->GetReadoutF();
   G4cerr << ">>>>> ReadoutFID.size = " << sciID.size() << G4endl;
   G4cerr << ">>>>> ReadoutF .size = " << vec.size()   << G4endl;
   G4cerr << ">>>>> Remmber position is defined whintin the mother volume" << G4endl;
   for (unsigned int i = 0; i<vec.size(); i++) {
      v3 = &vec[i];
      id = sciID[i];
      //G4cerr << "readout point no." << i << " fiber id = " << id
      //       << " z = " << v3->Z() << " y = " << v3->Y() << " x = " << v3->X()
      //       << G4endl;
      tree3->Fill();
   }

	#if 0
	G4AnalysisManager* ana = G4AnalysisManager::Instance(); 
	ana->Write();
	ana->CloseFile();
	#endif
   std::cerr << "write tree in file" << std::endl;
		fTreeEvtAct1->Write();std::cerr << "fTreeEvtAct1->Write();" << std::endl;
		fTreeEvtAct2->Write();std::cerr << "fTreeEvtAct2->Write();" << std::endl;
	//if(fTreeMLInfo)fTreeMLInfo->Write();
	//fTreeStkAct->Write();
		fTreeStpAct->Write();			std::cerr << "fTreeStpAct->Write();" << std::endl;
		tree2->Write();			std::cerr << "tree2->Write();" << std::endl;
		tree3->Write();			std::cerr << "tree3->Write();" << std::endl;
		if(fTreeInitialParticles){
			fTreeInitialParticles->Write();
			std::cerr << "fTreeInitialParticles->Write();" << std::endl; 
			std::cerr << this << std::endl;
		}
	std::cerr << "close file" << std::endl;
	fFile->Close();
}
