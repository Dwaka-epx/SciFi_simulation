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
// $Id: WLSPhysicsList.cc 78066 2013-12-03 11:08:36Z gcosmo $
//
/// \file optical/wls/src/WLSPhysicsList.cc
/// \brief Implementation of the WLSPhysicsList class
//
//
#include "WLSPhysicsList.hh"
#include "WLSPhysicsListMessenger.hh"

#include "WLSExtraPhysics.hh"
#include "WLSOpticalPhysics.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#include "G4PhysListFactory.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"


#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "WLSStepMax.hh"

#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"

#include "G4Scintillation.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"
#include "G4StepLimiterPhysics.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsList::WLSPhysicsList(G4String physName) 
	: G4VModularPhysicsList()
{
   G4LossTableManager::Instance();

   defaultCutValue  = 1. * mm;
	fCutForGamma     = defaultCutValue;
   fCutForElectron  = defaultCutValue;
   fCutForPositron  = defaultCutValue;
   fCutForProton    = defaultCutValue;

	//    G4PhysListFactory factory;
   G4VModularPhysicsList* phys = NULL;
   if        (physName=="QGSP_BERT_HP") {
       phys = new QGSP_BERT_HP;
   } else if (physName=="FTFP_BERT") {
       phys = new FTFP_BERT;
   } else if (physName=="G4EmStandardPhysics_option2") {
   	phys = new G4VModularPhysicsList();
	   //phys->RegisterPhysics(new G4EmStandardPhysics());
   	phys->RegisterPhysics(new G4EmStandardPhysics_option2());
      //phys->RegisterPhysics(new G4StepLimiterPhysics()); //追加。 
	}

//    if (factory.IsReferencePhysList(physName)) {
//       phys = factory.GetReferencePhysList(physName);
//       if(!phys)G4Exception("WLSPhysicsList::WLSPhysicsList","InvalidSetup",
//                            FatalException,"PhysicsList does not exist");
       fMessenger = new WLSPhysicsListMessenger(this);
//    }

	for (int i = 0; ; ++i) {
   	G4VPhysicsConstructor* elem 
      	= const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
      if (elem == NULL) break;
      G4cerr << "##### RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
      RegisterPhysics(elem);
 	}

   fAbsorptionOn = true;
    
   //This looks complex, but it is not:
   //Get from base-class the pointer of the phsyicsVector
   //to be used. Remember: G4VModularPhysicsList is now a split class.
   //Why G4VModularPhysicsList::RegisterPhysics method is not used instead?
   //If possible we can remove this...
   fPhysicsVector = GetSubInstanceManager().offset[GetInstanceID()].physicsVector;
    
   fPhysicsVector->push_back(new WLSExtraPhysics()); // call extra physcis process @ 21/01/28
   fPhysicsVector->push_back(new G4RadioactiveDecayPhysics());

	// call constrcucter for optical process
	// to hundle scintillation, cherenkov...
   G4cerr << "\n\n\n##### WLSOpticalPhysics is On/Off #####\n\n\n" << G4endl;
   //fPhysicsVector->push_back(fOpticalPhysics = new WLSOpticalPhysics(fAbsorptionOn));

   fStepMaxProcess = new WLSStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhysicsList::~WLSPhysicsList()
{
    delete fMessenger;

    delete fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::ClearPhysics()
{
    for (G4PhysConstVector::iterator p  = fPhysicsVector->begin();
                                     p != fPhysicsVector->end(); ++p) {
        delete (*p);
    }
    fPhysicsVector->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

void WLSPhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();
#if 1
    G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
    MuonPlusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu+",0.986));
    MuonPlusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
    //G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

    G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
    MuonMinusDecayTable -> Insert(new G4MuonDecayChannelWithSpin("mu-",0.986));
    MuonMinusDecayTable -> Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
    //G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);
#endif
  	//G4MuonPlus::MuonPlusDefinition();
  	//G4MuonMinus::MuonMinusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

	// out together with "new WLSOpticalPhysics"
	//SetOpticalPhysicsVerbose(0);
#if 1 
    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* pManager;
    pManager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    pManager = G4MuonMinus::MuonMinus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    pManager = G4PionPlus::PionPlus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(poldecay, idxPostStep);
      pManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    pManager = G4PionMinus::PionMinus()->GetProcessManager();

    if (pManager) {
      if (decay) pManager->RemoveProcess(decay);
      pManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pManager ->SetProcessOrdering(poldecay, idxPostStep);
      pManager ->SetProcessOrdering(poldecay, idxAtRest);
    }
#endif
    AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::RemoveFromPhysicsList(const G4String& name)
{
    G4bool success = false;
    for (G4PhysConstVector::iterator p  = fPhysicsVector->begin();
                                     p != fPhysicsVector->end(); ++p) {
        G4VPhysicsConstructor* e = (*p);
        if (e->GetPhysicsName() == name) {
           fPhysicsVector->erase(p);
           success = true;
           break;
        }
    }
    if (!success) {
       G4ExceptionDescription message;
       message << "PhysicsList::RemoveFromEMPhysicsList "<< name << "not found";
       G4Exception("example WLSPhysicsList::RemoveFromPhysicsList()",
       "ExamWLSPhysicsList01",FatalException,message);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetAbsorption(G4bool toggle)
{
       fAbsorptionOn = toggle;
       RemoveFromPhysicsList("Optical");
       fPhysicsVector->
          push_back(fOpticalPhysics = new WLSOpticalPhysics(toggle));
       fOpticalPhysics->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetCuts()
{
    if (verboseLevel >0) {
        G4cout << "WLSPhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length")
               << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(fCutForGamma, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");
    SetParticleCuts(fCutForProton, G4Proton::Proton()); // p

   //if (verboseLevel>0) DumpCutValuesTable();
   DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetCutForGamma(G4double cut)
{
    fCutForGamma = cut;
    SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

void WLSPhysicsList::SetCutForElectron(G4double cut)
{
    fCutForElectron = cut;
    SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

void WLSPhysicsList::SetCutForPositron(G4double cut)
{
    fCutForPositron = cut;
    SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetStepMax(G4double step)
{
  fStepMaxProcess->SetStepMax(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSStepMax* WLSPhysicsList::GetStepMaxProcess()
{
  return fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (fStepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
         if (pmanager) pmanager ->AddDiscreteProcess(fStepMaxProcess);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetNbOfPhotonsCerenkov(G4int maxNumber)
{
   fOpticalPhysics->SetNbOfPhotonsCerenkov(maxNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhysicsList::SetOpticalPhysicsVerbose(G4int verbose)
{
   fOpticalPhysics->GetCerenkovProcess()          ->SetVerboseLevel(verbose);
   fOpticalPhysics->GetScintillationProcess()     ->SetVerboseLevel(verbose);
   fOpticalPhysics->GetAbsorptionProcess()        ->SetVerboseLevel(verbose);
   fOpticalPhysics->GetRayleighScatteringProcess()->SetVerboseLevel(verbose);
   fOpticalPhysics->GetMieHGScatteringProcess()   ->SetVerboseLevel(verbose);
   fOpticalPhysics->GetBoundaryProcess()          ->SetVerboseLevel(verbose);
}
