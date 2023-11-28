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
// $Id: wls.cc 78066 2013-12-03 11:08:36Z gcosmo $
//
/// \file optical/wls/wls.cc
/// \brief Main program of the optical/wls example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef WIN32
#include <unistd.h>
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "WLSPhysicsList.hh"
#include "WLSDetectorConstruction.hh"
#include "WLSActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// argc holds the number of arguments (including the name) on the command line
// -> it is ONE when only the name is  given !!!
// argv[0] is always the name of the program
// argv[1] points to the first argument, and so on

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) 
{
#ifdef G4MULTITHREADED
  	G4MTRunManager * runManager = new G4MTRunManager;

#else
  	int seed = 123;
  	//if (argc >2) seed = atoi(argv[argc-1]);

   int NeutSttNum = 0; // Neut starting event number
   int NeutGenNum = 0; // Neut # of event generation 
   //G4String NeutName = "/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1_run_cos7/";
   //G4String NeutName = "/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1_run_cos7_210401/";
   //G4String NeutName = "/home/t2k/ogawat/myt2kwork/neut_generator/neut_5.4.0.1_run_cos7_200930/";
   G4String NeutName = "/home/fiberstudy/SciFi_simulation/build/input_neut/";

   //NeutName+= "run_neutrino_test/neut_5.4.0_600MeV_C.card.vect.root";
   //NeutName+= "run_num_600MeV_CC1PI0_ID12/neut_5.4.0_600MeV_C.card.vect.root";
   G4String condition_Name = "neut_5.4.0_675MeV_H2O_numu_1e5event";
   NeutName+= condition_Name+".card.vect.root";
   G4String OutName  = "./sim_output/mydata_"+condition_Name+".root";

	if (argc==3) {        // ./wavy run.mac {OutName}.root
      OutName = argv[2];
	} else if (argc==4) { // ./wavy run.mac *.root random_seed
		OutName = argv[2];
		seed    = atoi(argv[3]);
   	}else if (argc==7) { // ./wavy run.mac *.root random_seed NeutSttNum NeutGenNum SourceROOT
      OutName = argv[2];
      seed    = atoi(argv[3]);
      NeutSttNum= atoi(argv[4]);
      NeutGenNum= atoi(argv[5]);
      NeutName  = argv[6];
	} else {
	   //OutName = "mydata.root";
	   //G4String OutName  = "mydata_"+condition_Name+".root";
	}

	if (argv[1]!=NULL) G4cout << "input macro  name is " << argv[1] << G4endl;
	G4cout << "input rootfile name is " << OutName << G4endl;

   //std::system("ls -l "); // execute the UNIX command "ls -l >test.txt"
   std::system("rm g4_00.wrl"); // execute the UNIX command "ls -l >test.txt"
//   std::system("rm mydata.root"); // execute the UNIX command "ls -l >test.txt"
   std::system("rm mydet.gdml"); // execute the UNIX command "ls -l >test.txt"

  // ----- Choose the Random engine and set the seed
   G4cout << "\n\n\n\n ##### G4Random::setTheSeed : " << seed << " #####\n\n\n\n" << G4endl;
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(seed);

  G4RunManager *runManager = new G4RunManager;
#endif

  G4String physName = "FTFP_BERT";
  //G4String physName = "QGSP_BERT_HP";
  //G4String physName = "G4EmStandardPhysics_option2";

#ifndef WIN32
  G4int c = 0;
  while ((c=getopt(argc,argv,"p")) != -1) {
     switch (c) {
       case 'p':
         physName = optarg;
         G4cout << "Physics List used is " <<  physName << G4endl;
         break;
       case ':':       /* -p without operand */
         fprintf(stderr, "Option -%c requires an operand\n", optopt);
         break;
       case '?':
         fprintf(stderr, "Unrecognised option: -%c\n", optopt);
     }
  }
#endif

	//
  	// ----- Set mandatory initialization classes
	//
  	// Detector construction
  	WLSDetectorConstruction* detector = new WLSDetectorConstruction();
  	runManager->SetUserInitialization( detector );

  	// Physics list
  	runManager->SetUserInitialization(new WLSPhysicsList(physName));

  	// User action initialization
  	//runManager->SetUserInitialization(new WLSActionInitialization(detector,OutName));
	//WLSRunAction* runaction = new WLSRunAction(detector,OutName);
   	runManager->SetUserInitialization(new WLSActionInitialization(detector,OutName,NeutName,NeutSttNum,NeutGenNum));


#ifdef G4VIS_USE
  	// ----- Initialize visualization
  	G4VisManager* visManager = new G4VisExecutive;
  	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	G4cout << "\n\n\n\n#####";
   	G4cout << "If you get visualization, you have to turn on below";
   	G4cout << "#####\n\n\n\n" << G4endl;
  	visManager->Initialize();
#endif




	// ----- Get the pointer to the User Interface manager
  	G4UImanager * UImanager = G4UImanager::GetUIpointer();

	#ifndef WIN32
  	G4int optmax = argc;
  	if (argc > 2)  { optmax = optmax-1; }
  	if (optind < optmax) {
   	G4String command = "/control/execute ";
     	for ( ; optind < optmax; optind++) {
      	G4String macroFilename = argv[optind];
         UImanager->ApplyCommand(command+macroFilename);
     	}
  	}
	#else  // Simple UI for Windows runs, no possibility of additional arguments
  	if (argc!=1) {
   	G4String command = "/control/execute ";
     	G4String fileName = argv[1];
     	UImanager->ApplyCommand(command+fileName);
  	}
	#endif
  	else  {
   	// Define (G)UI terminal for interactive mode
		#ifdef G4UI_USE
     	G4UIExecutive * ui = new G4UIExecutive(argc,argv);
		#ifdef G4VIS_USE
   	   std::cerr << "\n\n########## G4VIS_USE: /control/execute init.in ##########\n\n" << std::endl;
     	UImanager->ApplyCommand("/control/execute init.in");
     	//UImanager->ApplyCommand("/control/execute vis.mac");
     	//UImanager->ApplyCommand("/control/execute init_vis.mac");
		#endif
     	if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
     	ui->SessionStart();
     	delete ui;
		#endif
  }

  	// ----- job termination
	#ifdef G4VIS_USE
  	delete visManager;
	#endif
  	delete runManager;

  	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
