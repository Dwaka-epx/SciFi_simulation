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
// $Id: WLSActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file optical/wls/include/WLSActionInitialization.hh
/// \brief Definition of the WLSActionInitialization class

#ifndef WLSActionInitialization_h
#define WLSActionInitialization_h 1

#include "G4String.hh"
#include "G4VUserActionInitialization.hh"
#include "WLSRunAction.hh"


class WLSDetectorConstruction;

class B4DetectorConstruction;

class G4GeneralParticleSource;

/// Action initialization class.
///

class WLSActionInitialization : public G4VUserActionInitialization
{
    public:
    //WLSActionInitialization(WLSDetectorConstruction*);
    //WLSActionInitialization(WLSDetectorConstruction*,G4String);
    //WLSActionInitialization(WLSDetectorConstruction*, G4String, G4String, int, int);
    WLSActionInitialization(WLSDetectorConstruction*, G4String, G4String, int, int);
    virtual ~WLSActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
 
    virtual G4VSteppingVerbose* InitializeSteppingVerbose() const;

    private:
    WLSDetectorConstruction* fDetector;
    G4String fName;
    G4String _NeutName;
    int _NeutSttNum;
    int _NeutGenNum;
};

#endif
