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
/// \file eventgenerator/HepMC/HepMCEx01/include/WLSExN04Field.hh
/// \brief Definition of the WLSExN04Field class
//
// $Id: WLSExN04Field.hh 77801 2013-11-28 13:33:20Z gcosmo $
//

#ifndef WLSExN04Field_H
#define WLSExN04Field_H 1

#include "G4MagneticField.hh"
#include "globals.hh"

class WLSExN04Field : public G4MagneticField {
	public:
  	WLSExN04Field();
  	~WLSExN04Field();

  	void GetFieldValue( const  double Point[3], double *Bfield ) const;

private:
  	G4double fBz;
  	G4double frmax_sq;
  	G4double fzmax;
};

#endif
