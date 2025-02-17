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
/// \file eventgenerator/HepMC/HepMCEx01/src/WLSExN04Field.cc
/// \brief Implementation of the WLSExN04Field class
//
// $Id: WLSExN04Field.cc 77801 2013-11-28 13:33:20Z gcosmo $
//

#include "G4SystemOfUnits.hh"
#include "WLSExN04Field.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WLSExN04Field::WLSExN04Field()
{
  fBz = 2.0*tesla;
  frmax_sq = sqr(1000.*cm);
  fzmax = 1000.*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WLSExN04Field::~WLSExN04Field()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WLSExN04Field::GetFieldValue(const double Point[3],double *Bfield) const
{
#if 0
  Bfield[0] = 0.;
  Bfield[1] = 0.;
  if ( std::abs(Point[2]) < fzmax &&
      (sqr(Point[0])+sqr(Point[1])) < frmax_sq ) {
    Bfield[2] = fBz;
  } else {
    Bfield[2] = 0.;
  }
#endif
	Bfield[0]=0.2*tesla;
	Bfield[1]=0.0*tesla;
	Bfield[2]=0.0*tesla;
   //G4cerr << "Magnetic filed is +Y direction" << G4endl;
	//if( Point[2] > - 100.*cm && Point[2] < 100.*cm){
	//Bfield[0] = 1.*tesla;
	//}
}
