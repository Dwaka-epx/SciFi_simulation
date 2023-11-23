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
// $Id: WLSDetectorConstruction.cc 84718 2014-10-20 07:40:45Z gcosmo $
//
/// \file optical/wls/src/WLSDetectorConstruction.cc
/// \brief Implementation of the WLSDetectorConstruction class
//
//
#include <sstream>

#include "G4ios.hh"
#include "globals.hh"

#include "WLSExN04Field.hh"
#include "G4FieldManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4EllipticalTube.hh"
#include "G4Para.hh" // add @21/06
#include "G4Trd.hh"  // add @21/06

#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"


#include "WLSDetectorConstruction.hh"
#include "WLSDetectorMessenger.hh"
#include "WLSMaterials.hh"
#include "WLSPhotonDetSD.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//
// GDML parser include
//
#include "G4GDMLParser.hh"

using namespace std;

WLSDetectorConstruction::WLSDetectorConstruction()
  	: fMaterials(NULL), fLogiWorld(NULL), fPhysWorld(NULL)
{
  fDetectorMessenger = new WLSDetectorMessenger(this);

  fNumOfCladLayers = 1;
  fMirrorToggle = 0;//true;
  fMirrorPolish = 1.;
  fMPPCPolish   = 1.;

  fWLSfiberZ   =   1*cm; // length is 10 cm for small exp.
  fWLSfiberRY  = 0.50*mm - 0.04*mm; // phi?
  fWLSfiberOrigin = 0.0;
 
  fMPPCShape = "Square";
  fMPPCHalfL = fWLSfiberRY;
  fMPPCDist  = 0.00*mm;
  fMPPCTheta = 0.0*deg;
  fMPPCZ     = 0.05*mm;
 
  fClrfiberZ  = fMPPCZ + 10.*nm;
  fMirrorZ    = 0.1*mm;

  fBarLength        = 1.*cm;
  fBarBase          = 1.*cm;
  fCoatingThickness = 0.080*mm;
  fCoatingRadius    = 0.010*mm;

 	fHoleRadius   = 0.70*mm;
  	fHoleLength  = 1*cm;//fBarLength;

	//
	// clear tree
	//
   fReadoutX.clear();
   fReadoutYID.clear();
   fReadoutY.clear();
   fReadoutFID.clear();
   fReadoutF.clear();

	#if 1 
   //
   // Magnetic field
   //
   static G4bool fieldIsInitialized = false;
   if( !fieldIsInitialized ) {
   	WLSExN04Field* myField = new WLSExN04Field;
   	G4FieldManager* fieldMgr
      	= G4TransportationManager::GetTransportationManager() ->GetFieldManager();
    	fieldMgr->SetDetectorField(myField);
    	fieldMgr->CreateChordFinder(myField);
    	fieldIsInitialized = true;
  	}
	#endif
}


WLSDetectorConstruction::~WLSDetectorConstruction()
{
  if (fDetectorMessenger) delete fDetectorMessenger;
  if (fMaterials)         delete fMaterials;
}


G4VPhysicalVolume* WLSDetectorConstruction::Construct()
{
  if ( fPhysWorld ) {
     G4GeometryManager    ::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore ::GetInstance()->Clean();
     G4SolidStore         ::GetInstance()->Clean();
     G4LogicalSkinSurface  ::CleanSurfaceTable();
     G4LogicalBorderSurface::CleanSurfaceTable();
  }

  fMaterials = WLSMaterials::GetInstance();

  UpdateGeometryParameters();

   //return ConstructDetector();
   G4VPhysicalVolume* fPhysWorld_ret = ConstructDetector();

	//
   // GDMLparser
	//
   G4GDMLParser fParser;
   fParser.Write("mydet.gdml", fPhysWorld_ret);

   return fPhysWorld_ret;
}



G4VPhysicalVolume* WLSDetectorConstruction::ConstructDetector()
{
	//
	// world
	//
   //G4String nameWorld = "G4_AIR";
   //G4String nameWorld = "G4_Galactic";
   G4String nameWorld = "Water";

   //float fWorldSize =  1000*cm; // 10 m
   float fWorldSize =  300*cm; // 3 m
   G4cerr << ">>>>> World: fWorldSize=" << fWorldSize << G4endl;

  	G4VSolid* solidWorld(0);
	solidWorld = new G4Box("World", fWorldSize, fWorldSize, fWorldSize);
  	fLogiWorld = new G4LogicalVolume(solidWorld, FindMaterial(nameWorld), "World");
  	fPhysWorld = new G4PVPlacement(0,G4ThreeVector(),fLogiWorld,"World",0,false,0);

   // scCoat
   double thickness = GetBarBase()/2 + fCoatingThickness;// + GetCoatingRadius();
   //G4cerr << "GetBarBase()=" << GetBarBase() << " GetBarLength()=" << GetBarLength() << G4endl;
   //G4cerr << ">>> ScCoat: thickness = " << thickness - GetCoatingThickness() << G4endl;
   //G4cerr << ">>> ScCube: thickness = " << GetBarBase()/2 * 2 << G4endl;

   G4VSolid* solidScCoat(0);
   solidScCoat = new G4Box("ScCoat",thickness,thickness,thickness);

   G4VSolid* solidScCube(0);
   solidScCube = new G4Box("ScCube", GetBarBase()/2,GetBarBase()/2,GetBarBase()/2);
   fLogiScCoat = new G4LogicalVolume(solidScCoat, FindMaterial("Coating"),"ScCoat");
   flogiScCube = new G4LogicalVolume(solidScCube, FindMaterial("Polystyrene"), "ScCube");

	//
	// define surface and table of surface properties table
	//
	/*
 	G4OpBoundaryProcess クラスを用いるときには model を設定してやる必要がある。
	この model には UNIFIED model と GLISUR model の 2 種類があり、
	GLISUR model は、Gean- t3.21(Geant4 はこの改訂版) の model であるが、
	UNIFIED model はこの GLISUR model より もリアルなシミュレーションを得ようとしており、
	表面の磨きや反射剤の全ての様子を取り扱 う model である。
	https://oxon.hatenablog.com/entry/20100121/1264056313 Nagoya, Okumura

	Lambertian reflectance ( based on Kikawa-san's optical simulation )

	groundfrontpainted 反射の仕方が常にランバート反射 
	   : 反射の仕方が常にランバート反射 
       polishとsigma_alphaの両パラメータは無視されます。

	groundbackpainted     
     polishedbackpaintedと計算過程はほぼ同じですが、反射材での反射は常にランバート反射です。
     その一方で、境界面での反射や屈折にはpolishもしくはsigma_alphaが使われます。
     従って、この境界面での反射の分布は、反射材でのランバート反射と境界面での表面粗さの組み合わさったものになります。
	*/
  	G4OpticalSurface* TiO2Surface(0);
	TiO2Surface = new G4OpticalSurface("TiO2Surface",
													// this largely changes transmission sistance 
                                      //glisur, ground, dielectric_metal, 1); // SetModel SetFinish SetType SetPolish
												  unified, polished, dielectric_metal, 0.98); // SetModel SetFinish SetType SetPolish
												  //unified, groundbackpainted, dielectric_metal, 0.9); // SetModel SetFinish SetType SetPolish
												  //unified, groundbackpainted, dielectric_metal, 0.98); // SetModel SetFinish SetType SetPolish

  	G4MaterialPropertiesTable* TiO2SurfaceProperty = new G4MaterialPropertiesTable();

  	G4double p_TiO2[] = {2.00*eV, 3.47*eV};
  	const int nbins = sizeof(p_TiO2)/sizeof(G4double);
  	G4double refl_TiO2[] = {0.97, 0.97};
  	assert(sizeof(refl_TiO2) == sizeof(p_TiO2));
  	G4double effi_TiO2[] = {0, 0};
  	assert(sizeof(effi_TiO2) == sizeof(p_TiO2));
   G4double rind_TiO2[] = {2.2, 2.2};
   assert(sizeof(effi_TiO2) == sizeof(p_TiO2));
	// Boundary Process
	// p210 ~ http://ftp.jaist.ac.jp/pub/Linux/Gentoo/distfiles/BookForAppliDev-4.10.1.pdf
 	// reflectivity : the ratio of reflection
  	// efficiency : the ratio of absorption & energy deposit on material
  	// At dielectric-metal border, first determine reflect or not with reflectivity,
  	// then determine deposit or not energy on material with efficiency.
  	// reflectivity= 0 means no reflection
  	// efficiency  = 0 means no energy deposite on material
  	TiO2SurfaceProperty -> AddProperty("REFLECTIVITY", p_TiO2,  refl_TiO2  ,nbins);
  	TiO2SurfaceProperty -> AddProperty("EFFICIENCY",   p_TiO2,  effi_TiO2,  nbins);
  	TiO2SurfaceProperty -> AddProperty("RINDEX",       p_TiO2,  rind_TiO2,  nbins);
	// SetSigmaAlpha : 
	// http://wiki.opengatecollaboration.org/index.php/Users_Guide:Generating_and_tracking_optical_photons
	// This parameter defines the standard deviation of the Gaussian distribution 
	//	of micro-facets around the average surface normal
  	TiO2Surface -> SetSigmaAlpha(0.0);
  	TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);

	//
	// ScCoat surface
	//
  	new G4LogicalSkinSurface("TiO2Surface", fLogiScCoat, TiO2Surface);

	//
	// Boundary Surface Properties between ScCube and 
	//
	G4OpticalSurface *opSurface = new G4OpticalSurface("RoughSurface",
                                      glisur, ground, dielectric_dielectric, 0.99); // SetModel SetFinish SetType SetPolish 
												  //unified, polished, dielectric_dielectric, 0.99); //  SetModel SetFinish SetType SetPolish
   //new G4LogicalBorderSurface("surfaceHoleXOt", fPhysWorld, fPhysScCoat, opSurface); 
   //new G4LogicalBorderSurface("surfaceHoleXIn", fPhysScCoat, fPhysWorld, opSurface); 



	//
	// from here, detector description @ 21/11/15 
	//
	const float fiber_length = 200 * mm; // mm : standard scinti layer
	//const float fiber_length = 500 * mm; // mm : standard scinti layer
	//const float fiber_length = 1000 * mm; // mm : standard scinti layer
	const float fiber_thickness = 1 * mm; // mm : standard scinti layer
	//const float fiber_thickness = 2 * mm; // mm : standard scinti layer
   const float coat_thickness = 0.05 * mm; // mm
   //const float gap_thickness = 100 * mm; // mm
   const float gap_thickness = 0.01 * mm; // mm

   const int ite_inside_unit= 1;
   
   const int kinds_of_layers = 4;
   const int &unit = kinds_of_layers;
   const int nLayersZ = 50;//12*unit; //number of xyuv layers arranged along with Z direction.

   G4RotationMatrix* rotMY= new G4RotationMatrix;
   G4RotationMatrix* rotMU= new G4RotationMatrix;
   G4RotationMatrix* rotMV= new G4RotationMatrix;
   
   rotMU-> rotateZ(-45.*deg);
   rotMY-> rotateZ(-90.*deg);
   rotMV-> rotateZ(-135.*deg);

	/*
   G4VSolid        *solFiberDetMother = new G4Box("FiberDetMother", 100*cm, 100*cm, fiber_thickness * nLayersZ /2); 
   G4LogicalVolume* logFiberDetMother = new G4LogicalVolume(solFiberDetMother, FindMaterial(nameWorld), "FiberDetMother");
   new G4PVPlacement(0,     G4ThreeVector(  0,  0, 0), logFiberDetMother, "FiberDetMother", fLogiWorld, false, 0);
	*/

   // set layer mother volume in turns
	// dx : tentatively large number
	// dy : tentatively large number
	// dz : 2mm (fiber_thickness) / 2 
   float xy = fiber_length * 1.2; // just set it bigger
   G4VSolid        *solFiberLayerMother       = new G4Box("FiberLayerMother",       xy/2, xy/2, fiber_thickness / 2.); 
   G4LogicalVolume* logFiberLayerMother       = new G4LogicalVolume(solFiberLayerMother, FindMaterial(nameWorld), "FiberLayerMother");

   G4VSolid        *solFiberLayerSparceMother = new G4Box("FiberLayerSparceMother", xy/2, xy/2, fiber_thickness / 2.); 
   G4LogicalVolume* logFiberLayerSparceMother = new G4LogicalVolume(solFiberLayerSparceMother, FindMaterial(nameWorld), "FiberLayerSparceMother");

   float distance_between_layers = 10;
   float preZ = 0; // offset ; position where xuyv layers are set.

   for (int i = 0; i < nLayersZ; ++i) { // + 2 is x/y in the end //fiber unit loop for xuyv
      //G4PVPlacement   *phyFiberLayerMother;
      //new G4PVPlacement(0, G4ThreeVector(0,0,200*mm), logFiberLayerMother, "FiberLayerMother-1", fLogiWorld, false, 0);
      //new G4PVPlacement(0, G4ThreeVector(0,0,600*mm), logFiberLayerMother, "FiberLayerMother-2", fLogiWorld, false, 0);
      stringstream name1;
      name1 << "FiberLayerMother." << i << ends;
      stringstream name2;
      name2 << "FiberLayerSparceMother." << i << ends;

		//
		// must remember  +90 degrees  around +Z for X-layre
                //

		#if 1 // if nLayersZ%4=0,then X layer. if nLayersZ%4=1,then U layer. if nLayersZ%4=2,then Y layer. if nLayersZ%4=3,then V layer.
      if ( i!=0 ) preZ = preZ + distance_between_layers;
      if      ( i%kinds_of_layers==0 ) new G4PVPlacement(0,     G4ThreeVector(     0,    0, preZ), logFiberLayerMother,       name1.str().data(), fLogiWorld, false, 0); // x layer
      else if ( i%kinds_of_layers==1 ) new G4PVPlacement(rotMU, G4ThreeVector(     0,    0, preZ), logFiberLayerMother,       name1.str().data(), fLogiWorld, false, 0);
      else if ( i%kinds_of_layers==2 ) new G4PVPlacement(rotMY, G4ThreeVector(     0,    0, preZ), logFiberLayerMother,       name1.str().data(), fLogiWorld, false, 0);
      else if ( i%kinds_of_layers==3 ) new G4PVPlacement(rotMV, G4ThreeVector(     0,    0, preZ), logFiberLayerMother,       name1.str().data(), fLogiWorld, false, 0);
      //      else                             new G4PVPlacement(rotMZ, G4ThreeVector(     0,    0, preZ), logFiberLayerSparceMother, name2.str().data(), fLogiWorld, false, 0);
		#endif
   }


	#if  1 
   //
	// straight fiber layers x/y
	//
	DefineLayerStraightFiber(logFiberLayerMother, 
                            fiber_length, fiber_thickness, coat_thickness);
   //
   // straight fiber layers x/y + straight fiber layers along x aligned to z (ND280.211124.pdf by Ichikawa) 
   //
	// DefineLayerSparceFiber(logFiberLayerSparceMother, fiber_length, fiber_thickness, coat_thickness);
   #endif
	#if 0
	//
	// waving fiber layers
	// this is not good now @ 21/11/25	
	DefineLayerWavingFiber(logFiberLayerMother, fiber_length, fiber_thickness, coat_thickness, gap_thickness, ite_inside_unit);
	#endif

	//
	// define  visualization
	//
   fLogiWorld         ->SetVisAttributes(G4VisAttributes::Invisible);
   logFiberLayerMother->SetVisAttributes(G4VisAttributes::Invisible);
   G4VisAttributes* va1 = new G4VisAttributes(G4Colour(0.9,0.0,0.0)); 
   G4VisAttributes* va2 = new G4VisAttributes(G4Colour(0.0,0.9,0.0));
   logFiberLayerMother->SetVisAttributes(va1);
   logFiberLayerSparceMother->SetVisAttributes(va2);
 	return fPhysWorld;
}



void WLSDetectorConstruction::DefineLayerStraightFiber(G4LogicalVolume* logFiberLayerMother, 
	double fiber_length, double fiber_thickness, double coat_thickness)
{
   //
   // straight fiber layers
   //
   G4VSolid *SoliFbrCoat = new G4Box("FbrCoat",  fiber_thickness/2.,                 fiber_length/2.,  fiber_thickness/2.);
   G4VSolid *SoliFbrScin = new G4Box("FbrScin", (fiber_thickness-coat_thickness)/2., fiber_length/2., (fiber_thickness-coat_thickness)/2.);

   // subtraction: 
   // this must be done before constructing logical volume
   //SoliFbrCoat = new G4SubtractionSolid("subtractSoliFbrCoat", SoliFbrCoat, SoliFbrScin, 0, G4ThreeVector(0,0,0));

   G4LogicalVolume   *LogiFbrCoat = new G4LogicalVolume(SoliFbrCoat, FindMaterial("Coating"),     "FbrCoat");
   G4LogicalVolume   *LogiFbrScin = new G4LogicalVolume(SoliFbrScin, FindMaterial("Polystyrene"), "FbrScin");
   //new G4LogicalSkinSurface("TiO2Surface_LogiFbrCoat", LogiFbrCoat, TiO2Surface);

   int nfibers = fiber_length / fiber_thickness; 

   for (int x=0; x<nfibers; x++) {  // x direction
      const int id = x;
      G4cerr << "DefineLayerStraightFiber fiber id = " << id << G4endl;
      // e.g.  area is 300 mm, / thickness (2mm) = 150 fibers
      //      
      float pos_x = - fiber_length / 2 + fiber_thickness * x;
      new G4PVPlacement(0, G4ThreeVector(pos_x, 0, 0), LogiFbrCoat, "FbrCoat", logFiberLayerMother, false, id);

      // to study/take dtata
      //new G4PVPlacement(0, G4ThreeVector(), LogiFbrScin, "FbrScin", LogiFbrCoat, false, id);
      // above one does not give proper id... why...
      // add: possible to access via placement tree
      //new G4PVPlacement(0, vec1, LogiFbrScin, "FbrScin", fLogiWorld, false, id);

      fReadoutYID.push_back(id);
      fReadoutY.push_back( TVector3(pos_x, 0, 0) );
   }

   // to visualize geometry in the display, "just draw or /vis/open VRML2FILE" turn on below G4PVPlacement and off above G4PVPlacement
   // to study/take dtata,                                                     turn on above G4PVPlacement and off below G4PVPlacement 
   new G4PVPlacement(0, G4ThreeVector(), LogiFbrScin, "FbrScin", LogiFbrCoat, false, 0);

   G4VisAttributes* va1 = new G4VisAttributes(G4Colour(0.0,0.9,0.9)); // cyan 
   G4VisAttributes* va2 = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
   LogiFbrCoat->SetVisAttributes(va1);
   LogiFbrScin->SetVisAttributes(va2);
}

/*
void WLSDetectorConstruction::DefineLayerSparceFiber(G4LogicalVolume* logFiberLayerSparceMother, 
	double fiber_length, double fiber_thickness, double coat_thickness)
{
   //
   // straight fiber layers
   //
   G4VSolid *SoliFbrCoat = new G4Box("FbrCoat",  fiber_thickness/2.,                 fiber_length/2.,  fiber_thickness/2.);
   G4VSolid *SoliFbrScin = new G4Box("FbrScin", (fiber_thickness-coat_thickness)/2., fiber_length/2., (fiber_thickness-coat_thickness)/2.);

   // subtraction: 
   // this must be done before constructing logical volume
   //SoliFbrCoat = new G4SubtractionSolid("subtractSoliFbrCoat", SoliFbrCoat, SoliFbrScin, 0, G4ThreeVector(0,0,0));

   G4LogicalVolume   *LogiFbrCoat = new G4LogicalVolume(SoliFbrCoat, FindMaterial("Coating"),     "FbrCoat");
   G4LogicalVolume   *LogiFbrScin = new G4LogicalVolume(SoliFbrScin, FindMaterial("Polystyrene"), "FbrScin");
   //new G4LogicalSkinSurface("TiO2Surface_LogiFbrCoat", LogiFbrCoat, TiO2Surface);

   int nfibers = fiber_length / 30 * mm; 

   for (int y=0; y< nfibers + 1; y++) {  // y direction
      const int id = y;
		G4cerr << "DefineLayerSparceFiber fiber id = " << id << G4endl;
      float pos_y = - fiber_length / 2  + 30 * mm * y;
      new G4PVPlacement(0, G4ThreeVector(pos_y, 0, 0), LogiFbrCoat, "FbrCoat", logFiberLayerSparceMother, false, id);
      fReadoutFID.push_back(id);
      fReadoutF.push_back( TVector3(pos_y, 0, 0) );
	}

   new G4PVPlacement(0, G4ThreeVector(), LogiFbrScin, "FbrScin", LogiFbrCoat, false, 0);
   G4VisAttributes* va1 = new G4VisAttributes(G4Colour(0.9,0.9,0.0)); // cyan 
   G4VisAttributes* va2 = new G4VisAttributes(G4Colour(0.8,0.0,0.8));
   LogiFbrCoat->SetVisAttributes(va1);
   LogiFbrScin->SetVisAttributes(va2);
}


void WLSDetectorConstruction::DefineLayerWavingFiber(G4LogicalVolume* logFiberLayerMother,
   double fiber_length, double fiber_thickness, double coat_thickness, double gap_thickness,
   int ite_inside_unit)
{
   //const float waveAmpli  = 10*mm; // mm
   const float waveAmpli  = 24*mm; // mm
   const float waveLengh  = 50*mm; // mm

   const float jointYHight= 1*mm; // mm : wavy scinti layer
   const float jointXWidth= 2*mm; // mm : wavy scinti layer
   const float jointZWidth= 2*mm; // mm : wavy scinti layer
   const float ObliBoxYHight = waveLengh/2 - jointYHight;

   // #of unit in x direction
   //const int layersX = 2; // ampli 50 mm, y 10 cm 
   //const int layersX = 6; // for ampli 50 mm, y 30 cm 
   //const int layersX = 9; // ampli 50 mm, y 50 cm 
   const int layersX = 18; // ampli 50 mm, y 100 cm

   // #of unit in y direction
   //const int layersY = 1; // 10 cm
   //const int layersY = 3; // for ampli 50 mm, y 30 cm
   //const int layersY = 5; // 50 cm
   const int layersY = 10; // ampli 50 mm, y 100 cm

   G4RotationMatrix* rotMZ= new G4RotationMatrix;
   G4RotationMatrix* rotMY= new G4RotationMatrix;
   G4RotationMatrix* rotMX= new G4RotationMatrix;
   rotMX-> rotateX(180.*deg);
   rotMY-> rotateX(180.*deg);
   rotMZ-> rotateZ( 90.*deg);

   G4RotationMatrix *rotMXY = new G4RotationMatrix; // redefine new one
   rotMXY-> rotateX(180.*deg);
   rotMXY-> rotateX(180.*deg);

	//
   // angle for tapezoidal box joint 
	//
   float angle     = atan( ObliBoxYHight/waveAmpli );
   float angleDegs = atan( ObliBoxYHight/waveAmpli ) / M_PI * 180;
   G4cerr << ">>>>> angle " << angleDegs << G4endl;
   G4cerr << ">>>>> angle " << (90-angleDegs) << G4endl;

   // box joint (including coating)
   G4VSolid *SoliJoint1 = new G4Box("Joint1", jointXWidth/2., jointYHight/2., jointZWidth/2.);

   // tapezoidal box joint (including coating)
   G4VSolid *SoliJoint2 = new G4Para("Joint2", jointXWidth/2., ObliBoxYHight/2., jointZWidth/2., (90-angleDegs)*deg, 0, 0);

   //G4VSolid *SoliSct1 = new G4Box("Sct1", (jointXWidth-coat_thickness)/2., (jointYHight-coat_thickness)/2., (jointZWidth-coat_thickness)/2.);
   G4VSolid *SoliSct1 = new G4Box("Sct1", (jointXWidth-coat_thickness)/2., jointYHight/2., (jointZWidth-coat_thickness)/2.);

   //G4VSolid *SoliSct2 = new G4Para("Sct2", (jointXWidth-coat_thickness)/2., (ObliBoxYHight-coat_thickness)/2., (jointZWidth-coat_thickness)/2., (90-angleDegs)*deg, 0, 0);
   G4VSolid *SoliSct2 = new G4Para("Sct2", (jointXWidth-coat_thickness)/2., ObliBoxYHight/2., (jointZWidth-coat_thickness)/2., (90-angleDegs)*deg, 0, 0);

   // subtraction: 
   // this must be done before constructing logical volume
   //SoliJoint1 = new G4SubtractionSolid("subtractSoliJoint1", SoliJoint1, SoliSct1, 0, G4ThreeVector(0,0,0));
   //SoliJoint2 = new G4SubtractionSolid("subtractSoliJoint2", SoliJoint2, SoliSct2, 0, G4ThreeVector(0,0,0));

   // then, put logical volume
   G4LogicalVolume   *LogiJoint1 = new G4LogicalVolume(SoliJoint1, FindMaterial("Coating"), "Joint1");
   G4LogicalVolume   *LogiJoint2 = new G4LogicalVolume(SoliJoint2, FindMaterial("Coating"), "Joint2");

   //new G4LogicalSkinSurface("TiO2Surface_LogiJoint1", LogiJoint1, TiO2Surface);
   //new G4LogicalSkinSurface("TiO2Surface_LogiJoint2", LogiJoint2, TiO2Surface);

   G4LogicalVolume   *LogiSct1 = new G4LogicalVolume(SoliSct1, FindMaterial("Polystyrene"), "Sct1");
   G4LogicalVolume   *LogiSct2 = new G4LogicalVolume(SoliSct2, FindMaterial("Polystyrene"), "Sct2");

   //for (int z=0; z < nLayersZ * ite_inside_unit; z++) {  // z direction
   for (int z=0; z < ite_inside_unit; z++) {  // z direction
      for (int x=0; x < layersX; x++) {  // x direction
         for (int y=0; y < layersY; y++) {  // y direction

         // y-direction is same id
         int id = 2*(layersX+1) * z + 2*x;
         //                        ( 50mm * 6 + 2mm * 6 ) / 2                   
			float global_Xshift  = - ( waveLengh*layersX + jointXWidth*layersX )/2;
         float ite_for_Xshift = jointXWidth + waveAmpli + waveAmpli + jointXWidth; 
         float ite_for_Yshift = jointYHight + waveLengh + jointYHight; 
         float ite_for_Zshift = fiber_thickness + gap_thickness + 
                                fiber_thickness + gap_thickness + 
                               (jointZWidth * ite_inside_unit)+ gap_thickness;

         // start (mother volume)  - ite_inside_unit * 2 / 2 *mm = -24mm
         //float positionZ = 2 * ( fiber_thickness + gap_thickness ) + 
         //                  jointZWidth * (z%ite_inside_unit) + jointZWidth/2 + 
         //                  ite_for_Zshift * (z/ite_inside_unit);
         //float positionZ = -24 *mm 
         float positionZ = -1 *mm  + 2*mm/2  + 2*mm * z;
      
         G4ThreeVector vec1(            jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                        jointYHight/2 + ite_for_Yshift * y, 
                                        positionZ );
 
         G4ThreeVector vec2(waveAmpli + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                        jointYHight + ObliBoxYHight + jointYHight/2 + ite_for_Yshift * y, 
                                        positionZ );

         G4ThreeVector vec3(waveAmpli + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                        jointYHight + ObliBoxYHight + jointYHight + jointYHight/2 + ite_for_Yshift * y,        
                                        positionZ );

         G4ThreeVector vec4(            jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                        jointYHight + ObliBoxYHight + jointYHight + jointYHight + ObliBoxYHight + jointYHight/2 + ite_for_Yshift * y, 
                                        positionZ );
 
         G4ThreeVector vec5((ObliBoxYHight/2)/tan(angle) + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift,                  
                                        jointYHight                                             + ObliBoxYHight/2 + ite_for_Yshift * y, 
                                        positionZ );

         G4ThreeVector vec6((ObliBoxYHight/2)/tan(angle) + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                        jointYHight + ObliBoxYHight + jointYHight + jointYHight + ObliBoxYHight/2 + ite_for_Yshift * y, 
                                        positionZ );

         G4cerr << "UPPER-NORMAL: nth-z (in each mother) is " << z  << " its object id is " << id << "   CONFIRM from ID:"
                << " nth-z = "          << id / ( (layersX+1)*2 ) 
                << " nth-x in nth-z = " << id % ( (layersX+1)*2 )
                << " pos: " << vec1
                << G4endl;
         // set the fiber components in mother volume
         new G4PVPlacement(0,     vec1, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec2, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec3, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec4, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec5, LogiJoint2, "Joint2", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMX, vec6, LogiJoint2, "Joint2", logFiberLayerMother, false, id);

         #if 1
         // 100mm: 100 fibers 
         // 1000mm: 1000 fibers 
         //for (int j = 0; j < nfibers_in_y_layer; ++j) {
         //   float sttx = ( nfibers_in_y_layer / 2 ) * -1 *mm;
         //   float x = float( sttx + (j * 1*mm) ); // start-x + no. x width 
         if (y==0) {
            fReadoutYID.push_back(id);
            fReadoutY.push_back( TVector3(jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift,
                                        jointYHight/2 + ite_for_Yshift * y,
                                        positionZ) );
         }
         #endif

         vec1.setY( vec1.getY()*-1 );
         vec2.setY( vec2.getY()*-1 );
         vec3.setY( vec3.getY()*-1 );
         vec4.setY( vec4.getY()*-1 );
         vec5.setY( vec5.getY()*-1 );
         vec6.setY( vec6.getY()*-1 );

         new G4PVPlacement(0,     vec1, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec2, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec3, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec4, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMX, vec5, LogiJoint2, "Joint2", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec6, LogiJoint2, "Joint2", logFiberLayerMother, false, id);
         }
      }
   }
	
   // turn around 
   //for (int z=0; z<nLayersZ*ite_inside_unit; z++) {  // z direction
   for (int z=0; z<ite_inside_unit; z++) {  // z direction
   //for (int z=0; z<5; z++) {  // z direction
      for (int x=0; x<layersX; x++) {  // x direction
         for (int y=0; y<layersY; y++) {  // y direction

         // y-direction is same id
         int id = 2*layersX * z + 2*x + 1;
         //                        ( 50mm * 6 + 2mm * 6 ) / 2                   
			float global_Xshift  = - ( waveLengh*layersX + jointXWidth*layersX )/2;
         float ite_for_Xshift = jointXWidth + waveAmpli + waveAmpli + jointXWidth;
         float ite_for_Yshift = jointYHight + waveLengh + jointYHight; 
         float ite_for_Zshift = fiber_thickness + gap_thickness +
                                fiber_thickness + gap_thickness +
                               (jointZWidth * ite_inside_unit)+ gap_thickness;

         // start (mother volume)  - ite_inside_unit * 2 / 2 *mm = -24mm
         //float positionZ = 2 * ( fiber_thickness + gap_thickness ) + 
         //                  jointZWidth * (z%ite_inside_unit) + jointZWidth/2 + 
         //                  ite_for_Zshift * (z/ite_inside_unit);
         //float positionZ = -24 *mm 
         float positionZ = -1 *mm  + 2*mm/2  + 2*mm * z;

         G4ThreeVector vec1(           jointXWidth + waveAmpli + waveAmpli + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight/2 + ite_for_Yshift * y, 
                                       positionZ );

         G4ThreeVector vec2(waveAmpli +jointXWidth                         + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight + ObliBoxYHight + jointYHight/2 + ite_for_Yshift * y, 
                                       positionZ );

         G4ThreeVector vec3(waveAmpli +jointXWidth                         + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight + ObliBoxYHight + jointYHight + jointYHight/2 + ite_for_Yshift * y, 
                                       positionZ );

         G4ThreeVector vec4(           jointXWidth + waveAmpli + waveAmpli + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight + ObliBoxYHight + jointYHight + jointYHight + ObliBoxYHight + jointYHight/2 + ite_for_Yshift * y, 
                                       positionZ );
   
         G4ThreeVector vec5(jointXWidth + waveAmpli + (ObliBoxYHight/2)/tan(angle) + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight                                             + ObliBoxYHight/2 + ite_for_Yshift * y,
                                       positionZ );

         G4ThreeVector vec6(jointXWidth + waveAmpli + (ObliBoxYHight/2)/tan(angle) + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift, 
                                       jointYHight + ObliBoxYHight + jointYHight + jointYHight + ObliBoxYHight/2 + ite_for_Yshift * y, 
                                       positionZ );

         G4cerr << "UPPER-REVERSE: nth-z (in each mother) is " << z  << " its object id is " << id << "   CONFIRM from ID:"
                << " nth-z = "          << id / ( (layersX)*2 )
                << " nth-x in nth-z = " << id % ( (layersX)*2 )
                << " pos: " << vec1
                << G4endl;

         new G4PVPlacement(0,     vec1, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec2, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec3, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec4, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMY, vec5, LogiJoint2, "Joint2", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMXY,vec6, LogiJoint2, "Joint2", logFiberLayerMother, false, id);

         #if 1
         // 100mm: 100 fibers 
         // 1000mm: 1000 fibers 
         //for (int j = 0; j < nfibers_in_y_layer; ++j) {
         //   float sttx = ( nfibers_in_y_layer / 2 ) * -1 *mm;
         //   float x = float( sttx + (j * 1*mm) ); // start-x + no. x width 
         if (y==0) {
            fReadoutYID.push_back(id);
            fReadoutY.push_back( TVector3(jointXWidth + waveAmpli + waveAmpli + jointXWidth/2 + ite_for_Xshift * x  +  global_Xshift,
                                       jointYHight/2 + ite_for_Yshift * y,
                                       positionZ) );
         }
         #endif

         vec1.setY( vec1.getY()*-1 );
         vec2.setY( vec2.getY()*-1 );
         vec3.setY( vec3.getY()*-1 );
         vec4.setY( vec4.getY()*-1 );
         vec5.setY( vec5.getY()*-1 );
         vec6.setY( vec6.getY()*-1 );

         new G4PVPlacement(0,     vec1, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec2, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec3, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(0,     vec4, LogiJoint1, "Joint1", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMXY,vec5, LogiJoint2, "Joint2", logFiberLayerMother, false, id);
         new G4PVPlacement(rotMX, vec6, LogiJoint2, "Joint2", logFiberLayerMother, false, id);

         }
      }
   }

   // to visualize geometry in the display, "just draw or /vis/open VRML2FILE" turn on below G4PVPlacement and off above G4PVPlacement
   // to study/take dtata,                                                     turn on above G4PVPlacement and off below G4PVPlacement 
   new G4PVPlacement(0, G4ThreeVector(), LogiSct1, "Sct1", LogiJoint1, false, 0);
   new G4PVPlacement(0, G4ThreeVector(), LogiSct2, "Sct2", LogiJoint2, false, 0);

   G4VisAttributes* va1 = new G4VisAttributes(G4Colour(0.8,0.0,0.8)); // Magenta 
   LogiJoint1->SetVisAttributes(va1);
   LogiJoint2->SetVisAttributes(va1);
   G4VisAttributes* va2 = new G4VisAttributes(G4Colour(0.8,0.8,0.8)); 
   LogiSct1->SetVisAttributes(va2);
   LogiSct2->SetVisAttributes(va2);
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::UpdateGeometryParameters()
{
  //fWLSfiberRX  = fXYRatio * fWLSfiberRY;
  fWLSfiberRX  = fWLSfiberRY;

  fClad1RX = fWLSfiberRX + 0.03*fWLSfiberRX;
  fClad1RY = fWLSfiberRY + 0.03*fWLSfiberRY;
  fClad1Z  = fWLSfiberZ;

  fClad2RX = fClad1RX + 0.03*fWLSfiberRX;
  fClad2RY = fClad1RY + 0.03*fWLSfiberRY;
  fClad2Z  = fWLSfiberZ;

  fWorldSizeX = fClad2RX   + fMPPCDist + fMPPCHalfL + 100.*cm;
  fWorldSizeY = fClad2RY   + fMPPCDist + fMPPCHalfL + 100.*cm;
  fWorldSizeZ = fWLSfiberZ + fMPPCDist + fMPPCHalfL + 100.*cm;
 
  fCoupleRX = fWorldSizeX;
  fCoupleRY = fWorldSizeY;
  fCoupleZ  = (fWorldSizeZ - fWLSfiberZ) / 2;
 
  fClrfiberHalfL = fMPPCHalfL;
 
  fMirrorRmax = fClad2RY;
 
  fCoupleOrigin = fWLSfiberOrigin + fWLSfiberZ + fCoupleZ;
  fMirrorOrigin = fWLSfiberOrigin - fWLSfiberZ - fMirrorZ;
  fMPPCOriginX  = std::sin(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
  fMPPCOriginZ  = -fCoupleZ + std::cos(fMPPCTheta) * (fMPPCDist + fClrfiberZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix
            WLSDetectorConstruction::StringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).

  G4RotationMatrix rot;

  unsigned int place = 0;

  while (place < rotation.size()) {

        G4double angle;
        char* p;

        const G4String tmpstring=rotation.substr(place+1);

        angle = strtod(tmpstring.c_str(),&p) * deg;
 
        if (!p || (*p != (char)',' && *p != (char)'\0')) {
           G4cerr << "Invalid rotation specification: " <<
                                                  rotation.c_str() << G4endl;
           return rot;
        }

        G4RotationMatrix thisRotation;

        switch(rotation.substr(place,1).c_str()[0]) {
              case 'X': case 'x':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
                break;
              case 'Y': case 'y':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
                break;
              case 'Z': case 'z':
                thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
                break;
              default:
                G4cerr << " Invalid rotation specification: "
                       << rotation << G4endl;
                return rot;
        }

       rot = thisRotation * rot;
       place = rotation.find(',',place);
       if (place > rotation.size()) break;
       ++place;
  }

  return rot;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetGeometry (G4String shape)
// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
{
  if (shape == "Circle" || shape == "Square" ) fMPPCShape = shape;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetNumberOfCladding(G4int num)
// Set the number of claddings
// Pre: 0 <= num <= 2
{
  fNumOfCladLayers = num;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetWLSLength (G4double length)
// Set the TOTAL length of the WLS fiber
{
  fWLSfiberZ = length;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetWLSRadius (G4double radius)
// Set the Y radius of WLS fiber
{
  fWLSfiberRY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetClad1Radius (G4double radius)
// Set the Y radius of Cladding 1
{
  fClad1RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetClad2Radius (G4double radius)
// Set the Y radius of Cladding 2
{
  fClad2RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
{
  fMPPCHalfL = halfL;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetGap (G4double gap)
// Set the distance between fiber end and PhotonDet
{ 
  fMPPCDist = gap;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetAlignment(G4double theta)
// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
{
  fMPPCTheta = theta;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetSurfaceRoughness(G4double roughness)
// Set the Surface Roughness between Cladding 1 and WLS fiber
// Pre: 0 < roughness <= 1
{
  fSurfaceRoughness = roughness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirrorPolish(G4double polish)
// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMirrorPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMirrorReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetPolish(G4double polish)
// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMPPCPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMPPCReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetMirror(G4bool flag)
// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
{
  fMirrorToggle = flag;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetXYRatio(G4double r)
// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
{
  fXYRatio = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetBarLength (G4double length)
// Set the length of the scintillator bar
{
  fBarLength = length;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetBarBase (G4double side)
// Set the side of the scintillator bar
{
  fBarBase = side;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetHoleRadius (G4double radius)
// Set the radius of the fiber hole
{
  fHoleRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetCoatingThickness (G4double thick)
// Set thickness of the coating on the bars
{
  fCoatingThickness = thick;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSDetectorConstruction::SetCoatingRadius (G4double radius)
// Set inner radius of the corner bar coating
{
  fCoatingRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

G4double WLSDetectorConstruction::GetWLSFiberLength() { return fWLSfiberZ;      }
G4double WLSDetectorConstruction::GetBarLength()      { return fBarLength;      }
G4double WLSDetectorConstruction::GetBarBase()        { return fBarBase;        }
G4double WLSDetectorConstruction::GetHoleRadius()     { return fHoleRadius;     }
G4double WLSDetectorConstruction::GetHoleLength()     { return fHoleLength;     }
G4double WLSDetectorConstruction::GetFiberRadius()    { return GetWLSFiberRMax(); }
G4double WLSDetectorConstruction::GetWLSFiberEnd()    { return fWLSfiberOrigin + fWLSfiberZ; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double WLSDetectorConstruction::GetWLSFiberRMax()
{
  if (fNumOfCladLayers == 2) return fClad2RY;
  if (fNumOfCladLayers == 1) return fClad1RY;
  return fWLSfiberRY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Return True if the fiber construction is ideal
G4bool WLSDetectorConstruction::IsPerfectFiber()
{
  return     fSurfaceRoughness == 1. && fXYRatio == 1.
             && (!fMirrorToggle    ||
             (fMirrorPolish    == 1. && fMirrorReflectivity == 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* WLSDetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
