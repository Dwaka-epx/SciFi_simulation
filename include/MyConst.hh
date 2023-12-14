
#ifndef MyConst_h

#define MyConst_h 1
	//****************** General ********************
  const int ndim =3;
	//****************** NEUT ********************
	#define USE_NEUT 0 // flag for NEUT STUDY or PROTON STUDY
	const	int nmaxParticles = 20;
  const bool flagProtonRejection = false; // flag for ML; erase proton.
  const bool flagOnlyProton = false; // flag for ML; erse other than proton
  const int NoutCCQE =3;//[Nmu,Nproton,(Npion)]
	//****************** Detector Params. ********************
  const int kindsOfLayers = 4;
  const int nlayers_dummy = 1200; 
  const int nfibers_dummy = 1200;
	const float fiberSheetSize = 300; // square sheet
	const int nActualLayers = 50;
	//****************** Particle Gun Params. ********************
	const double myEnergyRange = 100.;

#endif