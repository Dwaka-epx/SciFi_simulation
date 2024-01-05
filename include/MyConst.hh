
#ifndef MyConst_h

#define MyConst_h 1
	//****************** General ********************
  const int ndim =3;
	//****************** NEUT ********************
	#define USE_NEUT 0 // flag for NEUT STUDY or PROTON STUDY
	const	int nmaxParticles = 20;// Buffer number of out particle in neutrino interaction: should be larger enough than # of out particle(typically qround 10)
  const bool flagProtonRejection = false; // flag for ML; erase proton.
  const bool flagOnlyProton = false; // flag for ML; erse other than proton
  const int NoutCCQE =3;//[Nmu,Nproton,(Npion)]
	//****************** Detector Params. ********************
  const int kindsOfLayers = 4;
  const int nlayers_dummy = 1200; 
  const int nfibers_dummy = 1200;
	const float fiberSheetSize = 300; // square sheet
	const float layersPitch = 10.; //pitch*nActualLayers=500
	const double detectorSizeZ = 500.;//mm
	const int nActualLayers = static_cast<int>(detectorSizeZ/layersPitch);
	//****************** Particle Gun Params. ********************
	const double myMomentumRange = 1000.;//
	const double myProtonMass = 938.27208816;
	//****************** Drawing Histgrams ********************
	//const double nHitCutOff = 
#endif