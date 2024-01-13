/* Date: 23/11/18
 * pick up info which is necessary for machine learning.
 *
 *
*/
#include<iostream>
#include<stdio.h>
#include "../include/MyConst.hh"

#define DEBUG_ML 0

void ml(TString filepath = "./sim_output/MLinput/"
		,TString input="mydata_neut_5.4.0_675MeV_H2O_numu_1e6event"
		,TString pdfName= "figtest.pdf"
		,TString output = "ml_input"//"ml_input"
		, int NtotalEvents = 1e6
		, int neventsPerFile = 2e4
		){


for (int iEvent = 0; iEvent < NtotalEvents; iEvent++){
	if (iEvent % neventsPerFile ==0){

	//Open {wavy}output file
	TFile *finput = new TFile(filepath + input + ".root");
	if (!finput || finput->IsZombie()) {
	    std::cerr << "Error: Unable to open input file." << std::endl;
	    return;
	}
	const int nphase = 4;
	const int nlayers = 50;
	float chw[nlayers_dummy][nfibers_dummy];//energy deposit in each fibers.
	int NumOutParticle[NchargedGenerated];
	int eventFromEvtAct=0;

	//read TTree
	TTree *treeEvt = static_cast<TTree *>(finput->Get("treeEvtAct2"));
  TTree *treeMLInfo = (TTree*)finput->Get("treeInitialParticles");
  TTree *tree_yread = (TTree*)finput->Get("tree_yread");
	treeEvt->SetBranchAddress("chw", chw);

	treeMLInfo->SetBranchAddress("NumOutParticle", NumOutParticle);

	const int nfibers = tree_yread->GetEntries();

//****** set file# countor

	// make output file 
	TFile* foutput = new TFile(Form(filepath+output+"_file%d.root",iEvent/neventsPerFile),"recreate");
	TTree* outtree = new TTree("ml_input",Form("input data for Machine Leaning;evnt No. %d ~ %d",iEvent,iEvent+neventsPerFile-1));

//debug
#if DEBUG_ML
	//TH2D definition: just for viewing event
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	c1->Print(pdfName+"[","pdf");
	TH2D* view2d = new TH2D("view","view from x,u,y,v"
								,nlayers
								,0
								,nlayers
								,nfibers
								,0
								,nfibers
								);
#endif

// save dedx,Poisson,Poisson&Gauss as x[#of arrayed fiber][# of fiber sheets/phase][# of phase]
	const float MeV=1,g=1,cm=1,mm=0.1;
	const float dedxOfPositron= 2.53*MeV*g/std::pow(cm,2); //dE/dx of e+/e- in 675MeV/c
	const float rhoPolystylene = 1.05*g/pow(cm,3);
	const float fiberThickness= 0.96*mm;
	const float deTypical= dedxOfPositron*rhoPolystylene*fiberThickness;
	const float photonTypical = 29.4;//p.e/mm of 675MeV/c e+/e-
	const float energyToPhotonConversionFactor = photonTypical/deTypical;

	int jEventGlobal = 0;
	outtree->Branch("event",&jEventGlobal,"event/I");

	float dedx[nlayers][nfibers][nphase];
	float photon[nlayers][nfibers][nphase];//Poisson
	float electric[nlayers][nfibers][nphase];//Poisson&Gauss
	outtree->Branch("dedx",dedx,Form("dedx[%d][%d][%d]/F",nlayers,nfibers,nphase));
	outtree->Branch("photon",photon,Form("photon[%d][%d][%d]/F",nlayers,nfibers,nphase));
	outtree->Branch("electric",electric,Form("electric[%d][%d][%d]/F",nlayers,nfibers,nphase));

	for(int jEventLocal=0; jEventLocal < neventsPerFile; ++jEventLocal){
		jEventGlobal = jEventLocal + iEvent;
		treeEvt->GetEntry(jEventGlobal);
		cout <<"****************************"<<endl;
		cout <<"file No. ="<< iEvent/neventsPerFile <<endl;
		cout <<"jEventGlobal="<< jEventGlobal <<endl;
		cout <<"****************************"<<endl;

#if 1 //under  refactoring.
		for(int ilayer=0; ilayer < nlayers; ++ilayer){
			int iphase = ilayer%nphase; //0:x, 1:u, 2:y, 3:v, maybe...
			int	iphaselayer=ilayer/nphase;

			for(int ifiber=0; ifiber < nfibers; ++ifiber){
				const float edep = chw[ilayer][ifiber];
				dedx[iphaselayer][ifiber][iphase]=edep;
				const float pePoisson = gRandom->Poisson(edep*energyToPhotonConversionFactor);
				photon[iphaselayer][ifiber][iphase] = pePoisson;
				electric[iphaselayer][ifiber][iphase]=gRandom->Gaus(pePoisson);
#if DEBUG_ML
				if(iphase==0){
					view2d->Fill(ilayer,ifiber,electric[iphaselayer][ifiber][iphase]);
				}
#endif
			}

		}
#endif
		outtree->Fill();	
#if DEBUG_ML
		view2d->SetTitle(Form("ievent=%d phase=0;nlayer; nfiber",ievent));
		view2d->Draw("colz");
		c1->Update();
		c1->Print(pdfName,"pdf");
		view2d->Reset();
#endif
	}

#if DEBUG_ML
	c1->Print(pdfName+"]","pdf");
	delete c1;
#endif
//save initial particles momentum vector, if possible

//save as diffrent root file
	//view2d->Write();
	// outtree->Write();
	foutput->Write();
	foutput->Close();


	}
}
}

// //配列を引数にとるときに、安全性を確保するハードルが高かったため、一時保留。
// void calculatePhotonDistribution(int nphase,float dedx[][][], float chw[][][], float photon[][][], float electric[][][],  int nlayers, int nfibers){
	// const float MeV=1,g=1,cm=1,mm=0.1;
	// const float dedxOfPositron= 2.53*MeV*g/std::pow(cm,2); //dE/dx of e+/e- in 675MeV/c
	// const float rhoPolystylene = 1.05*g/pow(cm,3);
	// const float fiberThickness= 0.96*mm;
	// const float deTypical= dedxOfPositron*rhoPolystylene*fiberThickness;
	// const float photonTypical = 29.4;//p.e/mm of 675MeV/c e+/e-
	// const float energyToPhotonConversionFactor = photonTypical/deTypical;
// 
		// for(int ilayer=0; ilayer < nlayers; ++ilayer){
			// int iphase = ilayer%nphase; //0:x, 1:u, 2:y, 3:v, maybe...
			// int	iphaselayer=ilayer/nphase;
// 
			// for(int ifiber=0; ifiber < nfibers; ++ifiber){
				// const float edep = chw[ilayer][ifiber];
				// dedx[iphaselayer][ifiber][iphase]=edep;
				// const float pePoisson = gRandom->Poisson(edep*energyToPhotonConversionFactor);
				// photon[iphaselayer][ifiber][iphase] = pePoisson;
				// electric[iphaselayer][ifiber][iphase]=gRandom->Gaus(pePoisson);
			// }
		// }	
// }