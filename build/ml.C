/* Date: 23/11/18
 * pick up info which is necessary for machine learning.
 *
 *
*/
#include<iostream>
#include<stdio.h>
#include "../include/sizeOfFiberArray.hh"

void ml(TString filepath = "./sim_output/"
		,TString input="mydata_neut_5.4.0_675MeV_H2O_numu_1e5event"
		,TString pdfName= "figtest.pdf"
		,TString output = "ml_input"
		){

//Open {wavy}output file
 	TFile *finput = new TFile(filepath+input+".root");
	
	const int nevent = 10;
	const int nfibers = 300;
	const int nphase = 4;
	const int nlayers = 12*nphase ;
	float chw[nlayers_dummy][nfibers_dummy];//energy deposit in each fibers.

	//read TTree
	TTree *tree1 = static_cast<TTree *>(finput->Get("treeEvtAct2"));
  	TTree *tree2 = (TTree*)finput->Get("treeStpAct");
  	TTree *tree3 = (TTree*)finput->Get("tree_yread");
	tree1->SetBranchAddress("chw", chw);

	TFile* foutput = new TFile(filepath+output+".root","recreate");
	TTree* outtree = new TTree("ml_input","input data for ML");

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
	

// save dedx,Poisson,Poisson&Gauss as x[#of arrayed fiber][# of fiber sheets/phase][# of phase]

	float dedx[nlayers][nfibers][nphase];
	float photon[nlayers][nfibers][nphase];//Poisson
	float electric[nlayers][nfibers][nphase];//Poisson&Gauss
	outtree->Branch("dedx",dedx,Form("dedx[%d][%d][%d]/F",nlayers,nfibers,nphase));
	outtree->Branch("photon",photon,Form("dedx[%d][%d][%d]/F",nlayers,nfibers,nphase));
	outtree->Branch("electric",electric,Form("dedx[%d][%d][%d]/F",nlayers,nfibers,nphase));
	
	float MeV=1,g=1,cm=1,mm=0.1;
	float dedxOfPositron= 2.53*MeV*g/pow(cm,2); //dE/dx of e+/e- in 675MeV/c
	float rhoPolystylene = 1.05*g/pow(cm,3);
	float fiberThickness= 0.96*mm;
	float deTypical= dedxOfPositron*rhoPolystylene*fiberThickness;
	float photonTypical = 29.4;//p.e/mm of 675MeV/c e+/e-
	float energyToPhotonConversionFactor = photonTypical/deTypical;
	float pePoisson=0;

	for(int ievent=0; ievent < nevent; ++ievent){
		tree1->GetEntry(ievent);
		cout <<"****************************"<<endl;
		cout <<"ievent="<< ievent <<endl;
		cout <<"****************************"<<endl;

		for(int ilayer=0; ilayer < nlayers; ++ilayer){
			int iphase = ilayer%nphase; //0:x, 1:u, 2:y, 3:v, maybe...
			int	iphaselayer=ilayer/nphase;

			for(int ifiber=0; ifiber < nfibers; ++ifiber){
				float edep;
				edep = chw[ilayer][ifiber];
				dedx[iphaselayer][ifiber][iphase]=edep;
				pePoisson = photon[iphaselayer][ifiber][iphase]=gRandom->Poisson(edep*energyToPhotonConversionFactor);
				electric[iphaselayer][ifiber][iphase]=gRandom->Gaus(pePoisson);
				if(iphase==0){
					view2d->Fill(ilayer,ifiber,edep);
				}
			
			}

		}
		outtree->Fill();
		//view2d->SetTitle(Form("ievent=%d;nlayer; nfiber",ievent));
		//view2d->Draw("colz");
		//c1->Update();
		//c1->Print(pdfName,"pdf");
		//view2d->Reset();
	}
	c1->Print(pdfName+"]","pdf");



//save initial particles ID track number? track length


//save initial particles momentum

//save initial particles angle 

//save as diffrent root file
	//view2d->Write();
	outtree->Write();
	foutput->Write();
	foutput->Close();
	
}
