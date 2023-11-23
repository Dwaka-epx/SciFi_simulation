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
	float dedx[nlayers][nfibers][nphase];

	
	//read TTree
	TTree *tree1 = static_cast<TTree *>(finput->Get("treeEvtAct2"));
  	TTree *tree2 = (TTree*)finput->Get("treeStpAct");
  	TTree *tree3 = (TTree*)finput->Get("tree_yread");
	tree1->SetBranchAddress("chw", chw);

	TFile* foutput = new TFile(filepath+output+".root","recreate");
	TTree* outtree = new TTree("ml_input","input data for ML");

#if 0
	int fevent=0, flayer=0,ffiber=0, fphase=0;
	float de=0;
	outtree->Branch("fevent",&fevent,"fevent/I");
	outtree->Branch("flayer",&flayer,"flayer/I");
	outtree->Branch("fphase",&fphase,"fphase/I");
	outtree->Branch("ffiber",&ffiber,"ffiber/I");
	outtree->Branch("dedx",&de,"dedx/F");
#endif

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
	


//make x[nlayers/phase=20][nfibers=200][#phase=4]

	
	for(int ievent=0; ievent < nevent; ++ievent){
		tree1->GetEntry(ievent);
		cout <<"****************************"<<endl;
		cout <<"ievent="<< ievent <<endl;
		cout <<"****************************"<<endl;
		//fevent=ievent;

		for(int ilayer=0; ilayer < nlayers; ++ilayer){
			int iphase = ilayer%nphase; //0:x, 1:u, 2:y, 3:v, maybe...
			int	iphaselayer=ilayer/nphase;
			//fphase = iphase;
			//flayer = iphaselayer;

			for(int ifiber=0; ifiber < nfibers; ++ifiber){
				float edep;
				edep = chw[ilayer][ifiber];
				dedx[iphaselayer][ifiber][iphase]=edep;
				cout << ievent 
				<< "," << ilayer 
				<< "," <<ifiber << endl;
				//ffiber=ifiber;
				//outtree->Fill();
				if(iphase==0){
					view2d->Fill(ilayer,ifiber,edep);
				}
			
			}

		}
		view2d->SetTitle(Form("ievent=%d;nlayer; nfiber",ievent));
		view2d->Draw("colz");
		c1->Update();
		c1->Print(pdfName,"pdf");
		view2d->Reset();
	}
	c1->Print(pdfName+""]","pdf");

//	
//	// 10th event
	//for(int iphase=0; iphase<nphase; ++iphase){
	//	cout <<"****************************"<<endl;
	//	cout <<"iphase="<< iphase <<endl;
	//	cout <<"****************************"<<endl;		
	//	for(int ifiber=0; ifiber < nfibers; ++ifiber){
	//		std::cout << "ifiber =" << ifiber;
	//		for(int ilayer=0; ilayer < nlayers; ++ilayer){
	//			//std::cout << chw[ilayer][ifiber]<<",";
	//			cout << dedx[ilayer][ifiber][iphase] << ",";
	//		}
	//		std::cout << std::endl;
	//	}
	//
	//
	//}
//
//// save dedx as x[#of arrayed fiber][# of fiber sheets/phase][# of phase]
//
//// save poisson as x[#of arrayed fiber][# of fiber sheets/phase][# of phase]
//
//// save poisson+gauss as x[#of arrayed fiber][# of fiber sheets/phase][# of phase]
//
//
////save initial particles ID track number? 
//
////save initial particles momentum
//
////save initial particles angle 
//
////save as diffrent root file
	//view2d->Write();
	foutput->Write();
	foutput->Close();
	
}
