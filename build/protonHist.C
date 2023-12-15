/* Date: 23/11/18
 * pick up info which is necessary for machine learning.
 *
 *
*/
#include<iostream>
#include<stdio.h>
#include "../include/MyConst.hh"

void protonHist(TString filepath = "./sim_output/"
		,TString input="mydata_protonStudy"
		,TString output = "outputProtonHist")
{

//Open {wavy}output file
 	TFile *finput = new TFile(filepath+input+".root","READ");	
	
	const int nevent = 10000;

	//read TTree
	TTree *treeEvtAct = static_cast<TTree *>(finput->Get("treeEvtAct1"));
  TTree *treeStpAct = (TTree*)finput->Get("treeStpAct");
	TVector3 *pointer_Pinitial =0;
	TVector3 *pointer_initialPosition=0;
	int ievt=0;
	double fmyEnergyRange=0;
	treeEvtAct->SetBranchAddress("ievt", &ievt);
	treeEvtAct->SetBranchAddress("Pinitial", &pointer_Pinitial);
	treeEvtAct->SetBranchAddress("myEnergyRange", &fmyEnergyRange);
	treeEvtAct->SetBranchAddress("initialPosition", &pointer_initialPosition);


	TFile* foutput = new TFile(filepath+output+".root","recreate");
	TTree* outtree = new TTree("protonTree","momentum, cos(theta), # fiber hits");
	int nhit =0;
	outtree->Branch("Pinitial", &pointer_Pinitial);
	outtree->Branch("initialPosition", &pointer_initialPosition);
	outtree->Branch("nhit", &nhit);


	//TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	//c1->Divide(2,2);
	int Pbins =100;
	int nbins_cos = 100;
	double Pmax = 1800;
	double dnActualLayers = static_cast<double>(nActualLayers);

	TH2D* hHitMomentum = new TH2D("hHitMomentum","2D Hist: Momentum vs # of hit fibers"
				,Pbins,0.,Pmax
				,dnActualLayers,-0.5,dnActualLayers
				);
	TH2D* hHitCos = new TH2D("hHitCos","2D Hist: Cos vs # of hit fibers"
				,nbins_cos,-1.,1.
				,dnActualLayers,-0.5,dnActualLayers
				);

	for(int ievent=0; ievent < nevent; ++ievent){
		treeStpAct->GetEntry(ievent);
		treeEvtAct->GetEntry(ievent);
		TVector3 &Pinitial = *pointer_Pinitial; // jsut make a object to use Class.Member instead of Class->Member
		TVector3 &initialPosition = *pointer_initialPosition; // jsut make a object to use Class.Member instead of Class->Member

		std::cout <<"\n****************************"<<std::endl;
		std::cout <<"ievent="<< ievent <<std::endl;
		for (int idim= 0; idim < ndim; ++idim){
			std::cout<< Form("Pinitial[%d] = ",idim) << Pinitial[idim] <<std::endl;
			std::cout<< Form("initialPosition[%d] = ",idim) << initialPosition[idim] <<std::endl;
		} 
		nhit = treeStpAct->GetEntries(Form("(evt==%d) &&(detid!=0)&&(trackid==1)",ievent));
		std::cout << "nhit= "		<< nhit << endl;
		cout <<"****************************"<<endl;
		hHitMomentum->Fill(Pinitial.Mag(),nhit);
		hHitCos->Fill(Pinitial.CosTheta(),nhit);
		outtree->Fill();

	}
	//c1->cd(1);
	hHitMomentum->SetTitle(Form("EnergyRange[0,%lf] event;Momentum [MeV]; # of hits",fmyEnergyRange));
	//hHitMomentum->Draw("colz");
	//c1->cd(2);
	hHitCos->SetTitle(Form("EnergyRange[0,%lf]; cos(theta); # of hits",fmyEnergyRange));
	//hHitCos->Draw("colz");

//save initial momentum

//save initial angle 

//save as diffrent root file
	//c1->Print(pdfName+"]", "pdf");
	//hHitMomentum->Write();
	//c1->Write();
	//outtree->Write();
	foutput->Write();
	foutput->Close();
	
}
