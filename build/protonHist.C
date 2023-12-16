/* Date: 23/11/18
 * pick up info which is necessary for machine learning.
 *
 *
*/
#include<iostream>
#include<stdio.h>
#include <assert.h>
#include "../include/MyConst.hh"

void protonHist(TString filepath = "./sim_output/"
		,TString input="mydata_protonStudy"
		,TString output = "outputProtonHist")
{

	const int nevent = 100000;// should be same or smaller than treeEvtAct entries
//Open {wavy}output file
 	TFile *finput = new TFile(filepath+input+".root","READ");	
	
	//read TTree
	//********* treeEvtAct *********
	TTree *treeEvtAct = static_cast<TTree *>(finput->Get("treeEvtAct1"));
	TVector3 *pointer_Pinitial =0;
	TVector3 *pointer_initialPosition=0;
	int ievt=0;
	double fmyEnergyRange=0;
	treeEvtAct->SetBranchAddress("ievt", &ievt);
	treeEvtAct->SetBranchAddress("Pinitial", &pointer_Pinitial);
	treeEvtAct->SetBranchAddress("myEnergyRange", &fmyEnergyRange);
	treeEvtAct->SetBranchAddress("initialPosition", &pointer_initialPosition);
	assert(treeEvtAct->GetEntries()>=nevent);
	//********* treeStpAct *********
  TTree *treeStpAct = (TTree*)finput->Get("treeStpAct");
	int evtFromStp=0,detid=0,trackid=0;
	treeStpAct->SetBranchAddress("evt", &evtFromStp);
	treeStpAct->SetBranchAddress("detid", &detid);
	treeStpAct->SetBranchAddress("trackid", &trackid);
	
	//********* protonTree *********
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
	double Pmax = 1800, cosmax=1., cosmin=-1.;
	double dnActualLayers = static_cast<double>(nActualLayers);

	TH2D* hHitMomentum = new TH2D("hHitMomentum","2D Hist: Momentum vs # of hit fibers"
				,Pbins,0.,Pmax
				,dnActualLayers,-0.5,dnActualLayers
				);
	TH2D* hHitCos = new TH2D("hHitCos","2D Hist: Cos vs # of hit fibers"
				,nbins_cos,cosmin,cosmax
				,dnActualLayers,-0.5,dnActualLayers
				);
	cout << "make 2Dhsitgram" <<endl;
	// Get nhit by event
	int nhitArray[nevent]={};
	size_t size_nhitArray= sizeof(nhitArray)/sizeof(nhitArray[0]);//"sizeof" method returns 
	for (int ientry = 0; ientry < (treeStpAct->GetEntries()); ++ientry)
	{
		treeStpAct->GetEntry(ientry);
		cout << "ientry = " << ientry <<endl;
		assert(size_nhitArray);
		if(evtFromStp >= nevent)break;//nhitArray[0],...,nhitArray[99]
		if (detid!=0 && trackid==1){
			nhitArray[evtFromStp] +=1;
		}
	}
	cout << "Fill nhitArray" <<endl;

	for(int ievent=0; ievent < nevent; ++ievent){
		//treeStpAct->GetEntry(ievent);
		treeEvtAct->GetEntry(ievent);
		TVector3 &Pinitial = *pointer_Pinitial; // jsut make a object to use Class.Member instead of Class->Member
		TVector3 &initialPosition = *pointer_initialPosition; // jsut make a object to use Class.Member instead of Class->Member

		std::cout <<"\n****************************"<<std::endl;
		std::cout <<"ievent="<< ievent <<std::endl;
		for (int idim= 0; idim < ndim; ++idim){
			std::cout<< Form("Pinitial[%d] = ",idim) << Pinitial[idim] <<std::endl;
			std::cout<< Form("initialPosition[%d] = ",idim) << initialPosition[idim] <<std::endl;
		} 
		nhit = nhitArray[ievent];//treeStpAct->GetEntries(Form("(evt==%d) &&(detid!=0)&&(trackid==1)",ievent));
		std::cout << "nhit= "		<< nhit << endl;
		cout <<"****************************"<<endl;
		hHitMomentum->Fill(Pinitial.Mag(),nhit);
		hHitCos->Fill(Pinitial.CosTheta(),nhit);
		outtree->Fill();
	}
	for (int i = 0; i < size_nhitArray; ++i)
	{
		cout << Form("nhitArray[%d] = ",i) << nhitArray[i] <<endl;
	}
	
	//c1->cd(1);
	hHitMomentum->SetTitle(Form("EnergyRange[0,%1f] event;Momentum [MeV]; # of hits",fmyEnergyRange));
	//hHitMomentum->Draw("colz");
	//c1->cd(2);
	hHitCos->SetTitle(Form("EnergyRange[0,%1f]; cos(theta); # of hits",fmyEnergyRange));
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
