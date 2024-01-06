/* Date: 23/11/18
 * pick up info which is necessary for machine learning.
 * Draw 2D histgram of (Momentum,#hit) and (cos(theta),#hit) 
 *
*/
#include<iostream>
#include<stdio.h>
#include<algorithm>
#include <assert.h>
#include "../include/MyConst.hh"


#define USE_SECONDARY_HITS 1

void makeoutputProtonHist(float);
void protonHist(){
	makeoutputProtonHist(5.);
	// makeoutputProtonHist(10.);
	// makeoutputProtonHist(20.);
}

void makeoutputProtonHist(float fpitch){
	TString filepath = "./sim_output/";
	TString input	=	Form("mydata_protonStudy_pitch%.0fmm",fpitch);
	TString output = Form("outputProtonHist_pitch%.0fmm_SecondaryHits",fpitch);

	const int nevent = 100000;// should be same or smaller than treeEvtAct entries
//Open {wavy}output file
 	TFile *finput = new TFile(filepath+input+".root","READ");	
	
	//Read TTree
	//********* treeEvtAct *********
	TTree *treeEvtAct = static_cast<TTree *>(finput->Get("treeEvtAct1"));
	TVector3 *pointer_Pinitial =0;
	TVector3 *pointer_initialPosition=0;
	int ievt=0;
	double fmyMomentumRange=0;
	treeEvtAct->SetBranchAddress("ievt", &ievt);
	treeEvtAct->SetBranchAddress("Pinitial", &pointer_Pinitial);
	treeEvtAct->SetBranchAddress("myMomentumRange", &fmyMomentumRange);
	treeEvtAct->SetBranchAddress("initialPosition", &pointer_initialPosition);
	assert(treeEvtAct->GetEntries()>=nevent);
	//********* treeStpAct *********
  TTree *treeStpAct = (TTree*)finput->Get("treeStpAct");
	int evtFromStp=0,detid=0,trackid=0;
	float edep=0,edepCutoff=1e-3;
	treeStpAct->SetBranchAddress("evt", &evtFromStp);
	treeStpAct->SetBranchAddress("detid", &detid);
	treeStpAct->SetBranchAddress("trackid", &trackid);
	treeStpAct->SetBranchAddress("edep", &edep);
	
	//********* Output File *********
	TFile* foutput = new TFile(filepath+output+".root","recreate");
	TTree* outtree = new TTree("protonTree","momentum, cos(theta), # fiber hits");
	int nhit =0;
	outtree->Branch("Pinitial", &pointer_Pinitial);
	outtree->Branch("initialPosition", &pointer_initialPosition);
	outtree->Branch("nhit", &nhit);

	//************ Draw Histgram ************
	double myMargin = 0.12;
	gStyle->SetPadTopMargin(myMargin);
	gStyle->SetPadBottomMargin(myMargin);
	gStyle->SetPadLeftMargin(myMargin);
	gStyle->SetPadRightMargin(myMargin);
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	c1->Divide(2,1);
	int Pbins =100;
	int nbins_cos = 100;
	double Pmax = myMomentumRange, cosmax=1., cosmin=-1.;
	double dnActualLayers = static_cast<double>(detectorSizeZ/fpitch);

	TH2D* hHitMomentum = new TH2D("hHitMomentum","2D Hist: Momentum vs # of hit fibers"
				,Pbins,0.,Pmax
				,dnActualLayers,0,dnActualLayers
				);
	TH2D* hHitCos = new TH2D("hHitCos","2D Hist: Cos vs # of hit fibers"
				,nbins_cos,cosmin,cosmax
				,dnActualLayers,0,dnActualLayers
				);
	cout << "make 2Dhsitgram" <<endl;
	
	//********* Get nhit event by event *********
	int nhitArray[nevent]={};
	size_t size_nhitArray= sizeof(nhitArray)/sizeof(nhitArray[0]);//"sizeof" method returns size in unit of "byte".
	for (int ientry = 0; ientry < (treeStpAct->GetEntries()); ++ientry) //ientry corresponds istep.
	{
		treeStpAct->GetEntry(ientry);
		cout << "ientry = " << ientry <<endl;
		assert(size_nhitArray);
		if(evtFromStp >= nevent)break;//nhitArray[0],...,nhitArray[99]
#if (USE_SECONDARY_HITS ==0)
		if (detid!=0 && trackid==1){
			nhitArray[evtFromStp] +=1;
		}
#elif(USE_SECONDARY_HITS==1)
		if(detid!=0 && edep>edepCutoff){
			nhitArray[evtFromStp] +=1;
		}
#endif
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
	
	c1->cd(1);
	hHitMomentum->SetTitle(Form("Momentum Range[0,%.1f];Momentum [MeV]; # of hits",fmyMomentumRange));
	hHitMomentum->Draw("colz");
	c1->cd(2);
	hHitCos->SetTitle(Form("Momentum Range[0,%.1f]; cos(theta); # of hits",fmyMomentumRange));
	hHitCos->Draw("colz");

//save initial momentum

//save initial angle 


//save as diffrent root file
	//c1->Print(pdfName+"]", "pdf");
	//hHitMomentum->Write();
	c1->Write();
	//outtree->Write();
	foutput->Write();
	foutput->Close();
	cout << "\n ============ "<<output << " was made."<<" ============ \n"<<endl;
	
}
