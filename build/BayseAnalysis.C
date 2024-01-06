#include<stdio.h>
#include <assert.h>
#include <algorithm>

#include "../include/MyConst.hh"

void makemyBays(float fpitch);
void drawEfficiencySame();
void nhitDistribution(float fpitch, float lowerMomentum);
void momentumDistribution(float fpitch, int lowerHit);

void BayseAnalysis(){
// makemyBays(5.);
// makemyBays(10.);
// makemyBays(20.);

drawEfficiencySame();
}

void makemyBays(float fpitch){
	TString filepath = "./sim_output/";
	TString input=Form("outputProtonHist_pitch%.0fmm_onlyPrimaryProton",fpitch);
	TString output = Form("myBays_pitch%.0fmm_onlyPrimaryProton",fpitch);

	//********* Input File *********
 	TFile *finput = new TFile(filepath+input+".root","READ");	
	const TH2D* hHitMomentum = static_cast<TH2D*>(finput->Get("hHitMomentum"));
	int nbinsMomentum=0,nbinsHit=0;
	nbinsMomentum = hHitMomentum->GetNbinsX();
	nbinsHit = hHitMomentum->GetNbinsY();
	double Pmax = hHitMomentum->GetXaxis()->GetBinUpEdge(nbinsMomentum);

	//********* Output File *********
	TFile* foutput = new TFile(filepath+output+".root","recreate");
	//************ Draw Histgram ************
	TH2D* P_nhit_given_Ptrue = new TH2D("P_nhit_given_Ptrue"
		,"Probability of getting nhit given true momentum"
		,nbinsMomentum,0.,Pmax
		,nbinsHit,0,static_cast<double>(nbinsHit));// You only need normalization.
	TH2D* P_Preco_given_nhit = new TH2D("P_Preco_given_nhit"
		,"Probability of getting reco momentum given nhit"
		,nbinsHit,0,static_cast<double>(nbinsHit)
		,nbinsMomentum,0.,Pmax);// You have to transpose and normalize.

	//****************** P(nhit|P_true) ******************
	for (int ibinMomentum = 1; ibinMomentum < nbinsMomentum+1; ++ibinMomentum){//bin# starts from 1 to nbins
		cout << "ibinMomentum = " << ibinMomentum << endl;
		int nEntriesInBinx = 0;
		double nEntriesGivenPtrue = hHitMomentum->ProjectionY("histNhitGivenPtrue",ibinMomentum,ibinMomentum)->Integral();
		cout << "nEntriesGivenPtrue = " << nEntriesGivenPtrue << endl;
		if (nEntriesGivenPtrue >0 ){
			for (int jbinHit= 0; jbinHit < nbinsHit+1; ++jbinHit){
				double entriesAtij=0;				
				entriesAtij = hHitMomentum->GetBinContent(ibinMomentum,jbinHit);
				P_nhit_given_Ptrue->SetBinContent(ibinMomentum,jbinHit,entriesAtij/nEntriesGivenPtrue);
			}
		}else{
			continue;
		}

	}
	
	//****************** P(P_reco|nhit) ******************

	for (int ibinHit = 1; ibinHit < nbinsHit+1; ++ibinHit){//bin# starts from 1 to nbins
		cout << "ibinHit = " << ibinHit <<endl;
		double nEntriesGivenNhit = hHitMomentum->ProjectionX("histPrecoGivenNhit",ibinHit,ibinHit)->Integral();
		cout << "GetBinLowEdge(ibinHit) = " << hHitMomentum->GetXaxis()->GetBinLowEdge(ibinHit) <<endl;
		if(hHitMomentum->GetXaxis()->GetBinLowEdge(ibinHit)==0)continue; // need refactoring: this cut off should be done as follows: after iteration P(P_reco|nhit=0) =0.
		cout <<"nEntriesGivenNhit = "<<nEntriesGivenNhit<<endl;
		if (nEntriesGivenNhit >0){
			for (int jbinMomentum = 1; jbinMomentum < nbinsMomentum+1; ++jbinMomentum){
				double entriesAtij=0;
				entriesAtij = hHitMomentum->GetBinContent(jbinMomentum,ibinHit);
				P_Preco_given_nhit->SetBinContent(ibinHit,jbinMomentum,entriesAtij/nEntriesGivenNhit);
			}
		}
		else{
			continue;
		}
				
	}
	//************ Draw Probs ************
	double myMargin = 0.12;
	gStyle->SetPadTopMargin(myMargin);
	gStyle->SetPadBottomMargin(myMargin);
	gStyle->SetPadLeftMargin(myMargin);
	gStyle->SetPadRightMargin(myMargin);
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	c1->Divide(2,2);
	c1->cd(1);
	P_nhit_given_Ptrue->SetTitle("Probability of nhit given Ptrue;Momentum[MeV];nhit");
	P_nhit_given_Ptrue->Draw("colz");
	c1->cd(2);
	P_Preco_given_nhit->SetTitle("Probability of Preco given nhit;nhit;Momentum[MeV]");
	P_Preco_given_nhit->Draw("colz");
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	// ********************** Calculate P(Preco|Ptrue) **********************
	TH2D* P_Preco_given_Ptrue = new TH2D("P_Preco_given_Ptrue","Probability of Preco given Ptrue"
		,nbinsMomentum,0.,Pmax
		,nbinsMomentum,0.,Pmax);// You have to transpose and normalize.
	for (int ibinPtrue = 1; ibinPtrue < nbinsMomentum+1; ++ibinPtrue){
		for (int jbinPreco = 1; jbinPreco < nbinsMomentum+1; jbinPreco++){
			double entry_ij =0; 
			for (int kbinHit = 1; kbinHit < nbinsHit+1; ++kbinHit){
				entry_ij += (P_nhit_given_Ptrue->GetBinContent(ibinPtrue,kbinHit))*(P_Preco_given_nhit->GetBinContent(kbinHit,jbinPreco));
			}
			P_Preco_given_Ptrue->SetBinContent(ibinPtrue,jbinPreco,entry_ij);		
		}
	}

	//******** Draw 2D ********
	c1->cd(3);
	P_Preco_given_Ptrue->SetTitle("P(Preco|Pture);Ptrue[MeV];Preco[MeV]");
	P_Preco_given_Ptrue->Draw("colz");

	//cout << "******** Draw Ptrue Slice  ********" <<endl;
	c1->cd(4);
	float lowerMomentum =250.;
	float profileRange = 40.;// [lowerMomentum, lowerMomentum + profileRange]
	int lowerBin = P_Preco_given_Ptrue->GetXaxis()->FindBin(lowerMomentum);
	cout << "lowerMomentum= " << P_Preco_given_Ptrue->GetXaxis()->GetBinLowEdge(lowerBin) <<endl;
	int upperBin = P_Preco_given_Ptrue->GetXaxis()->FindBin(lowerMomentum+profileRange);
	cout << lowerBin <<" =< Bin =<" << upperBin << endl;
	float upperMomentum  = P_Preco_given_Ptrue->GetXaxis()->GetBinUpEdge(upperBin);
	cout << lowerMomentum<< "=< True Momentum <" << upperMomentum <<endl;
	TH1* p_Preco1d = P_Preco_given_Ptrue->ProjectionY("p_Preco1d",lowerBin,upperBin);
	p_Preco1d->SetTitle(Form(" (%.0f=<Pture<%.0f) pitch=%.0fmm;Preco [MeV];p(Preco|Ptrue)",lowerMomentum,upperMomentum,fpitch));
	p_Preco1d->SetStats(0);
	p_Preco1d->Draw();
	//******** Calculate 1sigma and 2sigma  ********
	TTree* treeQR = new TTree("treeQR" ,"quantile value of 68%, 95%, around 50%");
	const int nQuantile =5;
	double quantileValues[nQuantile]={};
	const double sigma1 = 0.34, sigma2 = 0.475;
	const Double_t quantiles[nQuantile]= {0.5-sigma2,0.5-sigma1,0.50,0.5+sigma1,0.5+sigma2};
	treeQR->Branch("quantileValues",quantileValues,Form("quantileValues[%d]/D",nQuantile));
	double quantileRange68=0, quantileRange95=0;
	treeQR->Branch("quantileRange68",&quantileRange68);//,"quantileRange68/D");
	treeQR->Branch("quantileRange95",&quantileRange95);//,"quantileRange95/D");
	p_Preco1d->GetQuantiles(nQuantile,quantileValues,quantiles);
	cout << "median = " << quantileValues[2] << endl; 
	cout << "quantile range 34%= " << Form("[%.2lf,%.2lf]",quantileValues[1],quantileValues[3]) << endl;
	cout << "quantile range 47.5%= " << Form("[%.2lf,%.2lf]",quantileValues[0],quantileValues[4]) << endl;
	quantileRange68 = (quantileValues[3]-quantileValues[1])/2.;
	quantileRange95 = (quantileValues[4]-quantileValues[0])/2.;
	cout << "quantileRange68 = " <<quantileRange68 <<endl;
	cout << "quantileRange95 = " <<quantileRange95 <<endl;

	// ********************** calcurate efficiency(n>=3) **********************
	const int nhitTrackableCutoff=3;
	const int nbinHitCutoff = hHitMomentum->GetYaxis()->FindBin(nhitTrackableCutoff);
	int nbinMaxHit = hHitMomentum->GetYaxis()->FindBin(nbinsHit);
	double momentumArray[nbinsMomentum], efficencyArray[nbinsMomentum] ;
	for (int ibinMomentum= 1; ibinMomentum <nbinsMomentum+1 ; ibinMomentum++){
		double efficiencyDenominator = hHitMomentum->ProjectionY("",ibinMomentum,ibinMomentum)->Integral();	
		double efficiencyNumerator = hHitMomentum->Integral(ibinMomentum,ibinMomentum,nbinHitCutoff,nbinMaxHit);	
		double efficiencyTrackable = efficiencyNumerator/efficiencyDenominator;

		momentumArray[ibinMomentum-1] = hHitMomentum->GetXaxis()->GetBinLowEdge(ibinMomentum);
		efficencyArray[ibinMomentum-1] = efficiencyTrackable;
	}
	treeQR->Branch("momentumArray",momentumArray,Form("momentumArray[%d]/D",nbinsMomentum));//,"quantileRange95/D");
	treeQR->Branch("efficencyArray",efficencyArray,Form("efficencyArray[%d]/D",nbinsMomentum));//,"quantileRange95/D");
	treeQR->Fill();

	TGraph* graphEfficiency = new TGraph(nbinsMomentum,momentumArray,efficencyArray);
	TCanvas* c2 = new TCanvas("c2","trackable efficiency",1000,600);
	graphEfficiency->SetTitle(Form("Trackable efficiency(nhit>=%d) pitch=%.0fmm; Momentum [MeV]; efficiency",nhitTrackableCutoff,fpitch));
	graphEfficiency->SetMarkerStyle(20);
	graphEfficiency->Draw("AP");
	

	// ********************** Calculate P(nhit|cos(theta)) ********************** 
	// MomentumDistributionCCQE = 
	// MomentumDistribution2p2h = 

	cout << "nbinsHit = "<<nbinsHit<<endl;
	cout << "Pmax = "<<Pmax<<endl;

	graphEfficiency->SetName("gr_trackable");
	graphEfficiency->Write();
	c1->Write();
	c2->Write();
	foutput->Write();
	foutput->Close();

}


//================================================================
void drawEfficiencySame(){
	TString filepath = "./sim_output/";
 	TFile *finput0 = new TFile(filepath+"myBays_pitch5mm_onlyPrimaryProton.root","READ");	
	TGraph *gr0 = static_cast<TGraph*>(finput0->Get("gr_trackable"));
	TTree *treeQR0 = static_cast<TTree *>(finput0->Get("treeQR"));
	TH2D* P_Preco_given_Ptrue = static_cast<TH2D*>(finput0->Get("P_Preco_given_Ptrue"));
	const int nbinsMomentum = P_Preco_given_Ptrue->GetNbinsX();
	double momentumArray[nbinsMomentum], efficencyArray0[nbinsMomentum];
	treeQR0->SetBranchAddress("momentumArray",momentumArray);
	treeQR0->SetBranchAddress("efficencyArray",efficencyArray0);
	treeQR0->GetEntry(0);

 	TFile *finput1 = new TFile(filepath+"myBays_pitch10mm_onlyPrimaryProton.root","READ");	
	TGraph *gr1 = static_cast<TGraph*>(finput1->Get("gr_trackable"));
	TTree *treeQR1 = static_cast<TTree *>(finput1->Get("treeQR"));
	double efficencyArray1[nbinsMomentum];
	treeQR1->SetBranchAddress("efficencyArray",efficencyArray1);
	treeQR1->GetEntry(0);

 	TFile *finput2 = new TFile(filepath+"myBays_pitch20mm_onlyPrimaryProton.root","READ");	
	TGraph *gr2 = static_cast<TGraph*>(finput2->Get("gr_trackable"));
	TTree *treeQR2 = static_cast<TTree *>(finput2->Get("treeQR"));
	double efficencyArray2[nbinsMomentum];
	treeQR2->SetBranchAddress("efficencyArray",efficencyArray2);
	treeQR2->GetEntry(0);

	TCanvas* c1 = new TCanvas("c1","trackable efficiency",1000,600);
	c1->Divide(2,1);
	c1->cd(1);
	gr1->SetMarkerStyle(20);
	gr2->SetMarkerStyle(20);
	gr1->SetMarkerColor(kRed);
	gr2->SetMarkerColor(kBlue);
	gr0->SetTitle("Trackable Efficiency (nhit>=3)");
	gr0->Draw("AP");
	gr1->Draw("P");
	gr2->Draw("P");
	TLegend *leg1 = new TLegend(0.65,0.13,0.85,0.35);
	leg1->AddEntry(gr0,"5mm","p");
	leg1->AddEntry(gr1,"10mm","p");
	leg1->AddEntry(gr2,"20mm","p");
	leg1->Draw();

	c1->cd(2);
	double max5over10=0,max10over20=0 ;
	double efficencyDif5over10[nbinsMomentum],efficencyDif10over20[nbinsMomentum];
	for (int ibinMomentum= 0; ibinMomentum < nbinsMomentum; ++ibinMomentum){
		if(efficencyArray1[ibinMomentum]>0){
			efficencyDif5over10[ibinMomentum] = efficencyArray0[ibinMomentum]/efficencyArray1[ibinMomentum];
			if(efficencyDif5over10[ibinMomentum] > max5over10){ max5over10=efficencyDif5over10[ibinMomentum];}
		}
		if(efficencyArray2[ibinMomentum>0]){
			efficencyDif10over20[ibinMomentum] = efficencyArray1[ibinMomentum]/efficencyArray2[ibinMomentum];
			if(efficencyDif10over20[ibinMomentum] > max10over20){ max10over20=efficencyDif10over20[ibinMomentum];}
		}
	}
	TGraph* graphEfficiencyRatio5over10 = new TGraph(nbinsMomentum,momentumArray,efficencyDif5over10);
	TGraph* graphEfficiencyRatio10over20 = new TGraph(nbinsMomentum,momentumArray,efficencyDif10over20);
	graphEfficiencyRatio5over10->SetTitle("effciency ratio; Momentum [MeV]; ratio");
	graphEfficiencyRatio5over10->SetMarkerStyle(20);
	graphEfficiencyRatio10over20->SetMarkerStyle(20);
	graphEfficiencyRatio5over10->SetMarkerColor(kRed);
	graphEfficiencyRatio10over20->SetMarkerColor(kBlue);
	graphEfficiencyRatio5over10->Draw("AP");
	graphEfficiencyRatio10over20->Draw("P");

	TLegend *leg2 = new TLegend(0.45,0.75,0.9,0.9);
	leg2->AddEntry(graphEfficiencyRatio5over10,"5mm/10mm","p");
	leg2->AddEntry(graphEfficiencyRatio10over20,"10mm/20mm","p");
	leg2->Draw();

	cout << "max5over10 = " << max5over10 <<endl;
	cout << "max10over20 = " << max10over20 <<endl;
	

}

//================================================================

void nhitDistribution(float fpitch, float lowerMomentum){
		TString filepath = "./sim_output/";
	TString input=Form("outputProtonHist_pitch%.0fmm_onlyPrimaryProton",fpitch);
 	TFile *finput = new TFile(filepath+input+".root","READ");
	const TH2D* hHitMomentum = static_cast<TH2D*>(finput->Get("hHitMomentum"));
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);

	float profileRange = 40.;// [lowerMomentum, lowerMomentum + profileRange]
	int lowerBin = hHitMomentum->GetXaxis()->FindBin(lowerMomentum);
	cout << "lowerMomentum= " << hHitMomentum->GetXaxis()->GetBinLowEdge(lowerBin) <<endl;
	int upperBin = hHitMomentum->GetXaxis()->FindBin(lowerMomentum+profileRange);
	cout << lowerBin <<" =< Bin =<" << upperBin << endl;
	float upperMomentum  = hHitMomentum->GetXaxis()->GetBinUpEdge(upperBin);
	cout << lowerMomentum<< "=< True Momentum <" << upperMomentum <<endl;
	TH1* p_nhit1d = hHitMomentum->ProjectionY("p_nhit1d",lowerBin,upperBin);
	p_nhit1d->SetTitle(Form("p(nhit|Ptrue) (%.0f=<Pture<%.0f) pitch=%.0fmm; nhit ; event",lowerMomentum,upperMomentum,fpitch));
	p_nhit1d->SetStats(0);
	p_nhit1d->Draw();
}

//================================================================


void momentumDistribution(float fpitch, int lowerHit){
		TString filepath = "./sim_output/";
	TString input=Form("outputProtonHist_pitch%.0fmm_onlyPrimaryProton",fpitch);
 	TFile *finput = new TFile(filepath+input+".root","READ");
	const TH2D* hHitMomentum = static_cast<TH2D*>(finput->Get("hHitMomentum"));
	TCanvas* c1 = new TCanvas("c1","c1",1000,600);
	int profileRange = 0.;
	int lowerBin = hHitMomentum->GetYaxis()->FindBin(lowerHit);
	cout << "lowerMomentum= " << hHitMomentum->GetXaxis()->GetBinLowEdge(lowerBin) <<endl;
	int upperBin = hHitMomentum->GetYaxis()->FindBin(lowerHit+profileRange);
	cout << lowerBin <<" =< Bin =<" << upperBin << endl;
	int upperHit  = hHitMomentum->GetYaxis()->GetBinUpEdge(upperBin);
	cout << lowerHit<< "=< nhit <" << upperHit <<endl;
	TH1* p_Preco1d = hHitMomentum->ProjectionX("p_Preco1d",lowerBin,upperBin);
	p_Preco1d->SetTitle(Form("p(Preco|nhit) (%d=<nhit<%d) pitch=%.0fmm; Preco [MeV] ; event",lowerHit,upperHit,fpitch));
	p_Preco1d->SetStats(0);
	p_Preco1d->Draw();
}