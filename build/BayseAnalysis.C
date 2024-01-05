#include<stdio.h>
#include <assert.h>
#include "../include/MyConst.hh"

void makemyBays(float);
void BayseAnalysis(){
	makemyBays(5.);
	makemyBays(10.);
	makemyBays(20.);
}

void makemyBays(float fpitch){
	TString filepath = "./sim_output/";
	TString input=Form("outputProtonHist_pitch%.0fmm",fpitch);
	TString output = Form("myBays_pitch%.0fmm",fpitch);

	//********* Input File *********
 	TFile *finput = new TFile(filepath+input+".root","READ");	
	const TH2D* hHitMomentum = static_cast<TH2D*>(finput->Get("hHitMomentum"));
	int nbinsMomentum=0,nbinsHit=0;
	nbinsMomentum = hHitMomentum->GetXaxis()->GetNbins();
	nbinsHit = hHitMomentum->GetYaxis()->GetNbins();
	double Pmax = hHitMomentum->GetXaxis()->GetBinUpEdge(hHitMomentum->GetNbinsX());

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
		if (nEntriesGivenPtrue >0 )
		{
			for (int jbinHit= 0; jbinHit < nbinsHit+1; ++jbinHit)
			{
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
		if(hHitMomentum->GetXaxis()->GetBinLowEdge(ibinHit)==0)continue;
		cout <<"nEntriesGivenNhit = "<<nEntriesGivenNhit<<endl;
		if (nEntriesGivenNhit >0)
		{
			for (int jbinMomentum = 1; jbinMomentum < nbinsMomentum+1; ++jbinMomentum){
				double entriesAtij=0;
				entriesAtij = hHitMomentum->GetBinContent(jbinMomentum,ibinHit);
				P_Preco_given_nhit->SetBinContent(ibinHit,jbinMomentum,entriesAtij/nEntriesGivenNhit);
			}
		}
		else
		{
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
	c1->Divide(2,1);
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

	TCanvas* c2 = new TCanvas("c2","c2",1000,600);
	c2->Divide(2,1);
	//******** Draw 2D ********
	c2->cd(1);
	P_Preco_given_Ptrue->SetTitle("P(Preco|Pture);Ptrue[MeV];Preco[MeV]");
	P_Preco_given_Ptrue->Draw("colz");

	//cout << "******** Draw Ptrue Slice  ********" <<endl;
	c2->cd(2);
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
	treeQR->Fill();

	// ********************** Calculate P(nhit|cos(theta)) ********************** 
	// MomentumDistributionCCQE = 
	// MomentumDistribution2p2h = 

	cout << "nbinsHit = "<<nbinsHit<<endl;
	cout << "Pmax = "<<Pmax<<endl;

	c1->Write();
	c2->Write();

	foutput->Write();
	//foutput->Close();


}