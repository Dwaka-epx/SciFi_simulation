/*


*/

TH1F *dedx;
TH1F *dedxXZ;


void SetTLegend(TLegend *leg, TGraph* gra);
void DrawTrueStepping(TCanvas *c1, TTree *tree2, int i);

void event2(TString fname="mydata_numu_600MeV_CCQE_Water")
{
   gStyle->SetPadRightMargin(0.12);
	gStyle->SetPalette(kCool);
	gStyle->SetOptStat(1);
	gStyle->SetLineScalePS(1); // to make narrow line for marker
   gStyle->SetTitleSize(0.05,"xyz");

	TFile *tf = new TFile(fname+".root");
	TString canName = "fig_"+fname+".pdf";
	tf->Print();// << endl;
	cerr << "canName = " << canName << endl;

	//
	// dedx in fiber in layer
	//
	int nevt;
   const float fiber_thickness = 1;
   const int unit = 32; // x.y.sx30 = 32
   const int actuallayers = unit * 20;
   const int nlayers_dummy = 1200; 
   const int nfibers_dummy = 1200;
   float chw[nlayers_dummy][nfibers_dummy];//={};
	
   TTree *tree1 = static_cast<TTree *>(tf->Get("treeEvtAct2"));
   tree1->SetBranchAddress("nevt", &nevt);
   //tree1->SetBranchAddress("chx", chx);
   //tree1->SetBranchAddress("chy", chy);
   tree1->SetBranchAddress("chw", chw);

	//
//	true trajectry tepping point
	//
   TTree *tree2 = (TTree*)tf->Get("treeStpAct");

	//
	// readout position
	//
   TTree *tree3 = (TTree*)tf->Get("tree_yread");

	int id;
   TVector3 *v3;// = new TVector3();
   tree3->SetBranchAddress("v3",&v3);
   tree3->SetBranchAddress("id",&id);

   int nReadOut = tree3->GetEntries();
   std::vector<TVector3> arrXZv3;
   std::vector<TVector3> arrYZv3;

   TVector3 zaxis(0,0,1);

   for (int i = 0; i < nReadOut; ++i) {
      tree3->GetEntry(i);
		//if ( !v3 ) continue;
      cerr << i << " id = " << id << " "; v3->Print();
      //arrXZv3[id] = *v3;
      arrXZv3.push_back(*v3);
      // rotation around z 
      // mother volume center is 0, and it is rotated around z
      v3->Rotate(TMath::Pi()/2., zaxis);
      arrYZv3.push_back(*v3);
   }

   // 
   // readout position for sparce layer
   // readout position for wavy layer 
   //
   TTree *tree4 = (TTree*)tf->Get("tree_fread");
   tree4->SetBranchAddress("v3",&v3);
   tree4->SetBranchAddress("id",&id);

   int nReadOut_sparce = tree4->GetEntries();
   std::vector<TVector3> arrYZv3_sparce;
   //TVector3 arrYZv3_sparce[999];

   for (int i = 0; i < nReadOut_sparce; ++i) {
      tree4->GetEntry(i);
      cerr << i << " id = " << id << " "; v3->Print();
      arrYZv3_sparce.push_back(*v3);
      //arrYZv3_sparce[ id ] = *v3;
   }
	


   int colN = TColor::GetNumberOfColors();
   cerr << "color max = " << colN << endl;
   // 0.3 is MIP position
   float maxdEdx1  = 0.50;
   double henkan = 40.0/0.138*0.707106781;
   double maxphoton = maxdEdx1*henkan;
   float colorStep1= maxdEdx1/colN;
   

	//
   // just dummy to get palette
	//
   TH2F *colpale = new TH2F("colpale","colpale",100,-1,1,100,-1,1);
   colpale->SetBinContent(0,0,maxdEdx1);
   colpale->GetZaxis()->SetRangeUser(0,maxphoton);
   colpale->GetZaxis()->SetTitle("p.e.");
   colpale->GetZaxis()->SetTitleOffset(0.55);
   colpale->GetZaxis()->SetLabelSize(0.030);
   colpale->GetZaxis()->SetLabelOffset(0.005);
   TCanvas *c2 = new TCanvas();
   colpale->Draw("colz");
   gPad->Update();
   TPaletteAxis *palette = (TPaletteAxis*)colpale->GetListOfFunctions()->FindObject("palette");

   double total_bin = 1000.;
   double max_bin = 3.;
   dedx = new TH1F("dedx","dedx",total_bin,0,max_bin);
   dedxXZ = new TH1F("dedxXZ","dedxXZ",total_bin,0,max_bin);
  
   TGraph *grXZ = new TGraph(); //readout point 
   TGraph *grYZ = new TGraph(); //readout point 
   
   grYZ->SetMarkerColor(kGray+1);
   grXZ->SetMarkerColor(kGray+1);
   grXZ->GetXaxis()->SetLimits(-20, fiber_thickness * actuallayers +20); // ok @21/09/14
   grYZ->GetXaxis()->SetLimits(-20, fiber_thickness * actuallayers +20); // ok @21/09/14

   grYZ->SetMarkerSize(8);// 21/12/20 

   TCanvas *c1 = new TCanvas("c1","c1",1000,600);
   c1->Divide(2,2);
        c1->cd(1); grYZ->SetTitle(";z axis [mm]; y axis [mm]");
        c1->cd(3); grYZ->SetTitle(";z axis [mm]; y axis [mm]");
       	c1->cd(2); grXZ->SetTitle(";z axis [mm]; x axis [mm]");
	c1->cd(4); grXZ->SetTitle(";z axis [mm]; x axis [mm]");
	c1->Print(canName+"[", "pdf");

   TLegend *leg = new TLegend(0.70,0.15,0.83,0.40);
   SetTLegend(leg, grXZ);

       
   //
   // readout point in x/y same layer
   //
   for (int l = 0; l < actuallayers; l++) {
   	for (int i = 0; i < arrXZv3.size(); ++i) { 
     	   const int layersX = 5; // ampli 50 mm, y 30 cm 
         int z_mother_position = fiber_thickness /* fiber thickness */ * l; 

         int wavyLayerZID;
         int wavyLayerXID;
         if ( i%2==0) { // even
            wavyLayerZID = i / ((layersX+1)*2);
            wavyLayerXID = i % ((layersX+1)*2);
         }
         if ( i%2==1) { // odd
            wavyLayerZID = i / ((layersX)*2);
            wavyLayerXID = i % ((layersX)*2);
         }

         if ( l%unit==0 ) grXZ->SetPoint(grXZ->GetN(), arrXZv3[i].Z() + z_mother_position, 
                                                       arrXZv3[i].X() /* y-axis on a plot is x poisitoon of the fiber */ );
         if ( l%unit==1 ) grYZ->SetPoint(grYZ->GetN(), arrYZv3[i].Z() + z_mother_position, 
                                                       arrYZv3[i].Y() /* y-axis on a plot is y poisitoon of the fiber */ );

	}
   }


  
   //
   // readout point in sparce layer
   // readout point in wavy layer
   //
   for (int l = 0; l < actuallayers; l++) {
      //for (int i = 0; i < arrYZv3_sparce.size(); ++i) { 
      for (int i = 0; i < nReadOut_sparce; ++i) { 

      int z_mother_position = fiber_thickness /* fiber thickness */ * l;
      if      ( l%unit==0 ) /* nothing */;
      else if ( l%unit==1 ) /* nothing */;	

      else  grYZ->SetPoint(grYZ->GetN(), arrYZv3_sparce[i].Z() + z_mother_position,
                                         arrYZv3_sparce[i].X() /* y-axis on a plot is y poisitoon of the fiber */ );
      }
   }


   TRandom3 *rnd = new TRandom3();
   //
   // event loop
   // 
     int ndisplay =0;
     for (int i=0; i<30; i++) {

       if(ndisplay>=20) break;
	tree1 ->GetEntry(i);
	#if 0 
	cerr << "nevt = " << nevt << endl;
	#endif
	
	    //2022/02/07
	    #if 0
	    int c;
	    cout << "Enter any character."<<endl;
	    c = getchar();
            #endif

	double totedep=0;
	int nhit=0;
	for (int l = 0; l < actuallayers; l++){ 
	  for (int ID= 0; ID< nfibers_dummy; ID++){ 
	    totedep += chw[l][ID];
	    if ( chw[l][ID]>1E-10 ) nhit++;
	  }
	}
		
		cout << "totedep " << totedep << " nhit= " << nhit << endl;
		if(nhit < 10.) continue;
		ndisplay++;
		cout<<"ndisplay ="<< ndisplay<<endl;


   	c1->cd(1); grYZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(3); grYZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(2); grXZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(4); grXZ->Draw("ap"); palette->Draw(); // ok @21/09/14 
		//return;

	       
   	//for (int l = 0; l < 30; l++) {
   	for (int l = 0; l < actuallayers; l++) {

      	// const int nfibers_dummy = 699; in WLSEventAction.hh
         // if exceeding it, you get strange value ?
	  for (int ID= 0; ID< nfibers_dummy; ID++) { 

   	    	if ( chw[l][ID]<1E-10 ) continue; // no edep
				dedx->Fill(chw[l][ID]);
            //cerr << "layer = " << l << " ID = " << ID << " dedx = " << chw[l][ID] << endl;

				// ID is even or odd
   			const int layersX = 5; // ampli 50 mm, y 30 cm 
				int wavyLayerZID, wavyLayerXID;
         	if ( ID%2==0) { // even
	         	wavyLayerZID = ID / ((layersX+1)*2); 
   	         wavyLayerXID = ID % ((layersX+1)*2);
      	   }
         	if ( ID%2==1) { // odd
	            wavyLayerZID = ID / ((layersX)*2); 
   	         wavyLayerXID = ID % ((layersX)*2);
      	   }
				#if 0
   	      cerr << "mother = " << l << " id = " << ID 
      	        << " wavyLayerZID = " << wavyLayerZID 
         	     << " wavyLayerXID = " << wavyLayerXID 
            	  << " edep = " << chw[l][ID] ; arrXZv3[ID].Print();
	              //<< endl;
				#endif
            int z_mother_position = fiber_thickness /* fiber thickness */ * l;


         
            //
            // x/y fiber layer, and sparce fiber layer
            //

	    double marker_size = 1;//0.4
	    double photon0 = chw[l][ID]*henkan;
	    double photon = rnd->Poisson(photon0);
	    float colorStep1= maxdEdx1/colN*henkan; // @21/12/24 maxdEdx1/colN	   
	    int min_photon = 2;

	    
            if      ( l%unit==0 ) { //
                    TMarker* mk = new TMarker( arrXZv3[ID].Z() + z_mother_position, 
                                               arrXZv3[ID].X(), 
                                               21);
		    mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( photon / colorStep1 )   ) );
                   mk->SetMarkerSize(marker_size);
		   mk->SetMarkerStyle(24); // 21/12/20
                   c1->cd(2) /* XZ */; 
		   if (photon >= min_photon) mk->Draw();
		   c1->cd(4) /* XZ */; //mk->Draw();
		   if (photon >= min_photon){
		     mk->Draw();
		     //DrawTrueStepping(c1, tree2, i);
		   }
		   
		   dedxXZ->Fill(chw[l][ID]); 				
	    }
            else if ( l%unit==1 ) { //
                   TMarker* mk = new TMarker( arrYZv3[ID].Z() + z_mother_position, 
                                             -arrYZv3[ID].Y(), 
                                              21);
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( photon / colorStep1 )  ) );
                   mk->SetMarkerSize(marker_size);
		   mk->SetMarkerStyle(24); // 21/12/2x0
                   c1->cd(1)/* YZ */; 
		   if (photon >= min_photon) mk->Draw();
		   c1->cd(3)/* YZ */; //mk->Draw();
   		   if (photon >= min_photon){
		     mk->Draw();
		     // DrawTrueStepping(c1, tree2, i);
		   }
 	    }
            else { 
                   TMarker* mk = new TMarker( arrYZv3_sparce[ID].Z() + z_mother_position,
                                             -arrYZv3_sparce[ID].X(), 
                                              21);
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( photon / colorStep1 )  ) );
                   mk->SetMarkerSize(marker_size);
		   mk->SetMarkerStyle(24); // 21/12/20
                   c1->cd(1)/* YZ */; 
		   if (photon >= min_photon) mk->Draw();
		   c1->cd(3)/* YZ */; //mk->Draw();
		   if (photon >= min_photon){
		     mk->Draw();
		     //DrawTrueStepping(c1, tree2, i);
		   }
	    }
	  }
	}
        
      
	DrawTrueStepping(c1, tree2, i);

        #if 1 
	c1->cd(3);leg->Draw();
	c1->cd(4);leg->Draw();
	c1->Print(canName, "pdf");
        #endif
	c1->Update();
	    //2022/02/07

     }
     
   c1->Print(canName+"]", "pdf");
   
   #if 0
   TCanvas *c0 = new TCanvas("c0","c0",600,600);
   c0->Divide(1,2);
   c0->cd(1); dedx->SetTitle(";Energy deposit [MeV]; Counts");
   c0->cd(2); dedxXZ->SetTitle(";Energy deposit [MeV]; Counts");
   c0->cd(1); dedx->Draw();
   c0->cd(2); dedxXZ->Draw();
   
	cout<<"hist sum "<< dedx->GetSum()<<endl;
	int  maxbin = dedx->GetMaximumBin();
	double peak_deposit = maxbin*max_bin/total_bin;
	cout<<"Peak energy deposit "<< peak_deposit <<endl;
   #endif
}


#define MU 13
#define ELEC 11
#define GAMMA 22
#define PROTON 2212
#define NEUTRON 2112
#define POSIT -11

#define MUCOL kGray+1
#define ELECCOL kBlue+2
#define POSITCOL kRed+2
#define GAMMACOL kGreen+2
#define PROTONCOL kYellow+1
#define NEUTRONCOL kPink+1
#define PIPLUSCOL kMagenta
#define PIMINUSCOL kBlue-10
#define PI0COL kSpring-2


void DrawTrueStepping(TCanvas *c1, TTree *tree2, int i)
{
      const int volLimX = 10000;//detwidth;
      const int volLimY = 10000;//detwidth;
      stringstream term0, term1, term2, term3, term4;
      term0 << " detid!=0 &&evt==" << i << ends;
      term1 << " detid!=0 && code==+"<< MU << "&&evt==" << i << ends;
      term2 << " detid!=0 && code=="<< ELEC <<"&&evt==" << i << ends;
      term3 << " detid!=0 && code=="<< GAMMA <<"&&evt==" << i << ends;
      term4 << " detid!=0 && code=="<<  POSIT <<"&&evt==" << i << ends; 
      stringstream termP;
      stringstream termN;
      termP << " detid!=0 && code==2212 &&evt==" << i << ends;
      termN << " detid!=0 && code==2112 &&evt==" << i << ends;
      stringstream termPIPLUS, termPIMINUS, termPI0;
      termPIPLUS << " detid!=0 && code==+211 &&evt==" << i << ends;
      termPIMINUS << " detid!=0 && code==-211 &&evt==" << i << ends;
      termPI0 << " detid!=0 && code==111 &&evt==" << i << ends;
      
      stringstream termA, termB;
      termA << "  y:z-1 " << ends;
      termB << "  x:z-1 " << ends;

      #if 1
      tree2->SetMarkerStyle(21);
      tree2->SetMarkerSize(0.6);
      tree2->SetLineWidth(1);

      c1->cd(3);
      tree2->SetMarkerColor(1); tree2->Draw(termA.str().data(),term0.str().data(),"same");
      tree2->SetMarkerColor(GAMMACOL);    tree2->Draw(termA.str().data(), term3.str().data(),"same"); //gamma
      tree2->SetMarkerColor(NEUTRONCOL);     tree2->Draw(termA.str().data(), termN.str().data(),"same");//neutron
      tree2->SetMarkerColor(MUCOL);     tree2->Draw(termA.str().data(), term1.str().data(),"same");
      tree2->SetMarkerColor(ELECCOL);     tree2->Draw(termA.str().data(), term2.str().data(),"same");//electron
      tree2->SetMarkerColor(POSITCOL);      tree2->Draw(termA.str().data(), term4.str().data(),"same");
      tree2->SetMarkerColor(PROTONCOL);   tree2->Draw(termA.str().data(), termP.str().data(),"same");
      tree2->SetMarkerColor(PIPLUSCOL);     tree2->Draw(termA.str().data(), termPIPLUS.str().data(),"same");
      tree2->SetMarkerColor(PIMINUSCOL);     tree2->Draw(termA.str().data(), termPIMINUS.str().data(),"same");
      tree2->SetMarkerColor(PI0COL);     tree2->Draw(termA.str().data(), termPI0.str().data(),"same");
      
      c1->cd(4);
      tree2->SetMarkerColor(1); tree2->Draw(termB.str().data(),term0.str().data(),"same");
      tree2->SetMarkerColor(GAMMACOL);    tree2->Draw(termB.str().data(), term3.str().data(),"same"); //gamma
      tree2->SetMarkerColor(NEUTRONCOL);     tree2->Draw(termB.str().data(), termN.str().data(),"same");//neutron
      tree2->SetMarkerColor(MUCOL);     tree2->Draw(termB.str().data(), term1.str().data(),"same");
      tree2->SetMarkerColor(ELECCOL);     tree2->Draw(termB.str().data(), term2.str().data(),"same");//electron
      tree2->SetMarkerColor(POSITCOL);      tree2->Draw(termB.str().data(), term4.str().data(),"same");
      tree2->SetMarkerColor(PROTONCOL);   tree2->Draw(termB.str().data(), termP.str().data(),"same");
      tree2->SetMarkerColor(PIPLUSCOL);     tree2->Draw(termB.str().data(), termPIPLUS.str().data(),"same");//pi+
      tree2->SetMarkerColor(PIMINUSCOL);     tree2->Draw(termB.str().data(), termPIMINUS.str().data(),"same");//pi-
      tree2->SetMarkerColor(PI0COL);     tree2->Draw(termB.str().data(), termPI0.str().data(),"same");
      #endif
}



void SetTLegend(TLegend *leg, TGraph *grXZ)
{
  /*
  TColor *col1 = gROOT->GetColor(POSITCOL);
   TColor *col2 = gROOT->GetColor(ELECCOL);
   TColor *col3 = gROOT->GetColor(GAMMACOL);
   TColor *col4 = gROOT->GetColor(MUCOL);
   TColor *col5 = gROOT->GetColor(MUCOL);
   col1->SetAlpha(0.4);
   col2->SetAlpha(0.4);
   col3->SetAlpha(0.4);
   col4->SetAlpha(0.4);
   col5->SetAlpha(0.4);
  */
  
   float makerSize = 0.4;
   TGraph *dum1 = new TGraph();
      dum1->SetMarkerColor(MUCOL);
      dum1->SetMarkerStyle(21);
      dum1->SetMarkerSize(makerSize*2);
      dum1->SetLineWidth(5);
   TGraph *dum2 = new TGraph();
      dum2->SetMarkerColor(POSITCOL);
      dum2->SetMarkerStyle(21);
      dum2->SetMarkerSize(makerSize*2);
      dum2->SetLineWidth(5);
   TGraph *dum3 = new TGraph();
      dum3->SetMarkerColor(ELECCOL);
      dum3->SetMarkerStyle(21);
      dum3->SetMarkerSize(makerSize*2);
      dum3->SetLineWidth(5);
   TGraph *dum4 = new TGraph();
      dum4->SetMarkerColor(GAMMACOL);
      dum4->SetMarkerStyle(21);
      dum4->SetMarkerSize(makerSize*2);
   TGraph *dum5 = new TGraph();
      dum5->SetMarkerColor(PROTONCOL);
      dum5->SetMarkerStyle(21);
      dum5->SetMarkerSize(makerSize*2);
   TGraph *dum6 = new TGraph();
      dum6->SetMarkerColor(NEUTRONCOL);
      dum6->SetMarkerStyle(21);
      dum6->SetMarkerSize(makerSize*2);
   TGraph *dum7 = new TGraph();
      dum7->SetMarkerColor(PIPLUSCOL);
      dum7->SetMarkerStyle(21);
      dum7->SetMarkerSize(makerSize*2);      
   TGraph *dum8 = new TGraph();
      dum8->SetMarkerColor(PIMINUSCOL);
      dum8->SetMarkerStyle(21);
      dum8->SetMarkerSize(makerSize*2);    
   TGraph *dum9 = new TGraph();
      dum9->SetMarkerColor(PI0COL);
      dum9->SetMarkerStyle(21);
      dum9->SetMarkerSize(makerSize*2);  
      
   leg->SetFillColor(10);
   leg->SetTextSize(0.03);
   leg->SetHeader("G4 Step points");
   leg->AddEntry(dum1,"#mu","p");
   leg->AddEntry(dum3,"e^{#minus}","p"); 
   leg->AddEntry(dum2,"e^{+}","p");
   leg->AddEntry(dum4,"#gamma","p");
   leg->AddEntry(dum5,"proton","p");
   leg->AddEntry(dum6,"neutron","p");
   leg->AddEntry(dum7,"#pi^{+}","p");
   leg->AddEntry(dum8,"#pi^{#minus}","p");
   leg->AddEntry(dum9,"#pi^{0}","p");
   leg->AddEntry(grXZ,"readout point","p");
}
