/*


*/

/*
TCanvas
TGraph *grYZ = new TGraph();
TGraph *grXZ = new TGraph();
*/
void SetTLegend(TLegend *leg, TGraph* gra);
void DrawTrueStepping(TCanvas *c1, TTree *tree2, int i);

void event2(TString fname="mydata")
{
   gStyle->SetPadRightMargin(0.12);
	gStyle->SetPalette(55);
	gStyle->SetOptStat(0);
	gStyle->SetLineScalePS(1); // to make narrow line for marker
   gStyle->SetTitleSize(0.05,"xyz");

	TFile *tf = new TFile(fname+".root");
	TString canName = "fig_"+fname+".pdf";
	tf->Print();// << endl;
	cerr << "canName = " << canName << endl;

	//
	// dedx of fiber in layer
	//
	int nevt;
   const float fiber_thickness = 1;
   //const float fiber_thickness = 2;
   // const int unit = 52; // x + y + wavy50 = 52
   const int unit = 32; // x.y.sx30 = 32
   //const int unit = 17; // x.y.sx15 = 17
   const int actuallayers = unit * 10;
   //const int actuallayers = unit * 10;
   const int nlayers_dummy = 1200; 
   const int nfibers_dummy = 1200;
   float chw[nlayers_dummy][nfibers_dummy];//={};
	
   TTree *tree1 = static_cast<TTree *>(tf->Get("treeEvtAct2"));
   tree1->SetBranchAddress("nevt", &nevt);
   //tree1->SetBranchAddress("chx", chx);
   //tree1->SetBranchAddress("chy", chy);
   tree1->SetBranchAddress("chw", chw);

	//
	//	true trajectry stepping point
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
   //std::vector<TVector3> arrYZv3_sparce;
   TVector3 arrYZv3_sparce[999];
	/*
   because the order of saving Vec3 is not normal for wavy structure, so I do not use vector
   0  id = 0  TVector3 A 3D physics vector (x,y,z)=(-152.500000,0.500000,0.000000) 
   1  id = 2  TVector3 A 3D physics vector (x,y,z)=(-102.500000,0.500000,0.000000) 
   ...
   5  id = 10 TVector3 A 3D physics vector (x,y,z)=(  97.500000,0.500000,0.000000) 
   6  id = 1  TVector3 A 3D physics vector (x,y,z)=(-103.500000,0.500000,0.000000) 
   ...
   11 id = 11 TVector3 A 3D physics vector (x,y,z)=( 146.500000,0.500000,0.000000) 
	*/
   for (int i = 0; i < nReadOut_sparce; ++i) {
      tree4->GetEntry(i);
      cerr << i << " id = " << id << " "; v3->Print();
      //arrYZv3_sparce.push_back(*v3);
      arrYZv3_sparce[ id ] = *v3;
   }
	//return;


   int colN = TColor::GetNumberOfColors();
   cerr << "color max = " << colN << endl;
   // 0.3 is MIP position
   float maxdEdx1  = 0.50;
   float maxdEdx2  = 16.0;
   float colorStep1= maxdEdx1/colN;
   float colorStep2= maxdEdx2/colN;

	//
   // just dummy to get palette
	//
   TH2F *colpale = new TH2F("colpale","colpale",100,-1,1,100,-1,1);
	colpale->SetBinContent(0,0,maxdEdx1);
   colpale->GetZaxis()->SetRangeUser(0,2);
   colpale->GetZaxis()->SetTitle("Energy deposit [MeV]");
   colpale->GetZaxis()->SetTitleOffset(0.55);
   colpale->GetZaxis()->SetLabelSize(0.030);
   colpale->GetZaxis()->SetLabelOffset(0.005);
   TCanvas *c0 = new TCanvas();
   colpale->Draw("colz");
   gPad->Update();
   TPaletteAxis *palette = (TPaletteAxis*)colpale->GetListOfFunctions()->FindObject("palette");


	TH1F *dedx = new TH1F("dedx","dedx",1000,0,5);
   TGraph *grYZ = new TGraph();
   TGraph *grXZ = new TGraph();
   grYZ->SetMarkerColor(kGray+1);
   grXZ->SetMarkerColor(kGray+1);
   grXZ->GetXaxis()->SetLimits(-20, fiber_thickness * actuallayers +20); // ok @21/09/14
   grYZ->GetXaxis()->SetLimits(-20, fiber_thickness * actuallayers +20); // ok @21/09/14

   TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
   c1->Divide(2,2);
	c1->cd(1); grYZ->SetTitle(";z axis; y axis");
	c1->cd(2); grXZ->SetTitle(";z axis; x axis");
   c1->Print(canName+"[", "pdf");

   TLegend *leg = new TLegend(0.20,0.19,0.33,0.35);
   SetTLegend(leg, grXZ);

	#if 1
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

         #if 0 // later this will be removed 
         if ( l%2==0 ) { // mother even, not rotated
		      grXZ->SetPoint(grXZ->GetN(), arrXZv3[i].Z() + z_mother_position,
                                         arrXZv3[i].X() /* y-axis on a plot is x poisitoon of the fiber */ );
		   }
   		if ( l%2==1 ) { // mother odd
            grYZ->SetPoint(grYZ->GetN(), arrYZv3[i].Z() + z_mother_position, 
                                         arrYZv3[i].Y() /* y-axis on a plot is y poisitoon of the fiber */ );
		   }
		   #endif

         if ( l%unit==0 ) grXZ->SetPoint(grXZ->GetN(), arrXZv3[i].Z() + z_mother_position, 
                                                       arrXZv3[i].X() /* y-axis on a plot is x poisitoon of the fiber */ );
         if ( l%unit==1 ) grYZ->SetPoint(grYZ->GetN(), arrYZv3[i].Z() + z_mother_position, 
                                                       arrYZv3[i].Y() /* y-axis on a plot is y poisitoon of the fiber */ );

		}
	}
	#endif

   #if 1
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

      #if 0 // sparce layer 
      else  grYZ->SetPoint(grYZ->GetN(), arrYZv3_sparce[i].Z() + z_mother_position,
                                         arrYZv3_sparce[i].X() /* y-axis on a plot is y poisitoon of the fiber */ );
		#endif
		#if 1 // wavy layer
      else  grXZ->SetPoint(grXZ->GetN(), arrYZv3_sparce[i].Z() + z_mother_position,
                                         arrYZv3_sparce[i].X() /* y-axis on a plot is y poisitoon of the fiber */ );
		#endif
		}
	}
	#endif


   //
   // event loop
   // 
	for (int i=0; i<2; i++) {
	//for (int i=0; i<50; i++) {
		tree1 ->GetEntry(i);
		cerr << "nevt = " << nevt << endl;
   	c1->cd(1); grYZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(2); grYZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(3); grXZ->Draw("ap"); palette->Draw(); // ok @21/09/14
   	c1->cd(4); grXZ->Draw("ap"); palette->Draw(); // ok @21/09/14 
		//return;

		#if 1
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

            #if 0 // later this will be removed
				if ( l%2==0 ) { // mother even, not rotated
				   int direction=0;
               if ( wavyLayerXID%2==0 ) direction=+1; // id = 1, 3, 5...
               if ( wavyLayerXID%2==1 ) direction=-1; // id = 1, 3, 5...
	            // if there is amplitude due to wavy structure 
	            for (int wavy = 0; wavy < 1; wavy++) { 
	            //for (int wavy = 0; wavy < 12; wavy++) { 
         		     TMarker* mk = new TMarker( arrXZv3[ID].Z() + z_mother_position, 
                                               arrXZv3[ID].X() + direction * wavy * 2, 
                                               21);
		 				 //cerr << " dedx = " << chw[l][ID] << " step = " << colorStep1
		 			    //     << " col = " << (int) ( chw[l][ID] / colorStep1 ) << endl;
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )   ) );
   	          	 mk->SetMarkerSize(0.4);
                   c1->cd(3) /* XZ */;	mk->Draw();
                   c1->cd(4) /* XZ */; mk->Draw();
					}
      		}
            if ( l%2==1 ) { // mother odd
					int direction=0;
               if ( wavyLayerXID%2==0 ) direction=+1; // id = 1, 3, 5...
               if ( wavyLayerXID%2==1 ) direction=-1; // id = 1, 3, 5...
               //if ( wavyLayerXID%2==0 ) direction=-1; // id = 1, 3, 5...
               //if ( wavyLayerXID%2==1 ) direction=+1; // id = 1, 3, 5...
	            // if there is amplitude due to wavy structure 
	            for (int wavy = 0; wavy < 1; wavy++) { 
	            //for (int wavy = 0; wavy < 12; wavy++) { 
                   TMarker* mk = new TMarker( arrYZv3[ID].Z() + z_mother_position, 
                                              arrYZv3[ID].Y() + direction * wavy * 2, 
                                              21);
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )  ) );
                   mk->SetMarkerSize(0.4);
                   c1->cd(1)/* YZ */; mk->Draw();
                   c1->cd(2)/* YZ */; mk->Draw();
					}
				}
            #endif


            #if 1
            //
            // x/y fiber layer, and sparce fiber layer
            //
            if      ( l%unit==0 ) { //
                    TMarker* mk = new TMarker( arrXZv3[ID].Z() + z_mother_position, 
                                               arrXZv3[ID].X(), 21);
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )   ) );
                   mk->SetMarkerSize(0.4);
                   c1->cd(3) /* XZ */; mk->Draw();
                   c1->cd(4) /* XZ */; mk->Draw();
 				}
            else if ( l%unit==1 ) { //
                   TMarker* mk = new TMarker( arrYZv3[ID].Z() + z_mother_position, 
                                             -arrYZv3[ID].Y(), 21);
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )  ) );
                   mk->SetMarkerSize(0.4);
                   c1->cd(1)/* YZ */; mk->Draw();
                   c1->cd(2)/* YZ */; mk->Draw();
 				}
				#if 0 // sparce x layer 
            else { 
                  TMarker* mk = new TMarker( arrYZv3_sparce[ID].Z() + z_mother_position,
                                            -arrYZv3_sparce[ID].X(), 21);
                  mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )  ) );
                  mk->SetMarkerSize(0.4);
                  c1->cd(1)/* YZ */; mk->Draw();
            		c1->cd(2)/* YZ */; mk->Draw();
				}
				#endif
				#if 1 // wavy layer
            else { 
                  int direction = 0;
                  if ( ID % 2 == 0 ) direction = +1;
                  if ( ID % 2 == 1 ) direction = -1;
                  // 
                  for (int wavystep = 0; wavystep < 24; wavystep++) {  
                  TMarker* mk = new TMarker( arrYZv3_sparce[ID].Z() + z_mother_position,
                                             arrYZv3_sparce[ID].X() + wavystep*direction, 21);
                  //cerr << " => " << arrYZv3_sparce[ID].Z() + z_mother_position << ", " << arrYZv3_sparce[ID].X() << endl; 
                  mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )  ) );
                  mk->SetMarkerSize(0.4);
                  c1->cd(3)/* YZ */; mk->Draw();
                  c1->cd(4)/* YZ */; mk->Draw();
						}
            }
				#endif
            #endif
			}
		}
		#endif

      DrawTrueStepping(c1, tree2, i);

		#if 1 
      c1->cd(2);leg->Draw();
      c1->cd(4);leg->Draw();
      c1->Print(canName, "pdf");
		#endif
	}
   c1->Print(canName+"]", "pdf");
   TCanvas *c2 = new TCanvas("c2","c2",600,600);
	dedx->Draw();
}



void DrawTrueStepping(TCanvas *c1, TTree *tree2, int i)
{
      const int volLimX = 10000;//detwidth;
      const int volLimY = 10000;//detwidth;
      stringstream term1, term2, term3, term4;
      term1 << " detid!=0 && code==+13 &&evt==" << i << ends;
      term2 << " detid!=0 && code==+11 &&evt==" << i << ends;
      term3 << " detid!=0 && code==+22 &&evt==" << i << ends;
      term4 << " detid!=0 && code==-11 &&evt==" << i << ends;
      stringstream termP;
      stringstream termN;
      termP << " detid!=0 && code==2212 &&evt==" << i << ends;
      termN << " detid!=0 && code==2112 &&evt==" << i << ends;

      stringstream termA, termB;
      termA << "  y:z-1 " << ends;
      termB << "  x:z-1 " << ends;

      #if 1
      tree2->SetMarkerStyle(24);
      tree2->SetMarkerSize(0.6);
      tree2->SetLineWidth(1);

      c1->cd(2);
      tree2->SetMarkerColor(kGray+3);     tree2->Draw(termA.str().data(), term1.str().data(),"same");
      tree2->SetMarkerColor(kBlue+2);     tree2->Draw(termA.str().data(), term2.str().data(),"same");
      //tree2->SetMarkerColor(kGreen+2);    tree2->Draw(termA.str().data(), term3.str().data(),"same");
      tree2->SetMarkerColor(kRed+2);      tree2->Draw(termA.str().data(), term4.str().data(),"same");
      tree2->SetMarkerColor(kYellow+1);   tree2->Draw(termA.str().data(), termP.str().data(),"same");
      tree2->SetMarkerColor(kPink+1);     tree2->Draw(termA.str().data(), termN.str().data(),"same");
      c1->cd(4);
      tree2->SetMarkerColor(kGray+3);     tree2->Draw(termB.str().data(), term1.str().data(),"same");
      tree2->SetMarkerColor(kBlue+2);     tree2->Draw(termB.str().data(), term2.str().data(),"same");
      //tree2->SetMarkerColor(kGreen+2);    tree2->Draw(termB.str().data(), term3.str().data(),"same");
      tree2->SetMarkerColor(kRed+2);      tree2->Draw(termB.str().data(), term4.str().data(),"same");
      tree2->SetMarkerColor(kYellow+1);   tree2->Draw(termB.str().data(), termP.str().data(),"same");
      tree2->SetMarkerColor(kPink+1);     tree2->Draw(termB.str().data(), termN.str().data(),"same");
      #endif
}



void SetTLegend(TLegend *leg, TGraph *grXZ)
{
   TColor *col1 = gROOT->GetColor(kRed+2);
   TColor *col2 = gROOT->GetColor(kBlue+2);
   TColor *col3 = gROOT->GetColor(kGreen+2);
   TColor *col4 = gROOT->GetColor(kGray+1);
   TColor *col5 = gROOT->GetColor(kGray);
   col1->SetAlpha(0.4);
   col2->SetAlpha(0.4);
   col3->SetAlpha(0.4);
   col4->SetAlpha(0.4);
   col5->SetAlpha(0.4);

   float makerSize = 0.5;
   TGraph *dum1 = new TGraph();
      dum1->SetMarkerColor(kGray+1);
      dum1->SetMarkerStyle(24);
      dum1->SetMarkerSize(makerSize*2);
      dum1->SetLineWidth(5);
   TGraph *dum2 = new TGraph();
      dum2->SetMarkerColor(kRed-10);
      dum2->SetMarkerStyle(24);
      dum2->SetMarkerSize(makerSize*2);
      dum2->SetLineWidth(5);
   TGraph *dum3 = new TGraph();
      dum3->SetMarkerColor(kBlue-10);
      dum3->SetMarkerStyle(24);
      dum3->SetMarkerSize(makerSize*2);
      dum3->SetLineWidth(5);
   TGraph *dum4 = new TGraph();
      dum4->SetMarkerColor(kGreen-10);
      dum4->SetMarkerStyle(24);
      dum4->SetMarkerSize(makerSize*2);
   TGraph *dum5 = new TGraph();
      dum5->SetMarkerColor(kYellow+1);
      dum5->SetMarkerStyle(24);
      dum5->SetMarkerSize(makerSize*2);
   TGraph *dum6 = new TGraph();
      dum6->SetMarkerColor(kPink+1);
      dum6->SetMarkerStyle(24);
      dum6->SetMarkerSize(makerSize*2);
   leg->SetFillColor(10);
   leg->SetTextSize(0.03);
   leg->SetHeader("G4 Step points");
   leg->AddEntry(dum1,"#mu","p");
   leg->AddEntry(dum3,"e^{#minus}","p");
   leg->AddEntry(dum2,"e^{+}","p");
   leg->AddEntry(dum5,"proton","p");
   leg->AddEntry(dum6,"neutron","p");
   leg->AddEntry(grXZ,"readout point","p");
}
