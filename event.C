/*

*/
void event(TString fname="mydata")
//void event_dwavy(TString fname="gps_p0.5GeV_mu_water_BisY_toX_1")
{
	gStyle->SetPalette(55);
	gStyle->SetOptStat(0);
	gStyle->SetLineScalePS(1); // to make narrow line for marker
   gStyle->SetTitleSize(0.05,"xyz");

	TFile *tf = new TFile(fname+".root");
	TString canName = "fig_"+fname+".pdf";
	tf->Print();// << endl;
	cerr << "canName = " << canName << endl;
	//return;

	int nevt;
   const int nlayers_dummy = 599;
   const int nfibers_dummy = 699;
   float chw[nlayers_dummy][nfibers_dummy];//={};

   TTree *tr = static_cast<TTree *>(tf->Get("treeEvtAct2"));
   tr->SetBranchAddress("nevt", &nevt);
   //tr->SetBranchAddress("chx", chx);
   //tr->SetBranchAddress("chy", chy);
   tr->SetBranchAddress("chw", chw);

	//	true stepping point
   TTree *tin4 = (TTree*)tf->Get("treeStpAct");

	// readout position
   TTree *tin1 = (TTree*)tf->Get("tree_yread");
	int id;
   TVector3 *v3;// = new TVector3();
   tin1->SetBranchAddress("v3",&v3);
   tin1->SetBranchAddress("id",&id);

   int XX = tin1->GetEntries();
   //std::vector<TVector3> arrXZv3;
	TVector3 arrXZv3[999];
	TVector3 arrYZv3[999];
	cerr << "XX = " << XX << endl;

   for (int i = 0; i < XX; ++i) {
      tin1->GetEntry(i);
      //arrXZv3.push_back(*v3);
		if ( !v3 ) continue;
      cerr << i << " id = " << id << " "; v3->Print();
      arrXZv3[id] = *v3;
      //if (arrXZv3[i].Y()==0) cerr << "fiber : " <<  arrXZv3[i].Z() << endl;
      //grXZ->SetPoint(grXZ->GetN(), arrXZv3[id].Z(), arrXZv3[id].X());
      //arrBox1.push_back( new TBox(arrXZv3[i].Z() - 30, arrXZv3[i].Y() -1.5 , arrXZv3[i].Z() + 30, arrXZv3[i].Y() + 1.5) );
   }
	//return;

   int colN = TColor::GetNumberOfColors();
   cerr << "color max = " << colN << endl;
   // 0.3 is MIP position
   float maxdEdx1  = 1.50;
   float maxdEdx2  = 16.0;
   float colorStep1= maxdEdx1/colN;
   float colorStep2= maxdEdx2/colN;

   TColor *col1 = gROOT->GetColor(kRed+2);
   TColor *col2 = gROOT->GetColor(kBlue+2);
   TColor *col3 = gROOT->GetColor(kGreen+2);
   TColor *col4 = gROOT->GetColor(kGray+3);
   TColor *col5 = gROOT->GetColor(kGray);
   col1->SetAlpha(0.3);
   col2->SetAlpha(0.3);
   col3->SetAlpha(0.3);
   col4->SetAlpha(0.3);
   col5->SetAlpha(0.2);

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


   TGraph *grYZ = new TGraph();
   TGraph *grXZ = new TGraph();

	TH1F *dedx = new TH1F("dedx","dedx",1000,0,5);

   TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
   c1->Divide(2,2);
   c1->Print(canName+"[", "pdf");

	c1->cd(1);
	grYZ->Draw("ap"); // ok @21/09/14
	grYZ->GetXaxis()->SetLimits(-20, 1400); // ok @21/09/14
   grYZ->SetTitle(";z axis; y axis");
	c1->cd(2);
	grXZ->Draw("ap"); // ok @21/09/14
	grXZ->GetXaxis()->SetLimits(-20, 1400); // ok @21/09/14
   grXZ->SetTitle(";z axis; x axis");


   TLegend *leg = new TLegend(0.84,0.16,0.97,0.32);
   leg->SetFillColor(10);
   leg->SetTextSize(0.03);
   leg->SetHeader("Step point");
   leg->AddEntry(dum1,"#mu","p");
   leg->AddEntry(dum3,"e^{#minus}","p");
   leg->AddEntry(dum2,"e^{+}","p");
   leg->AddEntry(dum5,"proton","p");
   leg->AddEntry(grXZ,"readout point","p");


	#if 1 
	//
	//
	//
	#if 0 
	float globalshift = -40;
	#else
	//float globalshift = -520;
	//float globalshift = +520;
	float globalshift = 0;
	#endif

   for (int l = 0; l < 50; l++) {
   //for (int l = 0; l < 599; l++) {
   	//for (int i = 0; i < 999; ++i) {
   	for (int i = 0; i < 300/2; ++i) {

		if ( !arrXZv3 ) continue;

		//arrXZv3[i].Print();
   	const int layersX = 5; // ampli 50 mm, y 30 cm 
      //int z_mother_position = 48 * l; // 48mm
      int z_mother_position = 2 /* fiber thickness */ * l; // 48mm
      int wavyLayerZID, wavyLayerXID;

      if ( i%2==0) { // even
         wavyLayerZID = i / ((layersX+1)*2);
         wavyLayerXID = i % ((layersX+1)*2);
      }
      if ( i%2==1) { // odd
         wavyLayerZID = i / ((layersX)*2);
         wavyLayerXID = i % ((layersX)*2);
      }

      if ( l%2==0 ) { // mother even, not rotated
			c1->cd(1);
		   grXZ->SetPoint(grXZ->GetN(), 
            	              arrXZv3[i].Z() + z_mother_position,
                             arrXZv3[i].X() /* y-axis on a plot is x poisitoon of the fiber */ );
			TBox *b = new TBox( arrXZv3[i].Z() + z_mother_position - 1, 
      		   	           arrXZv3[i].X() - 1 , 
         		   	        arrXZv3[i].Z() + z_mother_position + 1,
            		   	     arrXZv3[i].X() + 1 );
			b->SetLineColor(kGray);
			//b->Draw("l");
		}
		if ( l%2==1 ) { // mother odd
			c1->cd(2);
			TVector3 tmp = arrXZv3[i];

			TVector3 zaxis(0,0,1);
         // rotation around z 
			// mother volume center is 0, and it is rotated around z
         tmp.Rotate(TMath::Pi()/2., zaxis); 
         //tmp.SetY(tmp.Y()*-1);
			//tmp.RotateZ(90);
			arrYZv3[i] = tmp;
			//cerr << "i " << i << " "; arrXZv3[i].Print();
			//cerr << "i " << i << " "; arrYZv3[i].Print();

         grYZ->SetPoint(grYZ->GetN(), 
                             arrYZv3[i].Z() + z_mother_position, 
                             arrYZv3[i].Y() + globalshift /* y-axis on a plot is y poisitoon of the fiber */ );
			TBox *b = new TBox( arrYZv3[i].Z() + z_mother_position - 1, 
                             arrYZv3[i].Y() + globalshift - 1 , 
                             arrYZv3[i].Z() + z_mother_position + 1, 
                             arrYZv3[i].Y() + globalshift + 1 );
	      b->SetLineColor(kGray);
			//b->Draw("l");
		}
		}
	}
	#endif


   //return;


	#if 1
	for (int i=0; i<1; i++) {
	//for (int i=0; i<50; i++) {
		tr ->GetEntry(i);
		//tr2->GetEntry(i);
		cerr << "nevt = " << nevt << endl;

   	c1->cd(1);
   	grYZ->Draw("ap"); // ok @21/09/14
   	c1->cd(2);
   	grYZ->Draw("ap"); // ok @21/09/14
   	c1->cd(3);
   	grXZ->Draw("ap"); // ok @21/09/14
   	c1->cd(4);
   	grXZ->Draw("ap"); // ok @21/09/14
   	grXZ->GetXaxis()->SetLimits(-20, 120); // ok @21/09/14
   	grYZ->GetXaxis()->SetLimits(-20, 120); // ok @21/09/14

	#if 1
   	//for (int l = 0; l < 30; l++) {
   	for (int l = 0; l < 599; l++) {
	      //for (int ID= 0; ID< 999; ID++) {
	      for (int ID= 0; ID< 599; ID++) {

   	      if (chw[l][ID]<1E-10) continue;
				dedx->Fill(chw[l][ID]);

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
            	  << " edep = " << chw[l][ID] 
              	; arrXZv3[ID].Print();
	              //<< endl;
	#endif
				//int z_mother_position = 48 * l;
            int z_mother_position = 2 /* fiber thickness */ * l; // 48mm

				if ( l%2==0 ) { // mother even, not rotated
		
					int direction=0;
               if ( wavyLayerXID%2==0 ) direction=+1; // id = 1, 3, 5...
               if ( wavyLayerXID%2==1 ) direction=-1; // id = 1, 3, 5...

	            for (int wavy = 0; wavy < 1; wavy++) { // if there is amplitude due to wavy structure 
         		     TMarker* mk = new TMarker( arrXZv3[ID].Z() + z_mother_position, 
                                               arrXZv3[ID].X() + direction*wavy*2, 
                                               21);
	         		 //mk->SetMarkerColor( kRed-10 );
						 #if 0
		 				 cerr << " dedx = " << chw[l][ID] << " step = " << colorStep1
		 			         << " col = " << (int) ( chw[l][ID] / colorStep1 )
                        << endl;
						 #endif
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )   ) );
   	          	 mk->SetMarkerSize(0.4);
                   c1->cd(3) /* XZ */;	mk->Draw();
                   c1->cd(4) /* XZ */; mk->Draw();
					}
      		}
            if ( l%2==1 ) { // mother odd

					int direction=0;
               //if ( wavyLayerXID%2==0 ) direction=+1; // id = 1, 3, 5...
               //if ( wavyLayerXID%2==1 ) direction=-1; // id = 1, 3, 5...
               if ( wavyLayerXID%2==0 ) direction=-1; // id = 1, 3, 5...
               if ( wavyLayerXID%2==1 ) direction=+1; // id = 1, 3, 5...

	            for (int wavy = 0; wavy < 1; wavy++) { // if there is amplitude due to wavy structure 
                   TMarker* mk = new TMarker( arrYZv3[ID].Z() + z_mother_position, 
                                              arrYZv3[ID].Y() + globalshift + direction*wavy*2, 
                                              21);
                   //mk->SetMarkerColor( kRed-10 );
                   mk->SetMarkerColor( TColor::GetColorPalette(   (int) ( chw[l][ID] / colorStep1 )  ) );
                   mk->SetMarkerSize(0.4);
                   c1->cd(1)/* YZ */; mk->Draw();
                   c1->cd(2)/* YZ */; mk->Draw();
					}
				}
			}
		}
		#endif

      const int volLimX = 10000;//detwidth;
      const int volLimY = 10000;//detwidth;
      stringstream term1, term2, term3, term4;
      stringstream term5, term6, term7, term8;
#if 0
      term1 << " detid==1 && code==+13 &&evt==" << i << ends;
      term2 << " detid==1 && code==+11 &&evt==" << i << ends;
      term3 << " detid==1 && code==+22 &&evt==" << i << ends;
      term4 << " detid==1 && code==-11 &&evt==" << i << ends;
      term5 << " detid==2 && code==+13 &&evt==" << i << ends;
      term6 << " detid==2 && code==+11 &&evt==" << i << ends;
      term7 << " detid==2 && code==+22 &&evt==" << i << ends;
      term8 << " detid==2 && code==-11 &&evt==" << i << ends;
#else
      term1 << " detid!=0 && code==+13 &&evt==" << i << ends;
      term2 << " detid!=0 && code==+11 &&evt==" << i << ends;
      term3 << " detid!=0 && code==+22 &&evt==" << i << ends;
      term4 << " detid!=0 && code==-11 &&evt==" << i << ends;
      term5 << " detid!=0 && code==+13 &&evt==" << i << ends;
      term6 << " detid!=0 && code==+11 &&evt==" << i << ends;
      term7 << " detid!=0 && code==+22 &&evt==" << i << ends;
      term8 << " detid!=0 && code==-11 &&evt==" << i << ends;

		stringstream termP;
      termP << " detid!=0 && code==2212 &&evt==" << i << ends;
#endif
		stringstream termA, termB;
		//termA << " -y:z-1 " << ends; 
		termA << " -y:z " << ends; 
      termB << "  x:z " << ends;

#if 1 
      tin4->SetMarkerStyle(24);
      tin4->SetMarkerSize(0.6);
      tin4->SetLineWidth(1);

      c1->cd(2);
      tin4->SetMarkerColor(kGray+3);     tin4->Draw(termA.str().data(), term1.str().data(),"same");
      tin4->SetMarkerColor(kBlue+2);     tin4->Draw(termA.str().data(), term2.str().data(),"same");
      //tin4->SetMarkerColor(kGreen+2);    tin4->Draw(termA.str().data(), term3.str().data(),"same");
      tin4->SetMarkerColor(kRed+2);      tin4->Draw(termA.str().data(), term4.str().data(),"same");
      tin4->SetMarkerColor(kYellow+1);   tin4->Draw(termA.str().data(), termP.str().data(),"same");
      c1->cd(4);
      tin4->SetMarkerColor(kGray+3);     tin4->Draw(termB.str().data(), term1.str().data(),"same");
      tin4->SetMarkerColor(kBlue+2);     tin4->Draw(termB.str().data(), term2.str().data(),"same");
      //tin4->SetMarkerColor(kGreen+2);    tin4->Draw(termB.str().data(), term3.str().data(),"same");
      tin4->SetMarkerColor(kRed+2);      tin4->Draw(termB.str().data(), term4.str().data(),"same");
      tin4->SetMarkerColor(kYellow+1);   tin4->Draw(termB.str().data(), termP.str().data(),"same");
#endif

		#if 1 
      c1->cd(2);leg->Draw();
      c1->cd(4);leg->Draw();
      c1->Print(canName, "pdf");
		#endif
	}

	//dedx->Draw();

   c1->Print(canName+"]", "pdf");
#endif
}
