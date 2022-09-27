/*


*/

class MyMainFrame : public TGMainFrame 
{
	public:
   MyMainFrame(const TGWindow *p, TString detfile, TString datfile);
   //virtual ~MyMainFrame();
   //void DoSetlabel();
	void EventDraw();
   void DrawGLview();

	private:
   TGeoManager *gGeoManager;
   TGNumberEntry *guiEvent;
   TEvePointSet *pint_el;
   TEvePointSet *pint_po;
   TEvePointSet *pint_gm;
   TEvePointSet *pint_mu;
   TEvePointSet *pint_pr;
   TEvePointSet *pint_nt;
	TTree* tree;
   ClassDef(MyMainFrame, 0)
};


MyMainFrame::MyMainFrame(const TGWindow *p, TString detfile, TString datfile)
{
    //--- Definition of a simple geometry
    //   gSystem->Load("libGeom");
    //new TGeoManager("myGeo", "my geometry");
    gGeoManager = new TGeoManager("myGeo", "my geometry");
    gGeoManager->Import(detfile);
   // 0: world
	// 1: FiberDetMother
	// 2: FiberLayerMother

    #if 1
    //
    // draw the ROOT box via GL viewer
    //
    //gGeoManager->DefaultColors();
    //gGeoManager->CloseGeometry();
    gGeoManager->SetVisLevel(2); 
    gGeoManager->GetTopVolume()->Draw("ogl");
    #endif
    //return;

	#if 1 
    //gGeoManager->GetVolume("FiberDetMother")  ->SetTransparency(1);
    gGeoManager->GetVolume("FiberLayerMother")->SetTransparency(80);
    gGeoManager->GetVolume("FiberLayerSparceMother")->SetTransparency(80);
    gGeoManager->GetVolume("FiberLayerSparceMother")->SetLineColor(kRed+1);
	#endif

	#if 0 
    gGeoManager->GetVolume("FbrCoat")->Print();
    gGeoManager->GetVolume("FbrCoat")->SetTransparency(92);
    gGeoManager->GetVolume("FbrScin")->Print();
    gGeoManager->GetVolume("FbrScin")->SetLineColor(kRed-10);
    //gGeoManager->GetVolume("FbrScin")->SetTransparency(98);
	#endif
	#if 1 
    gGeoManager->GetVolume("Sct1")->Print();
    gGeoManager->GetVolume("Sct1")->SetLineColor(kGreen-2);
    gGeoManager->GetVolume("Sct2")->SetLineColor(kGreen-2);
    //gGeoManager->GetVolume("Sct1")->SetTransparency(2);
    //gGeoManager->GetVolume("Sct2")->SetTransparency(2);
    gGeoManager->GetVolume("Joint1")->Print();
    gGeoManager->GetVolume("Joint1")->SetLineColor(kMagenta-9);
    gGeoManager->GetVolume("Joint2")->SetLineColor(kMagenta-9);
    //gGeoManager->GetVolume("Joint1")->SetTransparency(2);
    //gGeoManager->GetVolume("Joint2")->SetTransparency(2);
	#endif
    //gGeoManager->DefaultColors();
    //gGeoManager->CloseGeometry();

	//
	// TEveManager viewer
	//
   TEveManager::Create();
	TGeoNode* node = gGeoManager->GetTopNode();
	TEveGeoTopNode* en = new TEveGeoTopNode(gGeoManager, node);
   en->SetVisLevel(1);
	gEve->AddGlobalElement(en);

	#if 1 
	//
	// draw track
	//
   TFile *tf = new TFile(datfile);
   tree = (TTree*)tf->Get("treeStpAct"); // trajectry data

   pint_el = new TEvePointSet("electron step");
   pint_po = new TEvePointSet("positron step");
   pint_gm = new TEvePointSet("gamma step");
   pint_mu = new TEvePointSet("muon step");
   pint_pr = new TEvePointSet("proton step");
   pint_nt = new TEvePointSet("neutron step");
	#endif


	//
	// Get broweser and make someting and close it
	//
	TEveBrowser* browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);
	{
   	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   	frmMain->SetWindowName("XX GUI");
   	frmMain->SetCleanup(kDeepCleanup);

   	//TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
      TGVerticalFrame *hf = new TGVerticalFrame(frmMain);
   	{
   		TGLabel *lbl1 = new TGLabel(hf, "Go to event: "); 
			hf->AddFrame(lbl1);

    		guiEvent = new TGNumberEntry(hf, 0, 9,999, TGNumberFormat::kNESInteger,
                                      TGNumberFormat::kNEANonNegative,
                                      TGNumberFormat::kNELLimitMinMax,
                                      0, 99999);
    		hf->AddFrame(guiEvent);
    		guiEvent->Connect("ValueSet(Long_t)", "MyMainFrame", this, "EventDraw()");

         TGLabel *lbl2 = new TGLabel(hf, "Draw Geom on GLview"); 
			hf->AddFrame(lbl2);

         TGTextButton *quit = new TGTextButton(hf, "DrawGeom");
         hf->AddFrame(quit, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
                                              0, 0, 0, 5));
         //quit->Connect("Pressed()", "TApplication", gApplication, "Terminate()");
         quit->Connect("Pressed()", "MyMainFrame", this, "DrawGLview()");
   	}
   	frmMain->AddFrame(hf);

   	frmMain->MapSubwindows();
   	frmMain->Resize();
   	frmMain->MapWindow();
	}
   browser->StopEmbedding();
   browser->SetTabTitle("EvtDraw", 0);


	#if 0 
	TString outname;
	outname = "genfitGDMLGeom.root";
    TFile *outfile = TFile::Open(outname,"RECREATE");
    gGeoManager->Write();
    outfile->Close();
	#endif
}

void MyMainFrame::DrawGLview()
{
    gGeoManager->SetVisLevel(2);
    gGeoManager->GetTopVolume()->Draw("ogl");

}

//void takeTraject(TTree *tree, int evt, int pdg, int col, TEveStraightLineSet* lineSet);
//void takePoints (TTree *tree, int evt, int pdg, int col, TEvePointSet*        pintSet);

void takeTraject(TTree *tree, int evt, int pdg, int col, TEveStraightLineSet* lineSet)
{
	stringstream cut;
	//cut << "code==" << pdg << "&&evt==" << evt << ends;
	cut << "edep>0.01 && code==" << pdg << "&&evt==" << evt << ends;
    //TCanvas *c = new TCanvas();
    //tree->Draw("x:y:z", "detid!=0 && code==+11 &&evt==0", "");
    //tree->Draw("x:y:z", "detid!=0 && code==+11 &&evt==0", "goff");
    tree->Draw("x:y:z", cut.str().data(), "goff");
    int     size = tree->GetSelectedRows();
    double *val0 = tree->GetV1(); //
    double *val1 = tree->GetV2(); // 
    double *val2 = tree->GetV3(); // 
    cerr << "size = " << size << endl;
 
    //TEveStraightLineSet* lineSet = new TEveStraightLineSet;
    // SetLine (int idx, Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2)
    //lineSet->AddLine (-10,-10,-10, 1,1,1);
    //lineSet->AddMarker(-10,-10,-10);
    for (int i=0; i<size; i++) {
        lineSet->AddMarker(val0[i]/10.,val1[i]/10.,val2[i]/10.);
    }
    lineSet->SetLineColor(4);
    lineSet->SetLineWidth(3);
    lineSet->SetMarkerColor(col);
    lineSet->SetMarkerSize(1.4);
}

void takePoints (TTree *tree, int evt, int pdg, int col, TEvePointSet*        pintSet)
{
    stringstream cut;
    //cut << "code==" << pdg << "&&evt==" << evt << ends;
    cut << "edep>0.02 && code==" << pdg << "&&evt==" << evt << ends;
    //TCanvas *c = new TCanvas();
    //tree->Draw("x:y:z", "detid!=0 && code==+11 &&evt==0", "");
    //tree->Draw("x:y:z", "detid!=0 && code==+11 &&evt==0", "goff");
    tree->Draw("x:y:z", cut.str().data(), "goff");
    int     size = tree->GetSelectedRows();
    double *val0 = tree->GetV1(); //
    double *val1 = tree->GetV2(); // 
    double *val2 = tree->GetV3(); // 
    cerr << "size = " << size << endl;

    pintSet->SetOwnIds(kTRUE);
    for (int i=0; i<size; i++) {
      pintSet->SetNextPoint(val0[i]/10.,val1[i]/10.,val2[i]/10.);
      pintSet->SetPointId(new TNamed(Form("Point %d", i), ""));
    }
    pintSet->SetMarkerColor(col);
    pintSet->SetMarkerSize(2);
    pintSet->SetMarkerStyle(4);
}

void MyMainFrame::EventDraw()
{
   int num = guiEvent->GetNumberEntry()->GetIntNumber();
	cerr << "EventDraw = " << num << endl;


   #if 0
    TEveStraightLineSet* line_el = new TEveStraightLineSet("electron");
    TEveStraightLineSet* line_po = new TEveStraightLineSet("positron");
    TEveStraightLineSet* line_gm = new TEveStraightLineSet("gamma");
    TEveStraightLineSet* line_mu = new TEveStraightLineSet("muon");
   takeTraject(tree, 0/*evt*/, +13, kGray+1, line_mu);
   takeTraject(tree, 0/*evt*/, +11,       4, line_el);
   takeTraject(tree, 0/*evt*/, -11,       2, line_po);
   takeTraject(tree, 0/*evt*/, +22, kGreen+2,line_gm);
   gEve->AddElement(line_el);
   gEve->AddElement(line_po);
   gEve->AddElement(line_gm);
   gEve->AddElement(line_mu);
   #endif

   #if 1
   pint_el->Reset();
   pint_po->Reset();
   pint_gm->Reset();
   pint_mu->Reset();
   pint_pr->Reset();
   pint_nt->Reset();
   takePoints(tree, num/*evt*/, +13,  kGray+1, pint_mu);
   takePoints(tree, num/*evt*/, +11,        4, pint_el);
   takePoints(tree, num/*evt*/, -11,       +2, pint_po);
   takePoints(tree, num/*evt*/, +22, kGreen+2, pint_gm);
   takePoints(tree, num/*evt*/, +2212, kYellow+2,pint_pr);
   takePoints(tree, num/*evt*/, +2112, kPink+2,  pint_nt);
	gEve->AddElement(pint_el);
   gEve->AddElement(pint_po);
   gEve->AddElement(pint_gm);
   gEve->AddElement(pint_mu);
   gEve->AddElement(pint_pr);
   gEve->AddElement(pint_nt);
   #endif


   #if 0 // just test
    float a=10; float d=5; float x=0; float y=0; float z=0;
    TEveBox* b = new TEveBox("Box", "Test Title");
    b->SetMainColor(kCyan);
    b->SetMainTransparency(60);
    b->SetVertex(0, x - a, y - a, z - a);
    b->SetVertex(1, x - a, y + a, z - a);
    b->SetVertex(2, x + a, y + a, z - a);
    b->SetVertex(3, x + a, y - a, z - a);
    b->SetVertex(4, x - a, y - a, z + a);
    b->SetVertex(5, x - a, y + a, z + a);
    b->SetVertex(6, x + a, y + a, z + a);
    b->SetVertex(7, x + a, y - a, z + a);
    gEve->AddElement(b);
   #endif

    double refPoint[3] = {0.,0.,0.};
    // Int_t axesType = 0(Off), 1(EDGE), 2(ORIGIN), Bool_t axesDepthTest, Bool_t referenceOn, const Double_t referencePos[3]
    gEve->GetDefaultViewer()->GetGLViewer()->SetGuideState(2, kFALSE, kFALSE, refPoint);
    gEve->Redraw3D(kTRUE);
    //gEve->FullRedraw3D(kTRUE); 
    //gEve->FullRedraw3D(kFALSE, kTRUE);  
    //gEve->Redraw3D(kFALSE, kTRUE);
    //gEve->Redraw3D(kTRUE, kFALSE);

}


void mydetGeom2()
{
	// starting up an event display, here
   new MyMainFrame(gClient->GetRoot(), "mydet.gdml", "mydata.root");
}

