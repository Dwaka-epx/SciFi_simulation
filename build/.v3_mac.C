

void v3_mac(TString fname="mydata")
{
  TFile *tf = new TFile(fname+".root");
  TTree *tree3 = (TTree*)tf->Get("tree_yread");

     int id;
     TVector3 *v3=NULL;// = new TVector3();
     tree3->SetBranchAddress("v3",&v3);
          
     int nReadOut = tree3->GetEntries();
    

     #if 1
     for (int i = 0; i < nReadOut; ++i) {
       tree3->GetEntry(i);
       //if ( !v3 ) continue;
       v3->Print();
       /*
       //arrXZv3[id] = *v3;
       arrXZv3.push_back(*v3);
       // rotation around z 
       // mother volume center is 0, and it is rotated around z
       v3->Rotate(TMath::Pi()/2., zaxis);
       arrYZv3.push_back(*v3);
       */
     }
     #endif

}
