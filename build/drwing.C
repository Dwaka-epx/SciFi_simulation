#define gamma_color kRed
#define e_color kBlue


void drwing(){
  TFile *file = new TFile("yprojection.root");

  TH1F* hist_gamma = (TH1F*)file->Get("Yprojection_mydata_gamma");
  TH1F* hist_e = (TH1F*)file->Get("Yprojection_mydata_e");
			
  hist_gamma->SetLineColor(gamma_color);
  hist_e->SetLineColor(kBlue);

  TCanvas *c1 = new TCanvas("c1","c1");
  hist_e->Draw();
    hist_gamma->Draw("same");
  
  c1->Update();
  TPaveStats *st1 = (TPaveStats*)hist_gamma->FindObject("stats");
  st1->SetTextColor(gamma_color);
  TPaveStats *st2 = (TPaveStats*)hist_e->FindObject("stats");
  st2->SetTextColor(e_color);
  // st2->SetX1NDC(0.18);
  // st2->SetX2NDC(0.38);
  // c1->Modified();

  int one_photo_bin = hist_gamma->GetXaxis()->FindBin(1.);
  int max_photo_bin = hist_gamma->GetXaxis()->FindBin(300.);

  for(int i=0;i++;i<300){
    double a= hist_gamma->Integral(one_photo_bin,i);
    double b= hist_e->Integral(i,max_photo_bin);
    if(abs(a-b)<10){
      cout << "|a-b|<10 when bin # =" << i <<endl;
    }
  }

  // cout << "one_photo_bin_g = "<<one_photo_bin_g <<endl;
  // cout << "one_photo_bin_e = "<<one_photo_bin_e <<endl;
  // cout << "max_photo_bin_g = "<<max_photo_bin_g <<endl;
  // cout << "max_photo_bin_e = "<<max_photo_bin_e <<endl;
}
