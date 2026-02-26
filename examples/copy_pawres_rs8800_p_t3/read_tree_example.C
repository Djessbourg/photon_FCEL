{
#include "Riostream.h"
  TFile *f = new TFile("ggdtest_d.root");
	TTree *t2 = (TTree*)f->Get("t2");
	//
  t2->Print();
  //
	Int_t iprov,ntrack;
	Double_t e[3],px[3],py[3],pz[3];
  Double_t x3;
	Double_t pt[3],y[3];
  Double_t x1,x2;
	Float_t pdf_weight[1000];
	Float_t weight;
  // we get the value stored into the header for the normalisation
  TList* list = t2->GetUserInfo();
  list->Print();
  TVectorT<float> &v = *(list->At(0)); 
  float& nb_evt = v[0];
  float& xsec = v[1];
  float& sqrt_s = v[2];
  float norma = xsec/nb_evt;
  //
	t2->SetBranchAddress("iprov",&iprov);
	t2->SetBranchAddress("ntrack",&ntrack);
	t2->SetBranchAddress("x3",&x3);
	t2->SetBranchAddress("energy",e);
	t2->SetBranchAddress("px",px);
	t2->SetBranchAddress("py",py);
	t2->SetBranchAddress("pz",pz);
	t2->SetBranchAddress("pdf_weight",pdf_weight);
  //
  Int_t bin_pt = 40;
  Int_t bin_y = 40;
  Double_t pt_min = 10.;
  Double_t pt_max = 100.;
  Double_t y_min = -10.;
  Double_t y_max = 10.;
  Double_t bin_size_pt,bin_size_y;
  bin_size_pt = (pt_max-pt_min)/(Double_t) bin_pt;
  bin_size_y = (y_max-y_min)/(Double_t) bin_y;
  //
	TH1D *hpt = new TH1D("pt","essai",bin_pt,pt_min,pt_max);
	TH1D *hy = new TH1D("y","essai",bin_y,y_min,y_max);
	TH1F *hx1 = new TH1F("x1","essai",100,0.,1.);
	TH1F *hx2 = new TH1F("x2","essai",100,0.,1.);
  int nbin = 21; 
  float xbin[nbin];
  float xmin = 1.e-8;
  float xmax = 1.;
  float temp;
  int j;
  for (j=0;j<nbin;j++) {
    temp = log(xmin) + j*(log(xmax)-log(xmin))/float(nbin);
    xbin[j] = exp(temp);
  }
  xbin[nbin] = xmax;
	TH1F *hx1 = new TH1F("x1","essai",nbin-1,xbin);
	TH1F *hx2 = new TH1F("x2","essai",nbin-1,xbin);
  //
	Int_t entries = (Int_t)t2->GetEntries();
  //
  for (Int_t i=0;i<entries;i++) {
    t2->GetEntry(i);
    // Pt and y are built, 0 is the photon, 1 is the hard parton and 2 is the soft one
    for (Int_t j=0;j<ntrack;j++) {
      pt[j] = sqrt(px[j]*px[j]+py[j]*py[j]);
      y[j] = log( (e[j]+pz[j])/(e[j]-pz[j]) )*0.5;
    }
    weight = pdf_weight[0];
    // we build now x1 and x2
    // for the 2 --> 3 case
    if ( (iprov == 33) || (iprov == 44) ) {
      x1 = (pt[0]/x3*exp(y[0])+pt[1]*exp(y[1])+pt[2]*exp(y[2]))/sqrt_s;
      x2 = (pt[0]/x3*exp(-y[0])+pt[1]*exp(-y[1])+pt[2]*exp(-y[2]))/sqrt_s;
    }
    // for the other cases : 1 is the recoiling particle against the photon
    else {
      x1 = pt[1]*(exp(y[0])+exp(y[1]))/sqrt_s;
      x2 = pt[1]*(exp(-y[0])+exp(-y[1]))/sqrt_s;
    }
    //
    hpt->Fill(pt[0],weight);
    hy->Fill(y[0],weight);
    hx1->Fill(x1,weight);
    hx2->Fill(x2,weight);
  }
  // normalisation of the histograms
  hpt->Scale(norma/bin_size_pt);
  hy->Scale(norma/bin_size_y);
  hy->Print();
	//
	TCanvas *c1 = new TCanvas("c1","Graph Draw Options",205,47,600,400);
	c1->SetFillColor(0);

  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();
	// pour eviter le titre
  hy->SetTitle("");
	// les titres des axes
  hy->GetXaxis()->SetTitle("y");
  hy->GetYaxis()->SetTitle("d #sigma/d y");
	//
  hy->Draw("AXIS");
  hy->Draw("");
  c1->cd(2);
  gPad->SetLogy();
	// pour eviter le titre
  hpt->SetTitle("");
  hpt->SetMaximum(200000);
  hpt->SetMinimum(5.);
	// les titres des axes
  hpt->GetXaxis()->SetTitle("P_{t #gamma}");
  hpt->GetYaxis()->SetTitle("d #sigma/d P_{t #gamma}");
	//
  hpt->Draw("AXIS");
  hpt->Draw("");
  c1->cd(3);
  gPad->SetLogy();
	// pour eviter le titre
  hx1->SetTitle("");
	// les titres des axes
  hx1->GetXaxis()->SetTitle("x_{2}");
  hx1->GetYaxis()->SetTitle("d #sigma/d x_{2}");
	//
  hx1->Draw("AXIS");
  hx1->Draw("");
  c1->cd(4);
  gPad->SetLogy();
	// pour eviter le titre
  hx2->SetTitle("");
	// les titres des axes
  hx2->GetXaxis()->SetTitle("x_{2}");
  hx2->GetYaxis()->SetTitle("d #sigma/d x_{2}");
	//
  hx2->Draw("AXIS");
  hx2->Draw("");
}
