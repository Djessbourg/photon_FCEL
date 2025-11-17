void LO_frag_ratios() {
    // Ouvre les fichiers ROOT
    TFile *f1 = TFile::Open("ggoLOdiag_1.root");
    TFile *f3 = TFile::Open("ggoLOdiag_3.root");
    TFile *f5 = TFile::Open("ggoLOdiag_5.root");
    TFile *f7 = TFile::Open("ggoLOdiag_7.root");
    TFile *f8 = TFile::Open("ggoLOdiag_8.root");
    TFile *f9 = TFile::Open("ggoLOdiag_9.root");
    TFile *f13a = TFile::Open("ggoLOdiag_13a.root");
    TFile *f13b = TFile::Open("ggoLOdiag_13b.root");
    TFile *f15 = TFile::Open("ggoLOdiag_15.root");
    TFile *f14a = TFile::Open("ggoLOdiag_14a.root");
    TFile *f14b = TFile::Open("ggoLOdiag_14b.root");
    TFile *f16 = TFile::Open("ggoLOdiag_16.root");

    // Récupère les histogrammes (remplace "histo" par le nom de ton histogramme)
    TH1D *h1 = (TH1D*)f1->Get("hp20");
    TH1D *h3 = (TH1D*)f3->Get("hp20");
    TH1D *h5 = (TH1D*)f5->Get("hp20");
    TH1D *h7 = (TH1D*)f7->Get("hp20");
    TH1D *h8 = (TH1D*)f8->Get("hp20");
    TH1D *h9 = (TH1D*)f9->Get("hp20");
    TH1D *h13a = (TH1D*)f13a->Get("hp20");
    TH1D *h13b = (TH1D*)f13b->Get("hp20");
    TH1D *h15 = (TH1D*)f15->Get("hp20");
    TH1D *h14a = (TH1D*)f14a->Get("hp20");
    TH1D *h14b = (TH1D*)f14b->Get("hp20");
    TH1D *h16 = (TH1D*)f16->Get("hp20");

    // Vérifie qu'ils existent
    if (!h1 || !h3 || !h5 || !h7 || !h8 || !h9 || !h13a || !h13b || !h14a || !h14b || !h15 || !h16 ) {
        std::cerr << "Histogramme non trouvé dans l'un des fichiers !" << std::endl;
        return;
    }

    // Clone h1 pour créer un histogramme de somme
    TH1D *hsum = (TH1D*)h1->Clone("hsum");
    hsum->Add(h3);
    hsum->Add(h5);
    hsum->Add(h7);
    hsum->Add(h8);
    hsum->Add(h9);
    hsum->Add(h13a);
    hsum->Add(h13b);
    hsum->Add(h15);
    hsum->Add(h14a);
    hsum->Add(h14b);
    hsum->Add(h16);
    
    // Clone h1 pour créer un histogramme de rapport h1 / hsum 
    TH1D *h1ratio = (TH1D*)h1->Clone("1");
    h1ratio->Divide(hsum);

    // Clone h3 pour créer un histogramme de rapport h3 / hsum 
    TH1D *h3ratio = (TH1D*)h3->Clone("3");
    h3ratio->Divide(hsum);
    
    // Clone h5 pour créer un histogramme de rapport h5 / hsum 
    TH1D *h5ratio = (TH1D*)h5->Clone("5");
    h5ratio->Divide(hsum);
    
    // Clone h7 pour créer un histogramme de rapport h7 / hsum 
    TH1D *h7ratio = (TH1D*)h7->Clone("7");
    h7ratio->Divide(hsum);
    
    // Clone h8 pour créer un histogramme de rapport h8 / hsum 
    TH1D *h8ratio = (TH1D*)h8->Clone("8");
    h8ratio->Divide(hsum);
    
    // Clone h9 pour créer un histogramme de rapport h9 / hsum 
    TH1D *h9ratio = (TH1D*)h9->Clone("9");
    h9ratio->Divide(hsum);
    
    // Clone h13a pour créer un histogramme de rapport h13a / hsum 
    TH1D *h13aratio = (TH1D*)h13a->Clone("13a");
    h13aratio->Divide(hsum);
    
    // Clone h13b pour créer un histogramme de rapport h13b / hsum 
    TH1D *h13bratio = (TH1D*)h13b->Clone("13b");
    h13bratio->Divide(hsum);
    
    // Clone h14a pour créer un histogramme de rapport h14a / hsum 
    TH1D *h14aratio = (TH1D*)h14a->Clone("14a");
    h14aratio->Divide(hsum);
    
    // Clone h14b pour créer un histogramme de rapport h14b / hsum 
    TH1D *h14bratio = (TH1D*)h14b->Clone("14b");
    h14bratio->Divide(hsum);
    
    // Clone h15 pour créer un histogramme de rapport h15 / hsum 
    TH1D *h15ratio = (TH1D*)h15->Clone("15");
    h15ratio->Divide(hsum);
    
    // Clone h16 pour créer un histogramme de rapport h16 / hsum 
    TH1D *h16ratio = (TH1D*)h16->Clone("16");
    h16ratio->Divide(hsum);
    

    // Sauvegarde dans un nouveau fichier
    TFile *fout = new TFile("LO_frag_ratios.root", "RECREATE");
    h1ratio->Write();
    h3ratio->Write();
    h5ratio->Write();
    h7ratio->Write();
    h8ratio->Write();
    h9ratio->Write();
    h13aratio->Write();
    h13bratio->Write();
    h14aratio->Write();
    h14bratio->Write();
    h15ratio->Write();
    h16ratio->Write();
    fout->Close();
}
