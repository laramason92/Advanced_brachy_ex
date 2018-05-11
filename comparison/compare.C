{
// Read reference data in granero.txt
//FILE *fg1=fopen("granero.txt", "r");
//Int_t n_points_granero =13;
//Float_t x1[n_points_granero], y1[n_points_granero];
//Float_t x, y;
//Int_t ncols_granero;
//Int_t nlines1 =0;

//while(1)
// {
//  ncols_granero = fscanf(fg1,"%f %f",&x, &y);
//  if (ncols_granero<0) break;
//  //std::cout << "x " << x << std::endl;
//  x1[nlines1]=x;
//  y1[nlines1]=y;
//  //std::cout << "y " << y << std::endl;
//  nlines1++;
//}

FILE *fg1=fopen("Donaire.txt", "r");
Int_t n_points_donaire =200;
Float_t x1[n_points_donaire], y1[n_points_donaire], z1[n_points_donaire];
Float_t x, y, z;
Int_t ncols_donaire;
Int_t nlines1 =0;

while(1)
 {
  ncols_donaire = fscanf(fg1,"%f %f %f",&x, &y,&z);
  if (ncols_donaire<0) break;
  //std::cout << "x " << x << std::endl;
  x1[nlines1]=x;
  y1[nlines1]=y/1.111;
  z1[nlines1]=z;
  std::cout << "x" << x << "y " << y1[nlines1] << "z" << z <<std::endl;
  nlines1++;
}

fclose(fg1);

// Read the results of the brachytherapy advanced example
// FlexiSorceMacro.mac with 280 M events
FILE *fg2=fopen("geant4_dose.txt", "r");
Int_t n_points_geant4 =398;
Float_t x2[n_points_geant4], y2[n_points_geant4];
Int_t ncols_geant4;
Int_t nlines2 =0;

while(1)
 {
  ncols_geant4 = fscanf(fg2,"%f %f",&x, &y);
  if (ncols_geant4<0) break;
 // std::cout << "x " << x << std::endl;
  x2[nlines2]=x;
  y2[nlines2]=y;
  nlines2++;
}

fclose(fg2);

//TGraph *gr1 = new TGraph (nlines1, x1, y1);
TGraphErrors *gr1 = new TGraphErrors(nlines1,x1,y1,z1,z1);
TGraph *gr2 = new TGraph (nlines2, x2, y2);

TCanvas *c1 = new TCanvas("c1","Graph Draw Options",
                             200,10,600,400);

gPad->SetLogy();

// draw the graph with axis, continuous line, and put
// a * at each point
gr1->SetTitle("Radial dose deposition");
gr1-> GetXaxis()->SetTitle("r (cm)");
gr1->GetYaxis()->SetTitle("Normalised dose rate (z=0)");
gr1->SetLineWidth(1);
gr1->SetMarkerColor(1);
gr1->SetMarkerStyle(2);
gr1->SetMarkerSize(0.5);
gr1->Draw("AP");

gr2->SetLineWidth(1);
gr2->SetMarkerColor(2);
gr2->SetMarkerStyle(2);
gr2->SetMarkerSize(0.5);
gr2->SetLineColor(2);
gr2->Draw("CP");

TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
leg->SetFillColor(0);
leg->AddEntry(gr1, "Reference data", "lp");
leg->AddEntry(gr2, "Geant4 - 280e6 events", "lp");
leg->Draw();


}
