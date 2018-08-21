//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
/*
Author: Susanna Guatelli
*/
// The class BrachyAnalysisManager creates and manages histograms and ntuples

// The analysis was included in this application following the extended Geant4
// example analysis/AnaEx01

#include <stdlib.h>
#include "G4SystemOfUnits.hh"
#include "BrachyAnalysisManager.hh"
#include "G4UnitsTable.hh"

BrachyAnalysisManager* BrachyAnalysisManager::instance = 0;

BrachyAnalysisManager::BrachyAnalysisManager()
{
#ifdef ANALYSIS_USE
 theTFile = 0;
 //histo6 =0;
 histo5 =0;
 histo4 =0;
 histo3 =0;
 //histo2_5cm =0;
 //histo2_4cm =0;
 //histo2_3cm =0;
 //histo2_2cm =0;
 //histo2_1cm =0;
 histo2 =0;
 histo = 0;
#endif
}

BrachyAnalysisManager::~BrachyAnalysisManager() 
{ 
#ifdef G4ANALYSIS_USE
 delete theTFile; theTFile = 0;
 //delete histo6; histo6 = 0;
 delete histo5; histo5 = 0;
 delete histo4; histo4 = 0;
 delete histo3; histo3 = 0;
 //delete histo2_5cm; histo2_5cm = 0;
 //delete histo2_4cm; histo2_4cm = 0;
 //delete histo2_3cm; histo2_3cm = 0;
 //delete histo2_2cm; histo2_2cm = 0;
 //delete histo2_1cm; histo2_1cm = 0;
 delete histo2; histo2 = 0;
 delete histo; histo = 0;
#endif
}

BrachyAnalysisManager* BrachyAnalysisManager::GetInstance()
{
  if (instance == 0) instance = new BrachyAnalysisManager;
  return instance;
}

void BrachyAnalysisManager::book() 
{  
#ifdef ANALYSIS_USE
 delete theTFile;
 theTFile = new TFile("brachytherapy.root", "RECREATE");
 
 histo = new TH1F("h10","energy spectrum", 1000, 0., 1000);
 histo2 = new TH2F("h20","edep2Dxy", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo2_1cm = new TH2F("h20_1cm","edep2Dxy_1cm", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo2_2cm = new TH2F("h20_2cm","edep2Dxy_2cm", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo2_3cm = new TH2F("h20_3cm","edep2Dxy_3cm", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo2_4cm = new TH2F("h20_4cm","edep2Dxy_4cm", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo2_5cm = new TH2F("h20_5cm","edep2Dxy_5cm", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm

 histo3 = new TH2F("h30","kerma2Dxy_voxelised_binning", 400, -50, 50, // binning, xmin, xmax, along x direction in cm WANT 1CM THICK BINS
				     400, 950, 1050);// binning, ymin, ymax, along y direction in cm

 histo5 = new TH2F("h50","kerma_mm_binning_2Dxy", 1, -50, 50,  // binning, xmin, xmax, along x direction in mm binning in the same way as the dose hist
				     1, 950, 1050);// binning, ymin, ymax, along y direction in mm

 //histo3 = new TH2F("h30","kerma2Dxy", 2001, -1000.5, 1000.5, // binning, xmin, xmax, along x direction in mm binning in the same way as the dose hist
 //				     2001, -1000.5, 1000.5);// binning, ymin, ymax, along y direction in mm
 histo4 = new TH2F("hgeom","edep2Dzy", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
 //histo6 = new TH3F("h60","edep3Dxyz", 801, -100.125, 100.125, // binning, xmin, xmax, along x direction in mm
 //				     801, -100.125, 100.125,// binning, ymin, ymax, along y direction in mm
 //				     801, -100.125, 100.125);// binning, ymin, ymax, along y direction in mm
#endif 
 }

#ifdef ANALYSIS_USE
void BrachyAnalysisManager::FillH2WithEnergyDeposition(G4double xx,
                                                     G4double yy, 
                                                     G4double energyDep)
{
  histo2 -> Fill(xx, yy,energyDep);
}
//void BrachyAnalysisManager::FillH2_1cm_WithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double energyDep)
//{
//  histo2_1cm -> Fill(xx, yy,energyDep);
//}
//void BrachyAnalysisManager::FillH2_2cm_WithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double energyDep)
//{
//  histo2_2cm -> Fill(xx, yy,energyDep);
//}
//void BrachyAnalysisManager::FillH2_3cm_WithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double energyDep)
//{
//  histo2_3cm -> Fill(xx, yy,energyDep);
//}
//void BrachyAnalysisManager::FillH2_4cm_WithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double energyDep)
//{
// histo2_4cm -> Fill(xx, yy,energyDep);
//}
//void BrachyAnalysisManager::FillH2_5cm_WithEnergyDeposition(G4double xx,
//                                                     G4double yy, 
//                                                     G4double energyDep)
//{
//  histo2_5cm -> Fill(xx, yy,energyDep);
//}

void BrachyAnalysisManager::FillH3WithKerma(G4double xx,
                                            G4double yy,
                                            G4double Kerma)
{
  histo3 -> Fill(xx, yy, Kerma);
}

void BrachyAnalysisManager::FillHGeomWithEnergyDeposition(G4double zz,
                                                     G4double yy, 
                                                     G4double energyDep)
{
  histo4 -> Fill(zz, yy,energyDep);
}

void BrachyAnalysisManager::FillH5WithKerma(G4double xx,
                                            G4double yy,
                                            G4double Kerma)
{
  histo5 -> Fill(xx, yy, Kerma);
}

//void BrachyAnalysisManager::FillH6WithEnergyDeposition(G4double xx,
//                                            G4double yy,
//                                            G4double zz,
//                                            G4double energyDep)
//{
//  histo6 -> Fill(xx, yy, zz, energyDep);
//}



void BrachyAnalysisManager::FillPrimaryParticleHistogram(G4double primaryParticleEnergy)
{
 // 1DHistogram: energy spectrum of primary particles  
  histo-> Fill(primaryParticleEnergy);
}
#endif 

void BrachyAnalysisManager::save() 
{  
#ifdef ANALYSIS_USE
 if (theTFile)
    {
	theTFile -> Write(); 
	theTFile -> Close();
    }
#endif
}












