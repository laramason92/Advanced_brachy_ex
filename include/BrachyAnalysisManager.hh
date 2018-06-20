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
//
/*
Author: Susanna Guatelli

The class BrachyAnalysisManager creates and manages histograms and ntuples
*/

#ifndef BrachyAnalysisManager_HH
#define BrachyAnalysisManager_HH

#include "globals.hh"

#ifdef ANALYSIS_USE
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TH3F.h"
#endif


class BrachyAnalysisManager
{
private:
  BrachyAnalysisManager();

public:
  ~BrachyAnalysisManager();
  static BrachyAnalysisManager* GetInstance();

  void book();
  // Create the output ROOT file 
  // Create the ntuple and histograms

#ifdef ANALYSIS_USE

  void FillH2WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);

  //void FillH2_1cm_WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);
  //void FillH2_2cm_WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);
  //void FillH2_3cm_WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);
  //void FillH2_4cm_WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);
  //void FillH2_5cm_WithEnergyDeposition(G4double xx,G4double yy, G4double energyDep);
  // Method to fill the 2D histogram with the energy deposition, integrated over a run, in each voxel
  // of the scoring mesh. The scoring mesh is in the plane containing the source.
  
  void FillH3WithKerma(G4double xx,G4double yy, G4double Kerma);

  void FillHGeomWithEnergyDeposition(G4double zz,G4double yy, G4double energyDep);
  // Method to fill the 2D histogram with the energy deposition, integrated over a run, in each voxel
  // of the scoring mesh. The scoring mesh is in the plane containing the source.
  
  void FillH5WithKerma(G4double xx,G4double yy, G4double Kerma);

  //void FillH6WithEnergyDeposition(G4double xx,G4double yy, G4double zz, G4double energyDep);

  void FillPrimaryParticleHistogram(G4double);
  // Energy spectrum of primary particles
#endif

  void save();
 // This method if called at the end of the run to store the 
 // results in the ROOT file

private:
    static BrachyAnalysisManager* instance;
 
#ifdef ANALYSIS_USE
    TFile* theTFile;
    TH1F* histo;
    TH2F* histo2;
    //TH2F* histo2_1cm;
    //TH2F* histo2_2cm;
    //TH2F* histo2_3cm;
    //TH2F* histo2_4cm;
    //TH2F* histo2_5cm;
    TH2F* histo3;
    TH2F* histo4;
    TH2F* histo5;
    //TH3F* histo6;

#endif
};
#endif



