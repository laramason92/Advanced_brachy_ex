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
// Code developed by:
// S. Guatelli, D. Cutajar, J. Poder, 
//Centre For Medical Radiation Physics, University of Wollongong
//
//    *******************************
//    *                             *
//    *    BrachyFactoryGammaMed       *
//    *                             *
//    *******************************
//
//
#include "globals.hh"
#include "BrachyFactoryGammaMed.hh"
#include "G4ParticleTable.hh"
#include "Randomize.hh"  
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh" 
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstructionGammaMed.hh"

BrachyFactoryGammaMed:: BrachyFactoryGammaMed()
{
  gammaMedSource = new  BrachyDetectorConstructionGammaMed(); 
}

BrachyFactoryGammaMed:: ~BrachyFactoryGammaMed()
{
  delete gammaMedSource;
}
 
void BrachyFactoryGammaMed::CreateSource(G4VPhysicalVolume* mother)
{
  gammaMedSource -> ConstructGammaMed(mother);
}

void BrachyFactoryGammaMed::MoveSource(G4double z)
{
 //gammaMedSource -> MoveSourceX(x);
 //gammaMedSource -> MoveSourceY(y);
 G4cout << "The move command has been passed to the Factory" << G4endl;
 gammaMedSource -> MoveSourceZ(z);
 
}

void BrachyFactoryGammaMed::CleanSource()
{
  gammaMedSource -> CleanGammaMed();
  gammaMedSource = 0;
}
