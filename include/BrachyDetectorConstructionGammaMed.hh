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
// Code developed by: S. Guatelli, D. Cutajar, J. Poder, 
//Centre For Medical Radiation Physics, University of Wollongong
//
//
//    ******************************************
//    *                                        *
//    *    BrachyDetectorConstructionGammaMed.hh  *
//    *                                        *
//    ******************************************
//
// 
//

#ifndef BrachyDetectorConstructionGammaMed_H
#define BrachyDetectorConstructionGammaMed_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class G4LogicalVolume;
class G4Tubs;
class G4Cons;
class G4Torus;
class G4Box;
class G4Sphere;
class G4VPhysicalVolume;
class BrachyMaterial;
class G4VisAttributes;
class G4GenericMessenger;

class BrachyDetectorConstructionGammaMed
{
public:
   BrachyDetectorConstructionGammaMed();
  ~BrachyDetectorConstructionGammaMed();

  void  ConstructGammaMed(G4VPhysicalVolume*);
  // Model the GammaMed iridium source

  void  MoveSourceX(G4double x);//, G4double y, G4double z);
  void  MoveSourceY(G4double y);//, G4double y, G4double z);
  void  MoveSourceZ(G4double z);//, G4double y, G4double z);

  void  CleanGammaMed(); 
  // Destroy the Iridium source in the experimental set-up

private:   

  //void  DefineCommands();

  G4GenericMessenger* fMessenger;

  G4double fSourceTransX;
  G4double fSourceTransY;
  G4double fSourceTransZ;

  //G4Box* crit_vol;
  //G4LogicalVolume* logical_crit_vol;
  //G4VPhysicalVolume* physical_crit_vol;

  G4Tubs* steel_shell;
  G4LogicalVolume* logical_steel_shell;    
  G4VPhysicalVolume* physical_steel_shell;

  G4Tubs* air_gap;
  G4LogicalVolume* logical_air_gap;
  G4VPhysicalVolume* physical_air_gap;

  G4Tubs* End1_steel_shell;
  G4LogicalVolume* logical_End1_steel_shell;
  G4VPhysicalVolume* physical_End1_steel_shell;

  G4Cons* End1cone_steel_shell;
  G4LogicalVolume* logical_End1cone_steel_shell;
  G4VPhysicalVolume* physical_End1cone_steel_shell;

  G4Cons* End2_steel_shell;
  G4LogicalVolume* logical_End2_steel_shell;
  G4VPhysicalVolume* physical_End2_steel_shell;

  G4Tubs* cable;
  G4LogicalVolume* logical_cable;
  G4VPhysicalVolume* physical_cable;
    
  G4Tubs* iridium_core;
  G4LogicalVolume* logical_iridium_core;
  G4VPhysicalVolume* physical_iridium_core;

  G4Torus* metal_ring;
  G4LogicalVolume* logical_metal_ring;
  G4VPhysicalVolume* physical_metal_ring;

  G4Torus* plas_ring;
  G4LogicalVolume* logical_plas_ring;
  G4VPhysicalVolume* physical_plas_ring;

  G4Torus* air_ring;
  G4LogicalVolume* logical_air_ring;
  G4VPhysicalVolume* physical_air_ring;

  G4Tubs* metal_rod1;
  G4LogicalVolume* logical_metal_rod1;
  G4VPhysicalVolume* physical_metal_rod1;

  G4Sphere* metal_rod1_end;
  G4LogicalVolume* logical_metal_rod1_end;
  G4VPhysicalVolume* physical_metal_rod1_end;

  G4Sphere* metal_rod2_end;
  G4LogicalVolume* logical_metal_rod2_end;
  G4VPhysicalVolume* physical_metal_rod2_end;

  G4Sphere* metal_rod2bent_end;
  G4LogicalVolume* logical_metal_rod2bent_end;
  G4VPhysicalVolume* physical_metal_rod2bent_end;

  G4Tubs* air_rod1;
  G4LogicalVolume* logical_air_rod1;
  G4VPhysicalVolume* physical_air_rod1;

  G4Tubs* metal_rod2;
  G4LogicalVolume* logical_metal_rod2;
  G4VPhysicalVolume* physical_metal_rod2;

  G4Tubs* air_rod2;
  G4LogicalVolume* logical_air_rod2;
  G4VPhysicalVolume* physical_air_rod2;

  G4Tubs* metal_rod2bent;
  G4LogicalVolume* logical_metal_rod2bent;
  G4VPhysicalVolume* physical_metal_rod2bent;

  G4Tubs* air_rod3;
  G4LogicalVolume* logical_air_rod3;
  G4VPhysicalVolume* physical_air_rod3;

  G4VisAttributes* steelAttributes;
  G4VisAttributes* endAttributes;
  G4VisAttributes* simpleIridiumVisAtt;
  G4VisAttributes* titaniumAttributes;
  G4VisAttributes* acetalAttributes;
  G4VisAttributes* airAttributes;
 
  BrachyMaterial* pMat;    
};
#endif








