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
//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Guatelli, D. Cutajar, J. Poder
// Centre For Medical Radiation Physics, University of Wollongong
// Edits by Lara Mason
// University of the Witwatersrand
// 
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionGammaMed.cc*
//    *                                      *
//    ****************************************
//
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionGammaMed.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "BrachyMaterial.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

BrachyDetectorConstructionGammaMed::BrachyDetectorConstructionGammaMed()
  : steel_shell(0),logical_steel_shell(0),air_gap(0), logical_air_gap(0), physical_air_gap(0),
    End1_steel_shell(0),logical_End1_steel_shell(0), physical_End1_steel_shell(0),
    End1cone_steel_shell(0),logical_End1cone_steel_shell(0), physical_End1cone_steel_shell(0),
    End2_steel_shell(0),logical_End2_steel_shell(0), physical_End2_steel_shell(0),
    cable(0),logical_cable(0),physical_cable(0),
    iridium_core(0),logical_iridium_core(0),physical_iridium_core(0),
    metal_ring(0),logical_metal_ring(0),physical_metal_ring(0),
    plas_ring(0),logical_plas_ring(0),physical_plas_ring(0),
    air_ring(0),logical_air_ring(0),physical_air_ring(0),
    metal_rod1(0),logical_metal_rod1(0),physical_metal_rod1(0),
    metal_rod1_end(0),logical_metal_rod1_end(0),physical_metal_rod1_end(0),
    metal_rod2_end(0),logical_metal_rod2_end(0),physical_metal_rod2_end(0),
    metal_rod2bent_end(0),logical_metal_rod2bent_end(0),physical_metal_rod2bent_end(0),
    air_rod1(0),logical_air_rod1(0),physical_air_rod1(0),
    metal_rod2(0),logical_metal_rod2(0),physical_metal_rod2(0),
    air_rod2(0),logical_air_rod2(0),physical_air_rod2(0),
    metal_rod2bent(0),logical_metal_rod2bent(0),physical_metal_rod2bent(0),
    air_rod3(0),logical_air_rod3(0),physical_air_rod3(0),
    steelAttributes(0), endAttributes(0), simpleIridiumVisAtt(0), titaniumAttributes(0)
{
  pMat = new BrachyMaterial();
}

BrachyDetectorConstructionGammaMed::~BrachyDetectorConstructionGammaMed()
{ 
  delete pMat; 
}

void BrachyDetectorConstructionGammaMed::ConstructGammaMed(G4VPhysicalVolume* mother)
{
  G4Material* steelMat = pMat -> GetMat("Stainless steel 304");
  G4Material* iridiumMat = pMat -> GetMat("Iridium");
  G4Material* airMat = pMat -> GetMat("Air");
  G4Material* titaniumMat = pMat -> GetMat("titanium");
  G4Material* acetalMat = pMat -> GetMat("acetal");

 //Define dimensions of the outer Steel shell around the solid source - not including the ends 

  G4double shellr_min = 0.00 * mm;
  G4double shellr_max = 0.9 * mm;
  G4double shell_length = 3.6 * mm; //4.96-0.15 
  G4double shelloffset_x = 0.0 * mm;
  G4double shelloffset_y = 0.0 * mm;
  G4double shelloffset_z = 0.05 * mm; 
  steel_shell = new G4Tubs("steel_shell",shellr_min, shellr_max/2, shell_length/2.,0.*deg,360.*deg);
  logical_steel_shell = new G4LogicalVolume(steel_shell, steelMat, "steel_shell_log", 0, 0, 0);
  physical_steel_shell = new G4PVPlacement(0,G4ThreeVector(shelloffset_x,shelloffset_y,shelloffset_z),"phys_steel_shell", logical_steel_shell, mother, false, 0, true);

//Define dimensions of the air gap between Steel shell and Iridium core
  G4double airr_min = 0.00 * mm;
  G4double airr_max = 0.7 * mm;
  G4double air_length = 3.6 * mm;
  G4double airoffset_x = 0.0 * mm;
  G4double airoffset_y = 0.0 * mm;
  G4double airoffset_z = 0.05 * mm; 
  air_gap = new G4Tubs("air_gap", airr_min, airr_max/2, air_length/2, 0.*deg, 360.*deg);
  logical_air_gap = new G4LogicalVolume(air_gap, airMat, "air_gap_log", 0, 0, 0);
  physical_air_gap = new G4PVPlacement(0, G4ThreeVector(airoffset_x,airoffset_y,airoffset_z), "phys_air_gap",logical_air_gap, physical_steel_shell, false, 0, true);

//Define the non-cable weld end of the Steel shell 
  G4double End1r_min = 0.0 * mm;
  G4double End1r_max = 0.9 * mm;
  G4double End1length = 0.465 * mm; 
  End1_steel_shell = new G4Tubs("End_1_steel_shell",End1r_min, End1r_max/2, End1length/2.,0.*deg,360.*deg);
  logical_End1_steel_shell = new G4LogicalVolume(End1_steel_shell, steelMat, "End1_steel_shell_log", 0, 0, 0);
  G4double end1offset_x = 0.0 * mm;
  G4double end1offset_y = 0.0 * mm;
  G4double end1offset_z = 2.0825 * mm;
  physical_End1_steel_shell = new G4PVPlacement(0,G4ThreeVector(end1offset_x,end1offset_y,end1offset_z),"phys_End1_steel_shell", logical_End1_steel_shell,mother, false, 0, true);

  G4double End1coner1_min = 0.0 * mm; 
  G4double End1coner1_max = 0.9 * mm;
  G4double End1coner2_min = 0.0 * mm;
  G4double End1coner2_max = 0.35 * mm;
  G4double End1conelength = 0.155 * mm; 
  End1cone_steel_shell = new G4Cons("End_1_cone_steel_shell",End1coner1_min, End1coner1_max/2,End1coner2_min, End1coner2_max/2, End1conelength/2.,0.*deg,360.*deg);
  logical_End1cone_steel_shell = new G4LogicalVolume(End1cone_steel_shell, steelMat, "End1cone_steel_shell_log", 0, 0, 0);
  G4double end1coneoffset_x = 0.0 * mm;
  G4double end1coneoffset_y = 0.0 * mm;
  G4double end1coneoffset_z = 2.3925 * mm;
  physical_End1cone_steel_shell = new G4PVPlacement(0,G4ThreeVector(end1coneoffset_x,end1coneoffset_y,end1coneoffset_z),"phys_End1cone_steel_shell", logical_End1cone_steel_shell,mother, false, 0, true);

//Define the cable weld end of the Steel shell 
  G4double End2r_min1 = 0.0 * mm;
  G4double End2r_max1 = 0.9 * mm;
  G4double End2r_min2 = 0.0 * mm;
  G4double End2r_max2 = 0.9 * mm;
  G4double End2length = 0.3 * mm;
  End2_steel_shell = new G4Cons("End_2_steel_shell",End2r_min2, End2r_max2/2, End2r_min1, End2r_max1/2, End2length/2.0, 0.0, 360.0*deg);
  logical_End2_steel_shell = new G4LogicalVolume(End2_steel_shell, steelMat, "End2_steel_shell_log", 0, 0, 0);
  G4double end2offset_x = 0.0 * mm;
  G4double end2offset_y = 0.0 * mm;
  G4double end2offset_z = -1.9 * mm; 
  physical_End2_steel_shell = new G4PVPlacement(0,G4ThreeVector(end2offset_x,end2offset_y,end2offset_z), "phys_End2_steel_shell", logical_End2_steel_shell,mother, false, 0, true);

//Define the cable 
  G4double cable_min = 0.0 * mm;
  G4double cable_max = 0.9 * mm;
  G4double cablelength = 5.0 * mm; 
  cable = new G4Tubs("cable",cable_min, cable_max/2, cablelength/2.,0.*deg,360.*deg);
  logical_cable = new G4LogicalVolume(cable, steelMat, "cable_log", 0, 0, 0);
  G4double cableoffset_x = 0.0 * mm;
  G4double cableoffset_y = 0.0 * mm;
  G4double cableoffset_z = -4.55 * mm;
  physical_cable = new G4PVPlacement(0,G4ThreeVector(cableoffset_x,cableoffset_y,cableoffset_z),"phys_cable", logical_cable, mother, false, 0, true);

// Define the Iridium core
  G4double corer_min = 0.0 * mm;	
  G4double corer_max = 0.6 * mm;
  G4double core_length = 3.5 * mm; 
  G4double iridiumcoreoffset_x = 0.0 * mm;
  G4double iridiumcoreoffset_y = 0.0 * mm;
  G4double iridiumcoreoffset_z = 0.0 * mm;
  iridium_core = new G4Tubs("iridium_core",corer_min, corer_max/2,core_length/2.,0.*deg,360.*deg);
  logical_iridium_core = new G4LogicalVolume(iridium_core, iridiumMat, "iridium_core_log", 0, 0, 0);
  physical_iridium_core = new G4PVPlacement(0,G4ThreeVector(iridiumcoreoffset_x,iridiumcoreoffset_y,iridiumcoreoffset_z), "phys_iridium_core", logical_iridium_core, physical_air_gap, false, 0, true);

// Define the first metal rod of applicator to do get correct measurements
  G4double rod1r_min = 0.0 * mm;	
  G4double rod1r_max = 3.0 * mm;
  G4double rod1_length = 60.0 * mm; 
  G4double rod1offset_x = 0.0 * mm;
  G4double rod1offset_y = 0.0 * mm;
  G4double rod1offset_z = 0.0 * mm;
  metal_rod1 = new G4Tubs("metal_rod1",rod1r_min, rod1r_max/2,rod1_length/2.,0.*deg,360.*deg);
  logical_metal_rod1 = new G4LogicalVolume(metal_rod1, titaniumMat, "metal_rod1_log", 0, 0, 0);
  G4RotationMatrix* rotationMatrixX90 = new G4RotationMatrix();
  rotationMatrixX90->rotateX(90.*deg);
  physical_metal_rod1 = new G4PVPlacement(0,G4ThreeVector(rod1offset_x,rod1offset_y,rod1offset_z), "phys_metal_rod1", logical_metal_rod1, mother, false, 0, true);//rot

// Define the end of first metal rod of applicator to do get correct measurements
  G4double rod1r_end_min = 0.0 * mm;	
  G4double rod1r_end_max = 3.0 * mm;
  G4double rod1_endoffset_x = 0.0 * mm;
  G4double rod1_endoffset_y = 0.0 * mm;
  G4double rod1_endoffset_z = -rod1_length/2.* mm;
  metal_rod1_end = new G4Sphere("metal_rod1_end",rod1r_end_min, rod1r_end_max/2,0.*deg,180.*deg,0.*deg,180.*deg);
  logical_metal_rod1_end = new G4LogicalVolume(metal_rod1_end, titaniumMat, "metal_rod1_end_log", 0, 0, 0);
  physical_metal_rod1_end = new G4PVPlacement(rotationMatrixX90,G4ThreeVector(rod1_endoffset_x,rod1_endoffset_y,rod1_endoffset_z), "phys_metal_rod1_end", logical_metal_rod1_end, mother, false, 0, true);//rot

// Define the air gap in the first metal rod
  G4double rod1r_air_min = 0.0 * mm;	
  G4double rod1r_air_max = 1.5 * mm;
  G4double rod1_air_length = 60.0 * mm; 
  G4double rod1offset_air_x = 0.0 * mm;
  G4double rod1offset_air_y = 0.0 * mm;
  G4double rod1offset_air_z = 0.0 * mm;
  air_rod1 = new G4Tubs("air_rod1",rod1r_air_min, rod1r_air_max/2,rod1_air_length/2.,0.*deg,360.*deg);
  logical_air_rod1 = new G4LogicalVolume(air_rod1, airMat, "air_rod1_log", 0, 0, 0);
  physical_air_rod1 = new G4PVPlacement(0,G4ThreeVector(rod1offset_air_x,rod1offset_air_y,rod1offset_air_z), "phys_air_rod1", logical_air_rod1, physical_metal_rod1, false, 0, true);

// Define the second metal rod of applicator to do get correct measurements
  G4double rod2r_min = 0.0 * mm;	
  G4double rod2r_max = 3.0 * mm;
  G4double rod2_length = 60.0 * mm; 
  G4double rod2offset_x = (16.*0.86602540378443871) * mm; //cos30, the 15 is half the diameter of the ring - change to variable name asap
  G4double rod2offset_y = 0.0 * mm;
  G4double rod2offset_z = 0.0 * mm;
  metal_rod2 = new G4Tubs("metal_rod2",rod2r_min, rod2r_max/2,rod2_length/2.,0.*deg,360.*deg);
  logical_metal_rod2 = new G4LogicalVolume(metal_rod2, titaniumMat, "metal_rod2_log", 0, 0, 0);
  physical_metal_rod2 = new G4PVPlacement(0,G4ThreeVector(rod2offset_x,rod2offset_y,rod2offset_z), "phys_metal_rod2", logical_metal_rod2, mother, false, 0, true);//rot

// Define the end of second metal rod of applicator to do get correct measurements
  G4double rod2r_end_min = 0.0 * mm;	
  G4double rod2r_end_max = 3.0 * mm;
  G4double rod2_endoffset_x = rod2offset_x * mm;
  G4double rod2_endoffset_y = 0.0 * mm;
  G4double rod2_endoffset_z = -rod2_length/2.* mm;
  metal_rod2_end = new G4Sphere("metal_rod2_end",rod2r_end_min, rod2r_end_max/2,0.*deg,180.*deg,0.*deg,180.*deg);
  logical_metal_rod2_end = new G4LogicalVolume(metal_rod2_end, titaniumMat, "metal_rod2_end_log", 0, 0, 0);
  physical_metal_rod2_end = new G4PVPlacement(rotationMatrixX90,G4ThreeVector(rod2_endoffset_x,rod2_endoffset_y,rod2_endoffset_z), "phys_metal_rod2_end", logical_metal_rod2_end, mother, false, 0, true);//rot

// Define the air gap in the second metal rod
  G4double rod2r_air_min = 0.0 * mm;	
  G4double rod2r_air_max = 1.5 * mm;
  G4double rod2_air_length = 60.0 * mm; 
  G4double rod2offset_air_x = 0.0 * mm;
  G4double rod2offset_air_y = 0.0 * mm;
  G4double rod2offset_air_z = 0.0 * mm;
  air_rod2 = new G4Tubs("air_rod2",rod2r_air_min, rod2r_air_max/2,rod2_air_length/2.,0.*deg,360.*deg);
  logical_air_rod2 = new G4LogicalVolume(air_rod2, airMat, "air_rod2_log", 0, 0, 0);
  physical_air_rod2 = new G4PVPlacement(0,G4ThreeVector(rod2offset_air_x,rod2offset_air_y,rod2offset_air_z), "phys_air_rod2", logical_air_rod2, physical_metal_rod2, false, 0, true);

  G4double rod2bent_length = 60.0 * mm; 
  G4double rod2bentoffset_x = rod2offset_x - (rod2bent_length*0.49999999999999994)/2.0 * mm;//sin30
  G4double rod2bentoffset_y = 0 * mm;//rod2_length/2. + (rod2bent_length*0.49999999999999994)/2 * mm;
  G4double rod2bentoffset_z = ((rod2bent_length)*0.86602540378443871)/2.+rod2_length/2. * mm;
  G4RotationMatrix* rotationMatrixY30 = new G4RotationMatrix();
  rotationMatrixY30->rotateY(30.*deg);
  metal_rod2bent = new G4Tubs("metal_rod2bent",rod2r_min, rod2r_max/2,rod2bent_length/2.,0.*deg,360.*deg);
  logical_metal_rod2bent = new G4LogicalVolume(metal_rod2bent, titaniumMat, "metal_rod2bent_log", 0, 0, 0);
  physical_metal_rod2bent = new G4PVPlacement(rotationMatrixY30,G4ThreeVector(rod2bentoffset_x,rod2bentoffset_y,rod2bentoffset_z), "phys_metal_rod2bent", logical_metal_rod2bent, mother, false, 0, true);

// Define the end of bent second metal rod of applicator to do get correct measurements
  G4double rod2bentr_end_min = 0.0 * mm;	
  G4double rod2bentr_end_max = 3.0 * mm;
  G4double rod2bent_endoffset_x = rod2bentoffset_x * mm;
  G4double rod2bent_endoffset_y = 0.0 * mm;
  G4double rod2bent_endoffset_z = rod2bentoffset_z + 0.86602540378443871*rod2bent_length/2.* mm;
  G4RotationMatrix* rotationMatrixX90Y30 = new G4RotationMatrix();
  rotationMatrixX90Y30->rotateX(-90.*deg);
  rotationMatrixX90Y30->rotateZ(-30.*deg); //for some reason once I've flipped the X axis 90 degrees a rotation around Y must be written as a negative rotation around Z??? Go figure
  G4ThreeVector rod2bentr_end_3vec(rod2bent_endoffset_x-0.5*rod2bent_length/2,rod2bent_endoffset_y,rod2bent_endoffset_z);
  //rod2bentr_end_3vec.rotateX(-90.*deg);
  //rod2bentr_end_3vec.rotateY(30.*deg);
  metal_rod2bent_end = new G4Sphere("metal_rod2bent_end",rod2bentr_end_min, rod2bentr_end_max/2,0.*deg,180.*deg,0.*deg,180.*deg);
  logical_metal_rod2bent_end = new G4LogicalVolume(metal_rod2bent_end, titaniumMat, "metal_rod2bent_end_log", 0, 0, 0);
  physical_metal_rod2bent_end = new G4PVPlacement(rotationMatrixX90Y30,rod2bentr_end_3vec, "phys_metal_rod2bent_end", logical_metal_rod2bent_end, mother, false, 0, true);//rot

// Define the air gap in the bent metal rod
  G4double rod3r_air_min = 0.0 * mm;	
  G4double rod3r_air_max = 1.5 * mm;
  G4double rod3_air_length = 60.0 * mm; 
  G4double rod3offset_air_x = 0.0 * mm;
  G4double rod3offset_air_y = 0.0 * mm;
  G4double rod3offset_air_z = 0.0 * mm;
  air_rod3 = new G4Tubs("air_rod3",rod3r_air_min, rod3r_air_max/2,rod3_air_length/2.,0.*deg,360.*deg);
  logical_air_rod3 = new G4LogicalVolume(air_rod3, airMat, "air_rod3_log", 0, 0, 0);
  physical_air_rod3 = new G4PVPlacement(0,G4ThreeVector(rod3offset_air_x,rod3offset_air_y,rod3offset_air_z), "phys_air_rod3", logical_air_rod3, physical_metal_rod2bent, false, 0, true);

// Define the plastic ring of applicator to do get correct measurements
  G4double plas_ringr_min = 0.0 * mm;	
  G4double plas_ringr_max = 7.5 * mm;
  G4double plas_ring_diameter = 32.0 * mm; 
  //G4double plas_ringphi_min = 0. *deg;
  //G4double plas_ringphi_max = 360. *deg;
  G4double plas_ringoffset_x = (plas_ring_diameter/2*0.86602540378443871) * mm;// 0.5 = sin30
  G4double plas_ringoffset_y = 0.0 * mm;
  G4double plas_ringoffset_z = rod1_length/2. + (plas_ring_diameter*0.5)/2 * mm;
  G4RotationMatrix* rotationMatrixYminus60 = new G4RotationMatrix();
  rotationMatrixYminus60->rotateY(30.*deg);
  plas_ring = new G4Torus("plas_ring",plas_ringr_min, plas_ringr_max/2,plas_ring_diameter/2,0.*deg,360.*deg);
  logical_plas_ring = new G4LogicalVolume(plas_ring, acetalMat, "plas_ring_log", 0, 0, 0);
  physical_plas_ring = new G4PVPlacement(rotationMatrixYminus60,G4ThreeVector(plas_ringoffset_x,plas_ringoffset_y,plas_ringoffset_z), "phys_plas_ring", logical_plas_ring, mother, false, 0, true);

// Define the metal ring of applicator to do get correct measurements
  G4double ringr_min = 0.0 * mm;	
  G4double ringr_max = 3.0 * mm;
  G4double ring_diameter = 32.0 * mm; 
  //G4double ringphi_min = 0. *deg;
  //G4double ringphi_max = 360. *deg;
  G4double ringoffset_x = 0.0 * mm;// 0.5 = sin30
  G4double ringoffset_y = 0.0 * mm;
  G4double ringoffset_z = 0.0 * mm;
  metal_ring = new G4Torus("metal_ring",ringr_min, ringr_max/2,ring_diameter/2,0.*deg,360.*deg);
  logical_metal_ring = new G4LogicalVolume(metal_ring, titaniumMat, "metal_ring_log", 0, 0, 0);
  physical_metal_ring = new G4PVPlacement(0,G4ThreeVector(ringoffset_x,ringoffset_y,ringoffset_z), "phys_metal_ring", logical_metal_ring, physical_plas_ring, false, 0, true); //Geant is clever - if this ring lives inside the plastic one then it automatically take on the rotation that I have applied to the plastic ring

// Define the air ring of applicator to do get correct measurements
  G4double air_ringr_min = 0.0 * mm;	
  G4double air_ringr_max = 1.5 * mm;
  G4double air_ring_diameter = 32.0 * mm; 
  //G4double ringphi_min = 0. *deg;
  //G4double ringphi_max = 360. *deg;
  G4double air_ringoffset_x = 0.0 * mm;// 0.5 = sin30
  G4double air_ringoffset_y = 0.0 * mm;
  G4double air_ringoffset_z = 0.0 * mm;
  air_ring = new G4Torus("air_ring",air_ringr_min, air_ringr_max/2,air_ring_diameter/2,0.*deg,360.*deg);
  logical_air_ring = new G4LogicalVolume(air_ring, airMat, "air_ring_log", 0, 0, 0);
  physical_air_ring = new G4PVPlacement(0,G4ThreeVector(air_ringoffset_x,air_ringoffset_y,air_ringoffset_z), "phys_air_ring", logical_air_ring, physical_metal_ring, false, 0, true); //Geant is clever - if this ring lives inside the plastic one then it automatically take on the rotation that I have applied to the plastic ring


// Visualisations

//Shell/cable attributes    
  steelAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  steelAttributes -> SetVisibility(true);
  steelAttributes -> SetForceAuxEdgeVisible(true);//SetForceWireframe(true);

  endAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  endAttributes -> SetVisibility(true);
  endAttributes -> SetForceAuxEdgeVisible(true);//SetForceWireframe(true);//SetForceAuxEdgeVisible(true);
  logical_steel_shell -> SetVisAttributes(steelAttributes);
  logical_End1_steel_shell -> SetVisAttributes(endAttributes);
  logical_End1cone_steel_shell -> SetVisAttributes(endAttributes);
  logical_End2_steel_shell -> SetVisAttributes(endAttributes);
  logical_cable -> SetVisAttributes(steelAttributes);

  titaniumAttributes = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); // red
  titaniumAttributes -> SetVisibility(true);
  titaniumAttributes -> SetForceAuxEdgeVisible(true);//SetForceWireframe(true);//SetForceAuxEdgeVisible(true);

  logical_metal_rod1 -> SetVisAttributes(titaniumAttributes);
  logical_metal_rod2 -> SetVisAttributes(titaniumAttributes);
  logical_metal_rod2bent -> SetVisAttributes(titaniumAttributes);
  logical_metal_ring -> SetVisAttributes(titaniumAttributes);
  logical_metal_rod1_end -> SetVisAttributes(titaniumAttributes);
  logical_metal_rod2_end -> SetVisAttributes(titaniumAttributes);
  logical_metal_rod2bent_end -> SetVisAttributes(titaniumAttributes);
 
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour white (1.0, 1.0, 1.0) ;

  simpleIridiumVisAtt = new G4VisAttributes(magenta);
  simpleIridiumVisAtt -> SetVisibility(true);
  simpleIridiumVisAtt -> SetForceWireframe(true);
  logical_iridium_core -> SetVisAttributes(simpleIridiumVisAtt);

  acetalAttributes = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); 
  acetalAttributes -> SetVisibility(true);
  acetalAttributes -> SetForceAuxEdgeVisible(true);//SetForceWireframe(true);//SetForceAuxEdgeVisible(true);
  logical_plas_ring -> SetVisAttributes(acetalAttributes);

  airAttributes = new G4VisAttributes(white);
  airAttributes -> SetVisibility(true);
  airAttributes -> SetForceAuxEdgeVisible(true);//SetForceWireframe(true);//
  logical_air_rod1 -> SetVisAttributes(airAttributes);
  logical_air_rod2 -> SetVisAttributes(airAttributes);
  logical_air_rod3 -> SetVisAttributes(airAttributes);
  logical_air_ring -> SetVisAttributes(airAttributes);

}

void BrachyDetectorConstructionGammaMed::CleanGammaMed()
{ 
  delete simpleIridiumVisAtt; 
  simpleIridiumVisAtt = 0;
  
  delete endAttributes; 
  endAttributes = 0;
  
  delete steelAttributes; 
  steelAttributes = 0;
  
  delete physical_iridium_core; 
  physical_iridium_core = 0 ;

  delete logical_iridium_core; 
  logical_iridium_core = 0;
  
  delete iridium_core; 
  iridium_core = 0;
  
  delete physical_cable;
  physical_cable = 0;
 
  delete logical_cable; 
  logical_cable = 0;

  delete cable; 
  cable = 0;
  
  delete physical_metal_ring;
  physical_metal_ring = 0;
 
  delete logical_metal_ring; 
  logical_metal_ring = 0;

  delete metal_ring; 
  metal_ring = 0;
  
  delete physical_plas_ring;
  physical_plas_ring = 0;
 
  delete logical_plas_ring; 
  logical_plas_ring = 0;

  delete plas_ring; 
  plas_ring = 0;
  
  delete physical_metal_rod1;
  physical_metal_rod1 = 0;
 
  delete logical_metal_rod1; 
  logical_metal_rod1 = 0;

  delete metal_rod1; 
  metal_rod1 = 0;
  
  delete physical_metal_rod1_end;
  physical_metal_rod1_end = 0;
 
  delete logical_metal_rod1_end; 
  logical_metal_rod1_end = 0;

  delete metal_rod1_end; 
  metal_rod1_end = 0;
  
  delete physical_metal_rod2_end;
  physical_metal_rod2_end = 0;
 
  delete logical_metal_rod2_end; 
  logical_metal_rod2_end = 0;

  delete metal_rod2_end; 
  metal_rod2_end = 0;
  
  delete physical_metal_rod2bent_end;
  physical_metal_rod2bent_end = 0;
 
  delete logical_metal_rod2bent_end; 
  logical_metal_rod2bent_end = 0;

  delete metal_rod2bent_end; 
  metal_rod2bent_end = 0;
  
  delete physical_air_rod1;
  physical_air_rod1 = 0;
 
  delete logical_air_rod1; 
  logical_air_rod1 = 0;

  delete air_rod1; 
  air_rod1 = 0;
  
  delete physical_metal_rod2;
  physical_metal_rod2 = 0;
 
  delete logical_metal_rod2; 
  logical_metal_rod2 = 0;

  delete metal_rod2; 
  metal_rod2 = 0;
  
  delete physical_air_rod2;
  physical_air_rod2 = 0;
 
  delete logical_air_rod2; 
  logical_air_rod2 = 0;

  delete air_rod2; 
  air_rod2 = 0;
  
  delete physical_metal_rod2bent;
  physical_metal_rod2bent = 0;
 
  delete logical_metal_rod2bent; 
  logical_metal_rod2bent = 0;

  delete metal_rod2bent; 
  metal_rod2bent = 0;
  
  delete physical_air_rod3;
  physical_air_rod3 = 0;
 
  delete logical_air_rod3; 
  logical_air_rod3 = 0;

  delete air_rod3; 
  air_rod3 = 0;
  
  delete physical_air_ring;
  physical_air_ring = 0;
 
  delete logical_air_ring; 
  logical_air_ring = 0;

  delete air_ring; 
  air_ring = 0;
  
  delete physical_End2_steel_shell; 
  physical_End2_steel_shell = 0;
   
  delete logical_End2_steel_shell; 
  logical_End2_steel_shell = 0;
  
  delete End2_steel_shell; 
  End2_steel_shell = 0;
  
  delete physical_End1_steel_shell; 
  physical_End1_steel_shell = 0;
   
  delete physical_End1cone_steel_shell; 
  physical_End1cone_steel_shell = 0;
   
  delete logical_End1_steel_shell; 
  logical_End1_steel_shell = 0;
  
  delete End1_steel_shell; 
  End1_steel_shell = 0;

  delete logical_End1cone_steel_shell; 
  logical_End1cone_steel_shell = 0;
  
  delete End1cone_steel_shell; 
  End1cone_steel_shell = 0;

  delete physical_air_gap;
  physical_air_gap = 0;

  delete logical_air_gap;
  logical_air_gap = 0;

  delete air_gap;
  air_gap = 0;

  delete physical_steel_shell;
  physical_steel_shell = 0;

  delete logical_steel_shell;
  logical_steel_shell = 0;
 
  delete steel_shell;
  steel_shell = 0;
 
 G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
