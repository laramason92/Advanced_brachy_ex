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
// $Id: SteppingAction.cc,v 1.10 2006/06/29 16:24:25 gunter Exp $
// GEANT4 tag $Name: geant4-09-02-ref-04 $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//

#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ParticleDefinition.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"

#ifdef ANALYSIS_USE
#include "BrachyAnalysisManager.hh"
#endif

#include "BrachySteppingAction.hh"
#include "G4SystemOfUnits.hh"

BrachySteppingAction::BrachySteppingAction()
{ 

}

BrachySteppingAction::~BrachySteppingAction()
{ 
}

void BrachySteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4SteppingManager*  steppingManager = fpSteppingManager;
  G4Track* theTrack = aStep-> GetTrack();


//  G4ThreeVector pos = theTrack->GetVertexPosition()/mm;
//  G4StepPoint* presteppt = aStep->GetPreStepPoint();

//  if ( (theTrack->GetParticleDefinition()->GetPDGCharge()!=0) && (theTrack->GetParentID()==1)){//is the charged child of a primary photon
          
//          G4String materialprestep = presteppt->GetMaterial()->GetName();

//          G4String materialvertex = theTrack->GetLogicalVolumeAtVertex()->GetMaterial()->GetName();

//          if ( (materialvertex == "DryAir")){

//		  G4double eKinVertex = theTrack->GetVertexKineticEnergy()/keV; // Ekin at vertex, divide by MeV to put the units as MeV
//		  G4double mass_g = presteppt->GetPhysicalVolume()->GetLogicalVolume()->GetMass()/g; //gives us the mass of the critical volume (10cm^3)

//		  G4double kerma = eKinVertex / mass_g;
	          
//                  if (materialvertex != materialprestep){
//                  G4cout << "materials dont match" << G4endl; 
//                  }

//		  G4double xpos_kerma = pos.x()/mm;
//		  G4double ypos_kerma = pos.y()/mm;
//		  G4double zpos_kerma = pos.z()/mm;

		  //G4cout << "particle id - " << id << G4endl;
		  //G4cout << "parent id - " << parentid << G4endl;
		  //G4cout << "charge - " << charge << G4endl;
		  //G4cout << "vertex position - " << pos << G4endl;
		  //G4cout <<  "vertex pos x - " << xpos_kerma << G4endl;
		  //G4cout <<  "material of prestep point- " << materialprestep << G4endl;   
		  // G4cout <<  "post step point - " << coordpost << G4endl;
		  //G4cout << "pre step point - " << coordpre << G4endl;
		  //G4cout << "mass - " << mass_g << G4endl;
		  //G4cout << "****************************" << G4endl; 



//#ifdef ANALYSIS_USE  
//            BrachyAnalysisManager* analysis = BrachyAnalysisManager::GetInstance();
//            if(zpos_kerma >= -50.0 *mm && zpos_kerma <= 50.0*mm) {analysis -> FillH3WithKerma(xpos_kerma,ypos_kerma,kerma);
//                        }
//            if(zpos_kerma >= -50.0 *mm && zpos_kerma <= 50.0*mm) {analysis -> FillH5WithKerma(xpos_kerma,ypos_kerma,kerma);
//             }                       
        
//#endif

//		  theTrack -> SetTrackStatus(fKillTrackAndSecondaries);
      
 
    //}
  //}   


  




  //if(theTrack-> GetTrackStatus() == fAlive) {return;}

  // G4cout << "Start secondariessss" << G4endl;  
  // Retrieve the secondary particles
  G4TrackVector* fSecondary = steppingManager -> GetfSecondary();
     
   for(size_t lp1=0;lp1<(*fSecondary).size(); lp1++)
   { 
     // Retrieve particle
     const G4ParticleDefinition* particleName = (*fSecondary)[lp1] -> GetDefinition();      

     if (particleName == G4Gamma::Definition())
     {
      G4String process = (*fSecondary)[lp1]-> GetCreatorProcess()-> GetProcessName();  
      
      // Retrieve the process originating it
      //G4cout << "creator process " << process << G4endl;
        if (process == "RadioactiveDecay")
         {
#ifdef ANALYSIS_USE  
          BrachyAnalysisManager* analysis = BrachyAnalysisManager::GetInstance();
          G4double energy = (*fSecondary)[lp1]  -> GetKineticEnergy();
          // Store the initial energy of particles in a 1D histogram
          analysis -> FillPrimaryParticleHistogram(energy/keV);
#endif

         }
      }
   } 
}
