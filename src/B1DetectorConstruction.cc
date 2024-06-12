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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
     
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

//     
// World
//  
G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
G4Box* solidWorld =    
    new G4Box("WorldBox",                       //its name
      3*m, 3*m,3*m);     //its size
      
G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "LogicalWorld");            //its name
                                   
G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                      
//Moon surface
G4Material* Oxygen = nist->FindOrBuildMaterial("G4_O");
G4Material* Sodium = nist->FindOrBuildMaterial("G4_Na");
G4Material* Magnesium = nist->FindOrBuildMaterial("G4_Mg");
G4Material* Aluminum = nist->FindOrBuildMaterial("G4_Al");
G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");       
G4Material* Calcium = nist->FindOrBuildMaterial("G4_Ca");       
G4Material* Titanium = nist->FindOrBuildMaterial("G4_Ti");       
G4Material* Manganese = nist->FindOrBuildMaterial("G4_Mn");       
G4Material* Iron = nist->FindOrBuildMaterial("G4_Fe");
G4Material* Hydrogen = nist->FindOrBuildMaterial("G4_H");                             
G4String name;
G4int ncomponents;
G4double fracMass;
G4Material* Layer_I_mat = new G4Material(name = "Layer_I_mat", 0.85*g/cm3, ncomponents = 10);
Layer_I_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_I_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_I_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_I_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_I_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_I_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_I_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_I_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_I_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_I_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_II_mat = new G4Material(name = "Layer_II_mat", 1.2*g/cm3, ncomponents = 10);
Layer_II_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_II_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_II_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_II_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_II_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_II_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_II_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_II_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_II_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_II_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_III_mat = new G4Material(name = "Layer_III_mat", 1.5*g/cm3, ncomponents = 10);
Layer_III_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_III_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_III_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_III_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_III_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_III_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_III_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_III_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_III_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_III_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_IV_mat = new G4Material(name = "Layer_IV_mat", 1.6*g/cm3, ncomponents = 10);
Layer_IV_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_IV_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_IV_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_IV_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_IV_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_IV_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_IV_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_IV_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_IV_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_IV_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_V_mat = new G4Material(name = "Layer_V_mat", 1.7*g/cm3, ncomponents = 10);
Layer_V_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_V_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_V_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_V_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_V_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_V_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_V_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_V_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_V_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_V_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_VI_mat = new G4Material(name = "Layer_VI_mat", 1.8*g/cm3, ncomponents = 10);
Layer_VI_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_VI_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_VI_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_VI_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_VI_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_VI_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_VI_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_VI_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_VI_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_VI_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
G4Material* Layer_VII_mat = new G4Material(name = "Layer_VII_mat", 1.9*g/cm3, ncomponents = 10);
Layer_VII_mat->AddMaterial(Oxygen, fracMass = 43.6 * perCent);
Layer_VII_mat->AddMaterial(Hydrogen, fracMass = 0.1 * perCent);
Layer_VII_mat->AddMaterial(Sodium, fracMass = 0.3 * perCent);
Layer_VII_mat->AddMaterial(Magnesium, fracMass = 5.6 * perCent);         
Layer_VII_mat->AddMaterial(Aluminum, fracMass = 9 * perCent);         
Layer_VII_mat->AddMaterial(Silicon, fracMass = 21.1 * perCent);         
Layer_VII_mat->AddMaterial(Calcium, fracMass = 8.5 * perCent);         
Layer_VII_mat->AddMaterial(Titanium, fracMass = 1.5 * perCent);
Layer_VII_mat->AddMaterial(Manganese, fracMass = 0.1 * perCent);
Layer_VII_mat->AddMaterial(Iron, fracMass = 10.2 * perCent);
//66-cm Layer_VII
G4Box* Layer_box =    
    new G4Box("Layer_box",                   
      1.5*m, 1*cm,1.5*m); 

G4LogicalVolume* Layer_VII =                         
    new G4LogicalVolume(Layer_box,        
                        Layer_VII_mat,         
                        "Layer_VII_logic");         
                                   
for (int i = 1; i < 34; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(-65+(2*i))*cm, 0*m), Layer_VII, "Layer_VII", logicWorld, false, i, checkOverlaps);                                    
  }  
//64-cm Layer_VI  
G4LogicalVolume* Layer_VI =                         
    new G4LogicalVolume(Layer_box,        
                        Layer_VI_mat,           
                        "Layer_VI_logic");        
for (int i = 1; i < 33; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(1+(2*i))*cm, 0*m), Layer_VI, "Layer_VI", logicWorld, false, 33+i, checkOverlaps);                                   
  }
//28-cm Layer_V
G4LogicalVolume* Layer_V =                         
    new G4LogicalVolume(Layer_box,        
                        Layer_V_mat,        
                        "Layer_V_logic");        
                                   
for (int i = 1; i < 15; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(65+(2*i))*cm, 0*m), Layer_V, "Layer_V", logicWorld, false, 65+i, checkOverlaps);                                    
  }                           

//26-cm Layer_IV    
G4LogicalVolume* Layer_IV =                         
    new G4LogicalVolume(Layer_box,        
                        Layer_IV_mat,         
                        "Layer_IV_logic");  
                                   
for (int i = 1; i < 14; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(93+(2*i))*cm, 0*m), Layer_IV, "Layer_IV", logicWorld, false, 79+i, checkOverlaps);                                  
  }
//6-cm Layer_III    
G4LogicalVolume* Layer_III =                         
    new G4LogicalVolume(Layer_box,         
                        Layer_III_mat,          
                        "Layer_III_logic");         
                                   
for (int i = 1; i < 4; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(119+(2*i))*cm, 0*m), Layer_III, "Layer_III", logicWorld, false, 92+i, checkOverlaps);                                  
  }
//6-cm Layer_II    
G4LogicalVolume* Layer_II =                         
    new G4LogicalVolume(Layer_box,        
                        Layer_II_mat,           
                        "Layer_II_logic");      
                                   
for (int i = 1; i < 4; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(125+(2*i))*cm, 0*m), Layer_II, "Layer_II", logicWorld, false, 95+i, checkOverlaps);                                  
  }
//4-cm Layer_I    
G4LogicalVolume* Layer_I =                         
    new G4LogicalVolume(Layer_box,         
                        Layer_I_mat,          
                        "Layer_I_logic");       
                                   
for (int i = 1; i < 3; i++){
                         new G4PVPlacement(0, G4ThreeVector(0*m,(131+(2*i))*cm, 0*m), Layer_I, "Layer_I", logicWorld, false, 98+i, checkOverlaps);                                     
  }      
  
//Pseudo detector
G4ThreeVector Pos_surface = G4ThreeVector(0*m, 1.361*m, 0*m);
G4Box* Surface_detector_box =    
    new G4Box("Surface_detector_box",                      
      1.5*m, 1*mm,1.5*m);     
      
G4LogicalVolume* Surface_detector =                         
    new G4LogicalVolume(Surface_detector_box,         
                        world_mat,          
                        "Surface_detector_logic");         
                                   
    new G4PVPlacement(0,                    
                      Pos_surface,      
                      Surface_detector,           
                      "Surface_detector",               
                      logicWorld,                     
                      false,                 
                      0,                     
                      checkOverlaps);
G4VisAttributes* visAttributesLayer_I = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
visAttributesLayer_I->SetVisibility(true);
Layer_I->SetVisAttributes(visAttributesLayer_I);
G4VisAttributes* visAttributesLayer_II = new G4VisAttributes(G4Colour(0.68, 0.85, 0.9));
visAttributesLayer_II->SetVisibility(true);
Layer_II->SetVisAttributes(visAttributesLayer_II);
G4VisAttributes* visAttributesLayer_III = new G4VisAttributes(G4Colour(0.0, 0.0, 0.55));
visAttributesLayer_III->SetVisibility(true);
Layer_III->SetVisAttributes(visAttributesLayer_III);
G4VisAttributes* visAttributesLayer_IV = new G4VisAttributes(G4Colour(1.0, 0.6, 0.6)); 
visAttributesLayer_IV->SetVisibility(true);
Layer_IV->SetVisAttributes(visAttributesLayer_IV);
G4VisAttributes* visAttributesLayer_V = new G4VisAttributes(G4Colour(0.55, 0.0, 0.0));
visAttributesLayer_V->SetVisibility(true);
Layer_V->SetVisAttributes(visAttributesLayer_V);
G4VisAttributes* visAttributesLayer_VI = new G4VisAttributes(G4Colour(0.82, 0.71, 0.55));
visAttributesLayer_VI->SetVisibility(true);
Layer_VI->SetVisAttributes(visAttributesLayer_VI);
G4VisAttributes* visAttributesLayer_VII = new G4VisAttributes(G4Colour(0.4, 0.2, 0.1));
Layer_VII->SetVisAttributes(visAttributesLayer_VII);                                                   
  //always return the physical World
  //  
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
