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
/// \file B4/B4e/src/DetectorConstruction.cc
/// \brief Implementation of the B4e::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "CellParameterisation.hh"
#include "CalorimeterSD.hh"
#include "Constants.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

namespace B4e
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;
  G4double density;
  new G4Material("liquidArgon", z = 18., a = 39.95 * g / mole, density = 1.390 * g / cm3);
  // The argon by NIST Manager is a gas with a different density
  G4Material* BGO = G4NistManager::Instance()->FindOrBuildMaterial("G4_BGO");

  G4cout << "Defined BGO material: " << *BGO << G4endl;

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  fNofLayers = 14;
  G4double absoThickness = 25. * mm;
  G4double gapThickness = 0.5 * mm;
  G4double calorSizeXY = 60. * cm;

  G4double barSizeX = 25. * mm;
  G4double barSizeY = 600. * mm;
  G4int fNofBars = 24;  // 每层 24 根 bars

  auto layerThickness = absoThickness + gapThickness;
  auto calorThickness = fNofLayers * layerThickness;
  SetCalorThickness(calorThickness);
  SetLayerThickness(layerThickness);
  auto worldSizeXY = 1.2 * calorSizeXY;
  auto worldSizeZ = 1.2 * calorThickness;

  // Get materials
  auto defaultMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  auto absorberMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_BGO");
  auto gapMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  if (!defaultMaterial || !absorberMaterial || !gapMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS = new G4Box("World",  // its name
                          worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);  // its size

  auto worldLV = new G4LogicalVolume(worldS,  // its solid
                                     defaultMaterial,  // its material
                                     "World");  // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
                                   G4ThreeVector(),  // at (0,0,0)
                                   worldLV,  // its logical volume
                                   "World",  // its name
                                   nullptr,  // its mother  volumezc
                                   false,  // no boolean operation
                                   0,  // copy number
                                   fCheckOverlaps);  // checking overlaps

  //
  // Calorimeter
  //
  auto calorimeterS = new G4Box("Calorimeter",  // its name
                                calorSizeXY / 2, calorSizeXY / 2, calorThickness / 2);  // its size

  auto calorLV = new G4LogicalVolume(calorimeterS,  // its solid
                                     defaultMaterial,  // its material
                                     "Calorimeter");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    calorLV,  // its logical volume
                    "Calorimeter",  // its name
                    worldLV,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    fCheckOverlaps);  // checking overlap
  //
  // Absorber
  //
  auto absorberS = new G4Box("Abso",  // its name
                barSizeX / 2, barSizeY / 2, absoThickness / 2);  // its size

  auto absorberLV = new G4LogicalVolume(absorberS,  // its solid
                                        absorberMaterial,  // its material
                                        "AbsoLV");  // its name

  G4VPVParameterisation* cellParam = new CellParameterisation();

  fAbsorberPV = new G4PVParameterised("AbsoLV",    // its name
                                      absorberLV,  // its logical volume
                                      calorLV, 
                                      kXAxis, 
                                      kNofEmCells,
                                      cellParam);
  //
  // print parameters
  //
  G4cout << G4endl << "------------------------------------------------------------" << G4endl
         << "---> The calorimeter is " << fNofLayers << " layers of: [ " << absoThickness / mm
         << "mm of " << absorberMaterial->GetName() << " + " << gapThickness / mm << "mm of "
         << gapMaterial->GetName() << " ] " << G4endl
         << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes(G4Colour::White()));
  absorberLV->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));
  // gapLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
  auto absoSD = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", kNofEmCells);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV", absoSD);

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4e
