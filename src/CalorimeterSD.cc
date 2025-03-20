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
/// \file B4/B4/src/CalorimeterSD.cc
/// \brief Implementation of the B4::CalorimeterSD class

#include "CalorimeterSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorimeterSD::CalorimeterSD(const G4String& name, const G4String& hitsCollectionName, G4int nofCells)
  : G4VSensitiveDetector(name), fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new CalorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);

  // Create hits
  // fNofCells for cells + one more for total sums
  for (G4int i = 0; i < fNofCells + 14 + 1; i++) {
    fHitsCollection->insert(new CalorHit());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool CalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  G4int fireTag = 0;

  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.) {
    stepLength = step->GetStepLength();
  }
  if (edep == 0. && stepLength == 0.) return false;
  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Get calorimeter cell id
  auto cellID = touchable->GetCopyNumber();

  // Get hit accounting data for this cell
  auto hit = (*fHitsCollection)[cellID];
  if (!hit) {
    G4ExceptionDescription msg;
    msg << "Cannot access hit " << cellID;
    G4Exception("CalorimeterSD::ProcessHits()", "MyCode0004", FatalException, msg);
  }
    // Add values
  hit->Add(edep, stepLength); 
  // Get hit for total accounting
  auto hitTotal = (*fHitsCollection)[fHitsCollection->entries() - 1];   hitTotal->Add(edep, stepLength);
  // Get hit for total accounting
  auto hitLayer0  = (*fHitsCollection)[fHitsCollection->entries() - 15]; if(cellID / 24 ==0  ) {hitLayer0 ->Add(edep, stepLength); }
  auto hitLayer1  = (*fHitsCollection)[fHitsCollection->entries() - 14]; if(cellID / 24 ==1  ) {hitLayer1 ->Add(edep, stepLength); }
  auto hitLayer2  = (*fHitsCollection)[fHitsCollection->entries() - 13]; if(cellID / 24 ==2  ) {hitLayer2 ->Add(edep, stepLength); }
  auto hitLayer3  = (*fHitsCollection)[fHitsCollection->entries() - 12]; if(cellID / 24 ==3  ) {hitLayer3 ->Add(edep, stepLength); }
  auto hitLayer4  = (*fHitsCollection)[fHitsCollection->entries() - 11]; if(cellID / 24 ==4  ) {hitLayer4 ->Add(edep, stepLength); }
  auto hitLayer5  = (*fHitsCollection)[fHitsCollection->entries() - 10]; if(cellID / 24 ==5  ) {hitLayer5 ->Add(edep, stepLength); }
  auto hitLayer6  = (*fHitsCollection)[fHitsCollection->entries() -  9]; if(cellID / 24 ==6  ) {hitLayer6 ->Add(edep, stepLength); }
  auto hitLayer7  = (*fHitsCollection)[fHitsCollection->entries() -  8]; if(cellID / 24 ==7  ) {hitLayer7 ->Add(edep, stepLength); }
  auto hitLayer8  = (*fHitsCollection)[fHitsCollection->entries() -  7]; if(cellID / 24 ==8  ) {hitLayer8 ->Add(edep, stepLength); }
  auto hitLayer9  = (*fHitsCollection)[fHitsCollection->entries() -  6]; if(cellID / 24 ==9  ) {hitLayer9 ->Add(edep, stepLength); }
  auto hitLayer10 = (*fHitsCollection)[fHitsCollection->entries() -  5]; if(cellID / 24 ==10 ) {hitLayer10->Add(edep, stepLength); }
  auto hitLayer11 = (*fHitsCollection)[fHitsCollection->entries() -  4]; if(cellID / 24 ==11 ) {hitLayer11->Add(edep, stepLength); }
  auto hitLayer12 = (*fHitsCollection)[fHitsCollection->entries() -  3]; if(cellID / 24 ==12 ) {hitLayer12->Add(edep, stepLength); }
  auto hitLayer13 = (*fHitsCollection)[fHitsCollection->entries() -  2]; if(cellID / 24 ==13 ) {hitLayer13->Add(edep, stepLength); }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    auto nofHits = fHitsCollection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << nofHits
           << " hits in the tracker chambers: " << G4endl;
    for (std::size_t i = 0; i < nofHits; ++i)
      (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4
