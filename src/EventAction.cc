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
/// \file B4/B4c/src/EventAction.cc
/// \brief Implementation of the B4c::EventAction class

#include "EventAction.hh"

#include "RunAction.hh"

#include "CalorHit.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include "G4RunManager.hh"
#include <iomanip>

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection* EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  auto hitsCollection = static_cast<CalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  if (!hitsCollection) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                                       G4double gapEdep, G4double gapTrackLength) const
{
  // print event statistics
  G4cout << "   Absorber: total energy: " << std::setw(7) << G4BestUnit(absoEdep, "Energy")
         << "       total track length: " << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
         << G4endl << "        Gap: total energy: " << std::setw(7) << G4BestUnit(gapEdep, "Energy")
         << "       total track length: " << std::setw(7) << G4BestUnit(gapTrackLength, "Length")
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/) 
{
    fFirstInteraction = G4ThreeVector(-1e9, -1e9, -1e9);  
    fInteractionDepth = -1;
    fInteractionLayer = -1;
    fSecondaries = -1;
    fInteractionType = -1;
    fHadrInteractionDepth = -1;
    fHadrInteractionLayer = -1;
    fHadrSecondaries = -1;
    fHadronicTag=-1; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Get hits collections IDs (only once)
  if (fAbsHCID == -1) {
    fAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitsCollection");
    fGapHCID = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }

  // Get hits collections
  auto absoHC = GetHitsCollection(fAbsHCID, event);
  auto gapHC  = GetHitsCollection(fGapHCID, event);

  // Get primary particle information
  const G4PrimaryParticle* primary = event->GetPrimaryVertex(0)->GetPrimary();  
  G4String particleName = primary->GetParticleDefinition()->GetParticleName();
  G4int particleID = 0;
  if (particleName == "proton") {particleID = 1; }
  else if (particleName == "deuteron") {particleID = 2; } 
  else if (particleName == "He3") { particleID = 3; } 
  else if (particleName == "alpha") { particleID = 4; } 
  else if (particleName == "C12") { particleID = 6; } 
  else if (particleName == "e-") { particleID = 0; } 
  else {particleID = -1;  } 

  G4double energy = primary->GetKineticEnergy();
  // G4cout << "!! particleName: " << particleName << " ,Kinetic Energy = " << energy << "MeV" << G4endl;

  // Print per event (modulo n)
  auto eventID = event->GetEventID();

  G4RunManager* rm = G4RunManager::GetRunManager(); 
  auto printModulo = rm->GetPrintProgress();
  
  if ((printModulo > 0) && (eventID % printModulo == 0)) {
    G4cout << "--> End of event: " << eventID << "\n" << G4endl;
  }

  // Fill histograms, ntuple
  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleIColumn(0, particleID);
  analysisManager->FillNtupleDColumn(1, energy / CLHEP::GeV);

  for (int i = 0; i < absoHC->entries(); ++i) 
  {
    auto absoperHit = (*absoHC)[i];
    if (absoperHit) {
      // Get energy deposition and layer number for this hit  Fill ntuple for each layer's energy deposit
      analysisManager->FillNtupleDColumn(i+2, absoperHit->GetEdep() / CLHEP::GeV);     // Energy deposit in this layer
      analysisManager->FillNtupleDColumn(i+17, absoperHit->GetTrackLength() / CLHEP::m);       // Layer number
      // index == absoHC->entries()-1 return the total value
    }
  }
  // G4cout << "fInteractionDepth = " << fInteractionDepth << G4endl;
  analysisManager->FillNtupleDColumn(32, fInteractionDepth); // First Interaction Depth
  analysisManager->FillNtupleIColumn(33, fInteractionLayer); // First Interaction Layer
  analysisManager->FillNtupleIColumn(34, fSecondaries);      // First Interaction No Secondaries
  analysisManager->FillNtupleIColumn(35, fInteractionType);  // First Interaction Type

  analysisManager->FillNtupleDColumn(36, fHadrInteractionDepth); // First Hadronic Interaction Depth
  analysisManager->FillNtupleIColumn(37, fHadrInteractionLayer); // First Hadronic Interaction Layer
  analysisManager->FillNtupleIColumn(38, fHadrSecondaries);      // First Hadronic Interaction No Secondaries
  analysisManager->FillNtupleIColumn(39, fHadronicTag);          // First Hadronic Type

  analysisManager->AddNtupleRow();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4c
