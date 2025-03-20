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
/// \file B4/B4/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(B4::EventAction* eventAction) : fEventAction(eventAction)
{ // Default filename
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("B4");
  // Creating 2D histograms
  analysisManager->CreateH2("Chamber1 XY", "Drift Chamber 1 X vs Z",  50, -1000., 1000, 50, -300., 300.);
  analysisManager->CreateH2("Chamber2 XY", "Drift Chamber 2 Y vs Z",  50, -1500., 1500, 50, -300., 300.);
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  analysisManager->CreateNtupleIColumn("Particle");
  analysisManager->CreateNtupleDColumn("Energy"); // (MeV)
  analysisManager->CreateNtupleDColumn("First_Depth"); // Depth ()
  analysisManager->CreateNtupleIColumn("First_Layer"); // Layer (0,1,2...)
  analysisManager->CreateNtupleIColumn("First_Second"); // Number of Secondaries
  analysisManager->CreateNtupleIColumn("First_Type"); // Type (0-EM, 1-HD, 2-Others)
  analysisManager->CreateNtupleDColumn("First_Had_Depth"); // Depth ()
  analysisManager->CreateNtupleIColumn("First_Had_Layer"); // Layer (0,1,2...)
  analysisManager->CreateNtupleIColumn("First_Had_Second");  // Number of Secondaries
  analysisManager->CreateNtupleIColumn("First_Had_Type");    // -1 initial; 0 exclude Inelasitc ; 1 -Inelastic 

  analysisManager->CreateNtupleIColumn("LayerHitsVector"  , fEventAction->GetLayerHits()); 
  analysisManager->CreateNtupleDColumn("LayerEnergyVector", fEventAction->GetLayerEdep());
  analysisManager->CreateNtupleDColumn("LayerLengthVector", fEventAction->GetLayerLeng());

  analysisManager->CreateNtupleDColumn("BarEnergyVector", fEventAction->GetEmCalEdep());
  analysisManager->CreateNtupleDColumn("BarLengthVector", fEventAction->GetEmCalLeng());

  analysisManager->CreateNtupleIColumn("Nhits");
  analysisManager->CreateNtupleDColumn("Total_E");
  analysisManager->CreateNtupleDColumn("Total_L");
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  // inform the runManager to save random number seed
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  auto analysisManager = G4AnalysisManager::Instance();
  // Reset histograms from previous run
  // analysisManager->Reset();
  // Open an output file
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  auto analysisManager = G4AnalysisManager::Instance();
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4
