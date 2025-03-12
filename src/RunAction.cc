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
/// \file B4/B4c/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{ // Default filename
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  // Create directories
  // analysisManager->SetHistoDirectoryName("histograms");
  // analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("B4");
  // Note: merging ntuples is available only with Root output
  analysisManager->CreateNtuple("B4", "Edep and TrackL");
  analysisManager->CreateNtupleIColumn("Particle");
  analysisManager->CreateNtupleDColumn("Energy"); // (MeV)
  analysisManager->CreateNtupleDColumn("L0_E");
  analysisManager->CreateNtupleDColumn("L1_E");
  analysisManager->CreateNtupleDColumn("L2_E");
  analysisManager->CreateNtupleDColumn("L3_E");
  analysisManager->CreateNtupleDColumn("L4_E");
  analysisManager->CreateNtupleDColumn("L5_E");
  analysisManager->CreateNtupleDColumn("L6_E");
  analysisManager->CreateNtupleDColumn("L7_E");
  analysisManager->CreateNtupleDColumn("L8_E");
  analysisManager->CreateNtupleDColumn("L9_E");
  analysisManager->CreateNtupleDColumn("L10_E");
  analysisManager->CreateNtupleDColumn("L11_E");
  analysisManager->CreateNtupleDColumn("L12_E");
  analysisManager->CreateNtupleDColumn("L13_E");
  analysisManager->CreateNtupleDColumn("Total_E");
  analysisManager->CreateNtupleDColumn("L0_L"); 
  analysisManager->CreateNtupleDColumn("L1_L");
  analysisManager->CreateNtupleDColumn("L2_L");
  analysisManager->CreateNtupleDColumn("L3_L");
  analysisManager->CreateNtupleDColumn("L4_L");
  analysisManager->CreateNtupleDColumn("L5_L");
  analysisManager->CreateNtupleDColumn("L6_L");
  analysisManager->CreateNtupleDColumn("L7_L");
  analysisManager->CreateNtupleDColumn("L8_L");
  analysisManager->CreateNtupleDColumn("L9_L");
  analysisManager->CreateNtupleDColumn("L10_L");
  analysisManager->CreateNtupleDColumn("L11_L");
  analysisManager->CreateNtupleDColumn("L12_L");
  analysisManager->CreateNtupleDColumn("L13_L");
  analysisManager->CreateNtupleDColumn("Total_L");
  analysisManager->CreateNtupleDColumn("First_Depth"); // Depth ()
  analysisManager->CreateNtupleIColumn("First_Layer"); // Layer (0,1,2...)
  analysisManager->CreateNtupleIColumn("First_Second"); // Number of Secondaries
  analysisManager->CreateNtupleIColumn("First_Type"); // Type (0-EM, 1-HD, 2-Others)
  analysisManager->CreateNtupleDColumn("First_Had_Depth"); // Depth ()
  analysisManager->CreateNtupleIColumn("First_Had_Layer"); // Layer (0,1,2...)
  analysisManager->CreateNtupleIColumn("First_Had_Second");  // Number of Secondaries
  analysisManager->CreateNtupleIColumn("First_Had_Type");    // -1 initial; 0 exclude Inelasitc ; 1 -Inelastic 

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
