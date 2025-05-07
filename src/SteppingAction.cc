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
/// \file B4/B4a/src/SteppingAction.cc
/// \brief Implementation of the B4a::SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcessType.hh"
#include "G4Step.hh"
#include "globals.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detConstruction, EventAction* eventAction)
  : fDetConstruction(detConstruction), fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Collect energy and track length step by step
  // G4cout<< "SteppingAction::UserSteppingAction" << G4endl;  
  if (fDetConstruction) {
      fDetConstruction->GetCalorThickness(); // G4cout << "Calorimeter Thickness: " << fDetConstruction->GetCalorThickness() << " mm" << G4endl;
      fDetConstruction->GetLayerThickness(); // G4cout << "Layer Thickness: " << fDetConstruction->GetLayerThickness() << " mm" << G4endl;
  }
  G4Track* track = step->GetTrack();
  G4int parentID = track->GetParentID(); // 0 表示Primary
  // Only Primary
  if (parentID != 0) return;
  if (step->GetSecondaryInCurrentStep()->size() == 0) return;
  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();         // G4cout << "Volume Name: " << volume->GetName() << G4endl;
  if (!process || !volume) return;
  if (volume != fDetConstruction->GetAbsorberPV()) return;
  G4ThreeVector position = step->GetPostStepPoint()->GetPosition(); // G4cout << "Position X: " << position.x() << G4endl; G4cout << "Position Y: " << position.y() << G4endl;  G4cout << "Position Z: " << position.z() << G4endl;
  G4double interactionDepth = position.z() + fDetConstruction->GetCalorThickness() / 2;  // Depth in detector G4cout << "Depth Z: " << position.z() + fDetConstruction->GetCalorThickness() / 2  << G4endl;
  G4int interactionLayer = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() / 22;  // Layer number, while GetCopyNumber(depth=0) find its copynumber, 
  G4int nSecondaries = 0;
  const auto* secondaries = step->GetSecondaryInCurrentStep();
  for (const auto& sec : *secondaries) {
    if (sec->GetDefinition()->GetPDGCharge() != 0.0) {
      ++nSecondaries;
    }
  }
  G4ProcessType processType = process->GetProcessType();
  G4int processSubType = process->GetProcessSubType();

  
  // G4cout << "First Interaction Position: " << fEventAction->GetInteractionDepth() << G4endl; // G4cout << "First Interaction Layer: " << fEventAction->GetInteractionLayer() << G4endl;
 
   if (fEventAction) 
  {
    // 记录第一次 interaction（EM / Had / Other）
    if (fEventAction->GetInteractionType() == -1) {
      G4int typeCode = (processType == fElectromagnetic) ? 0 :
                       (processType == fHadronic)        ? 1 : 2;
      fEventAction->SetInteractionType(typeCode);
      fEventAction->SetSecondaries(nSecondaries);
      fEventAction->SetInteractionDepth(interactionDepth);
      fEventAction->SetInteractionLayer(interactionLayer);
    }
  
    // 记录第一次强相互作用（无论是否为第一次 overall interaction）
    if (processType == fHadronic && fEventAction->GetHadrTag() == -1) {
      fEventAction->SetHadrInteractionDepth(interactionDepth);
      fEventAction->SetHadrInteractionLayer(interactionLayer);
      fEventAction->SetHadrSecondaries(nSecondaries);
      if (processSubType == fHadronInelastic)          {         fEventAction->SetHadrTag(1);        } // 非弹性强相互作用
      else if (processSubType == fHadronElastic)       {         fEventAction->SetHadrTag(2);        } // 弹性强相互作用
      else                                             {         fEventAction->SetHadrTag(0);        } // 其他强相互作用
    }
  
    // 记录第一次非弹性强相互作用
    if (processType == fHadronic && processSubType == fHadronInelastic && fEventAction->GetFirstIneTag() == -1) 
    {
      fEventAction->SetFirstIneDepth(interactionDepth);
      fEventAction->SetFirstIneLayer(interactionLayer);
      fEventAction->SetFirstIneSecondaries(nSecondaries);
      fEventAction->SetFirstIneTag(1);
    }
  }
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4a
