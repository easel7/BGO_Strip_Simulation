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
/// \file B4/B4/include/EventAction.hh
/// \brief Definition of the B4::EventAction class

#ifndef B4EventAction_h
#define B4EventAction_h 1

#include "G4UserEventAction.hh"

#include "Constants.hh"

#include "CalorHit.hh"

#include "globals.hh"

#include <vector>

class G4Event;

namespace B4
{

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in the hits collections.

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void SetFirstInteraction(G4ThreeVector pos) { fFirstInteraction = pos; }
    G4ThreeVector GetFirstInteraction() const { return fFirstInteraction; }

    void SetInteractionType(G4int Type) { fInteractionType = Type; }
    G4int GetInteractionType() const { return fInteractionType; }

    void SetSecondaries(G4int Secondaries) { fSecondaries = Secondaries;}
    G4int GetSecondaries() const { return fSecondaries; }

    void SetInteractionLayer(G4int Layer) { fInteractionLayer = Layer; }
    G4int GetInteractionLayer() const { return fInteractionLayer; }

    void SetInteractionDepth(G4double InteractionDepth) { fInteractionDepth = InteractionDepth; }
    G4double GetInteractionDepth() const { return fInteractionDepth; }

    void SetHadrSecondaries(G4int Secondaries) { fHadrSecondaries = Secondaries;}
    G4int GetHadrSecondaries() const { return fHadrSecondaries; }

    void SetHadrInteractionLayer(G4int Layer) { fHadrInteractionLayer = Layer; }
    G4int GetHadrInteractionLayer() const { return fHadrInteractionLayer; }

    void SetHadrInteractionDepth(G4double InteractionDepth) { fHadrInteractionDepth = InteractionDepth; }
    G4double GetHadrInteractionDepth() const { return fHadrInteractionDepth; }

    void SetHadrTag(G4int Tag) { fHadronicTag = Tag; }
    G4int GetHadrTag() const { return fHadronicTag; }

    std::vector<G4double>& GetEmCalEdep() { return fCalEdep; }
    std::vector<G4double>& GetEmCalLeng() { return fCalLeng; }
    std::vector<G4int>&    GetLayerHits() { return fLayHits; }
    std::vector<G4double>& GetLayerEdep() { return fLayEdep; }
    std::vector<G4double>& GetLayerLeng() { return fLayLeng; }

  private:
    // methods
    CalorHitsCollection* GetHitsCollection(G4int hcID, const G4Event* event) const;
    // data members
    G4int fAbsHCID = -1;
    G4int fGapHCID = -1;

    G4ThreeVector fFirstInteraction = G4ThreeVector(-1e9, -1e9, -1e9);  // 记录第一个电磁相互作用点
    G4double fInteractionDepth = -1;
    G4int fInteractionLayer = -1;
    G4int fSecondaries = -1;
    G4int fInteractionType = -1;
    
    G4double fHadrInteractionDepth = -1;
    G4int fHadrInteractionLayer = -1;
    G4int fHadrSecondaries = -1;

    G4int fHadronicTag=-1; 

    std::vector<G4double> fCalEdep{ std::vector<G4double>(kNofEmCells, 0.) };
    std::vector<G4double> fCalLeng{ std::vector<G4double>(kNofEmCells, 0.) };
    
    std::vector<G4int>    fLayHits{ std::vector<G4int>   (kNofEmLayers, 0) };
    std::vector<G4double> fLayEdep{ std::vector<G4double>(kNofEmLayers, 0.) };
    std::vector<G4double> fLayLeng{ std::vector<G4double>(kNofEmLayers, 0.) };
};

}  // namespace B4

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
