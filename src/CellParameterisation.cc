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
/// \file B4/src/CellParameterisation.cc
/// \brief Implementation of the B4::CellParameterisation class

#include "CellParameterisation.hh"

#include "Constants.hh"

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CellParameterisation::CellParameterisation()
{
  for (auto copyNo = 0; copyNo < kNofEmCells; copyNo++) 
  {
    auto iz   = copyNo / (kNofEmRows * kNofEmColumns); 
    auto iy = (copyNo % (kNofEmRows * kNofEmColumns)) / kNofEmRows;
    auto ix = (copyNo % (kNofEmRows * kNofEmColumns)) % kNofEmRows;

    if (iz % 2 == 1)  // 如果 iz 是奇数
    {
      fXCell[copyNo] = 0 ;
      fYCell[copyNo] = (iy - 11) *  2.5 * cm - 1.25 * cm;
      fZCell[copyNo] = (iz - 6)  *  2.55 * cm - 1.30 * cm;
    }
    else
    {
      fXCell[copyNo] = (iy - 11) * 2.5 * cm - 1.25 * cm;
      fYCell[copyNo] = 0 ;
      fZCell[copyNo] = (iz - 6)  * 2.55 * cm - 1.30 * cm;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void CellParameterisation::ComputeTransformation(const G4int copyNo,
                                                 G4VPhysicalVolume* physVol) const
{
  
  G4int iz = copyNo / (kNofEmRows * kNofEmColumns); // 计算层编号
  G4ThreeVector position = G4ThreeVector(fXCell[copyNo], fYCell[copyNo], fZCell[copyNo]);

  // 创建旋转矩阵
  G4RotationMatrix* rotation = nullptr;

  if (iz % 2 == 1) { // 如果 iz 是奇数
      rotation = new G4RotationMatrix();  
      rotation->rotateZ(90 * deg); // 绕 Z 轴旋转 90°
  }

  // 设置变换
  physVol->SetTranslation(position);
  physVol->SetRotation(rotation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B4
