#include "DetSegmentation/FCCSWGridModuleThetaMerged.h"

#include <iostream>

namespace dd4hep {
namespace DDSegmentation {

/// default constructor using an encoding string
FCCSWGridModuleThetaMerged::FCCSWGridModuleThetaMerged(const std::string& cellEncoding) : GridTheta(cellEncoding) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for readout module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged modules per layer", m_mergedModules, std::vector<int>());
  registerParameter("nModules", "Number of modules", m_nModules, 1545);
  //m_nModules = 1545;
  //m_alpha = 50./180.*M_PI;
  //m_rLayers = {217.28, 218.78, 222.28, 225.78, 229.28, 232.78, 236.28, 239.78, 243.28, 246.78, 250.28, 253.78, 257.33};
}

FCCSWGridModuleThetaMerged::FCCSWGridModuleThetaMerged(const BitFieldCoder* decoder) : GridTheta(decoder) {
  // define type and description
  _type = "FCCSWGridModuleThetaMerged";
  _description = "Module-theta segmentation with per-layer merging along theta and/or module";

  // register all necessary parameters (additional to those registered in GridTheta)
  registerIdentifier("identifier_layer", "Cell ID identifier for layer", m_layerID, "layer");
  registerIdentifier("identifier_module", "Cell ID identifier for module", m_moduleID, "module");
  registerParameter("mergedCells_Theta", "Numbers of merged cells in theta per layer", m_mergedCellsTheta, std::vector<int>());
  registerParameter("mergedModules", "Numbers of merged cells in phi per layer", m_mergedModules, std::vector<int>());
  registerParameter("nModules", "Number of modules", m_nModules, 1545);

  //m_nModules = 1545;
  //m_alpha = 50./180.*M_PI;
  //m_rLayers = {217.28, 218.78, 222.28, 225.78, 229.28, 232.78, 236.28, 239.78, 243.28, 246.78, 250.28, 253.78, 257.33};
}

/// determine the local position based on the cell ID
Vector3D FCCSWGridModuleThetaMerged::position(const CellID& cID) const {

  //debug
  std::cout << "cellID: " << cID << std::endl;

  // cannot return for R=0 otherwise will lose phi info, return for R=1
  return positionFromRThetaPhi(1.0, theta(cID), phi(cID));

  // version 1: determine radius/theta/phi from cellID
  // return positionFromRThetaPhi(radius(cID), theta(cID), phi(cID));

  // version 2: use detector position 
  // read in outGlobal the x-y-z coordinates of the 1st of the merged group of cells
  /*
  VolumeID vID = cID;
  _decoder->set(vID, m_thetaID, 0);
  auto detelement = m_volman.lookupDetElement(vID);
  const auto& transformMatrix = detelement.nominal().worldTransformation();
  double outGlobal[3];
  double inLocal[] = {0, 0, 0};
  transformMatrix.LocalToMaster(inLocal, outGlobal);
  // now get R-phi
  Vector3D position(outGlobal[0],outGlobal[1],outGlobal[2]);
  double _radius = std::sqrt(outGlobal[0]*outGlobal[0] + outGlobal[1]*outGlobal[1]);
  double _phi = std::atan2(outGlobal[1], outGlobal[0]);
  // now add shifts in phi due to merging of modules and theta cells
  int layer = _decoder->get(cID, m_layerID);
  if (m_mergedModules[layer]>1) {
    _phi += (m_mergedModules[layer]-1) * M_PI / m_nModules ;
    // should we check (and adjust) that phi is in -pi..pi? or 0..2pi?
  }
  // get theta from grid
  double _theta = theta(cID);
  return positionFromRThetaPhi(_radius, _theta, _phi);
  */
}

/// determine the cell ID based on the global position
CellID FCCSWGridModuleThetaMerged::cellID(const Vector3D& /* localPosition */, const Vector3D& globalPosition,
                          const VolumeID& vID) const {
  CellID cID = vID;

  // retrieve layer (since merging depends on layer)
  int layer = _decoder->get(vID, m_layerID);

  // retrieve theta
  double lTheta = thetaFromXYZ(globalPosition);

  // calculate theta bin with original segmentation
  int thetaBin = positionToBin(lTheta, m_gridSizeTheta, m_offsetTheta);

  // adjust theta bin if cells are merged along theta in this layer
  if (m_mergedCellsTheta[layer]>1) {
    if ((thetaBin % m_mergedCellsTheta[layer])>0)
      thetaBin -= (thetaBin % m_mergedCellsTheta[layer]);
  }

  // set theta field of cellID
  _decoder->set(cID, m_thetaID, thetaBin);

  // retrieve module number
  int module = _decoder->get(vID, m_moduleID);

  // adjust module number if modules are merged in this layer
  if (m_mergedModules[layer]>1) {
    if ((module % m_mergedModules[layer])>0)
      module -= (module % m_mergedModules[layer]);
  }

  // set module field of cellID
  _decoder->set(cID, m_moduleID, module);

  return cID;
}

/// determine the radius based on the cell ID
/*
double FCCSWGridModuleThetaMerged::radius(const CellID& cID) const {

  int layer = _decoder->get(cID, m_layerID);

  return (m_rLayers[layer+1]+m_rLayers[layer])/2.0;
}
*/

/// determine the azimuth based on the cell ID
double FCCSWGridModuleThetaMerged::phi(const CellID& cID) const {

  // retrieve layer, radius and module
  int layer = _decoder->get(cID, m_layerID);

  /*
  // calculate phi of volume
  int module = _decoder->get(cID, m_moduleID);
  double r = radius(cID);

  // calculate initial phi
  // not sure this is fully consistent with Geant4 positioning!!!!
  // I need to find out which phi is module=0
  double phi = 2.0*M_PI/m_nModules * (module) - m_alpha;

  // calculate phi offset due to electrode inclination
  double r0 = m_rLayers[0];
  double L = -r0*cos(m_alpha) + std::sqrt(r*r-r0*r0*std::sin(m_alpha)*std::sin(m_alpha));
  double dphi = std::acos((r*r+r0*r0-L*L)/(2*r*r0));
  phi += dphi;
  */

  // calculate phi offset due to merging 
  double phi = 0.0;
  if (m_mergedModules[layer]>1) {
    phi += (m_mergedModules[layer]-1) * M_PI / m_nModules ;
    // should we check (and adjust) that phi is in -pi..pi? or 0..2pi?
  }

  //debug
  std::cout << "layer: " << layer << std::endl;
  std::cout << "merged modules: " << m_mergedModules[layer] << std::endl;
  std::cout << "phi: " << phi << std::endl;

  return phi;
}

double FCCSWGridModuleThetaMerged::theta(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // retrieve theta bin from cellID and determine theta position
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  // debug
  std::cout << "theta bin: " << thetaValue << std::endl;
  double _theta = binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);
  std::cout << "gridSizeTheta, offsetTheta: " << m_gridSizeTheta << " , " << m_offsetTheta << std::endl;
  // adjust return value if cells are merged along theta in this layer
  // shift by (N-1)*half theta grid size
  if (m_mergedCellsTheta[layer]>1) {
    _theta += (m_mergedCellsTheta[layer]-1) * m_gridSizeTheta / 2.0 ;
  }

  std::cout << "theta: " << _theta << std::endl;
  return _theta;
}

}
}
