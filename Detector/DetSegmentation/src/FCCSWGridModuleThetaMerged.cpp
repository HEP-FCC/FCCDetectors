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
}

/// determine the local position based on the cell ID
Vector3D FCCSWGridModuleThetaMerged::position(const CellID& cID) const {

  // debug
  // std::cout << "cellID: " << cID << std::endl;

  // cannot return for R=0 otherwise will lose phi info, return for R=1
  return positionFromRThetaPhi(1.0, theta(cID), phi(cID));
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

/// determine the azimuth based on the cell ID
/// the value returned is the relative shift in phi
/// with respect to the first module in the group of
/// merged ones - which will be then added on top of
/// the phi of the volume containing the first cell
/// by the positioning tool
double FCCSWGridModuleThetaMerged::phi(const CellID& cID) const {

  // retrieve layer, radius and module
  int layer = _decoder->get(cID, m_layerID);

  // calculate phi offset due to merging 
  double phi = 0.0;
  if (m_mergedModules[layer]>1) {
    phi += (m_mergedModules[layer]-1) * M_PI / m_nModules ;
  }

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "merged modules: " << m_mergedModules[layer] << std::endl;
  // std::cout << "phi: " << phi << std::endl;

  return phi;
}

/// determine the polar angle based on the cell ID and the
/// number of merged theta cells
double FCCSWGridModuleThetaMerged::theta(const CellID& cID) const {

  // retrieve layer
  int layer = _decoder->get(cID, m_layerID);

  // retrieve theta bin from cellID and determine theta position
  CellID thetaValue = _decoder->get(cID, m_thetaID);
  double _theta = binToPosition(thetaValue, m_gridSizeTheta, m_offsetTheta);

  // adjust return value if cells are merged along theta in this layer
  // shift by (N-1)*half theta grid size
  if (m_mergedCellsTheta[layer]>1) {
    _theta += (m_mergedCellsTheta[layer]-1) * m_gridSizeTheta / 2.0 ;
  }

  // debug
  // std::cout << "layer: " << layer << std::endl;
  // std::cout << "theta bin: " << thetaValue << std::endl;
  // std::cout << "gridSizeTheta, offsetTheta: " << m_gridSizeTheta << " , " << m_offsetTheta << std::endl;
  // std::cout << "merged cells: " << m_mergedCellsTheta[layer] << std::endl;
  // std::cout << "theta: " << _theta << std::endl;

  return _theta;
}

}
}
