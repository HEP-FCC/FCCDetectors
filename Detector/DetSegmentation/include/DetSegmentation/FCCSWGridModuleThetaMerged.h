#ifndef DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H
#define DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H

// FCCSW
#include "DetSegmentation/GridTheta.h"
#include "DD4hep/VolumeManager.h"

/** FCCSWGridModuleThetaMerged Detector/DetSegmentation/DetSegmentation/FCCSWGridModuleThetaMerged.h FCCSWGridModuleThetaMerged.h
 *
 *  Segmentation in theta and module.
 *  Based on GridTheta, merges modules and theta cells based on layer number
 *
 */

namespace dd4hep {
namespace DDSegmentation {
class FCCSWGridModuleThetaMerged : public GridTheta {
public:
  /// default constructor using an arbitrary type
  FCCSWGridModuleThetaMerged(const std::string& aCellEncoding);
  /// Default constructor used by derived classes passing an existing decoder
  FCCSWGridModuleThetaMerged(const BitFieldCoder* decoder);

  /// destructor
  virtual ~FCCSWGridModuleThetaMerged() = default;

  /**  Determine the local position based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Position (relative to R, phi of Geant4 volume it belongs to, scaled for R=1).
   */
  virtual Vector3D position(const CellID& aCellID) const;
  /**  Determine the cell ID based on the position.
   *   @param[in] aLocalPosition (not used).
   *   @param[in] aGlobalPosition
   *   @param[in] aVolumeId ID of the Geant4 volume
   *   return Cell ID.
   */
  virtual CellID cellID(const Vector3D& aLocalPosition, const Vector3D& aGlobalPosition,
                        const VolumeID& aVolumeID) const;
  /**  Determine the azimuthal angle (relative to the G4 volume) based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Phi.
   */
  double phi(const CellID& aCellID) const;
  /**  Determine the polar angle (relative to the G4 volume) based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Theta.
   */
  double theta(const CellID& aCellID) const;
  /**  Determine the radius based on the cell ID.
   *   @param[in] aCellId ID of a cell.
   *   return Radius.
   */
  // double radius(const CellID& aCellID) const;
  /**  Get the number of merged cells in theta for given layer
   *   @param[in] layer
   *   return The number of merged cells in theta
   */
  inline int mergedThetaCells(const unsigned int layer) const { return m_mergedCellsTheta[layer]; }
  /**  Get the number of merged modules (inclined in phi)
   *   @param[in] layer
   *   return The number of merged cells in theta
   */
  inline int mergedModules(const unsigned int layer) const { return m_mergedModules[layer]; }
  /**  Get the number of modules
   *   return The number of modules (as it was set by the user in the xml file..)
   */
  inline int nModules() const { return m_nModules; }
  /**  Set the number of modules
   *   return The number of modules (as it was set by the user in the xml file..)
   */
  inline int setNModules(const int nModules) { m_nModules = nModules; }
  /**  Get the field name used for the layer
   *   return The field name for the layer.
   */
  inline const std::string& fieldNameLayer() const { return m_layerID; }
  /**  Get the field name used for the module number
   *   return The field name for the module.
   */
  inline const std::string& fieldNameModule() const { return m_moduleID; }

protected:
  /// the field name used for layer
  std::string m_layerID;
  /// the field name used for the read-out module (can differ from module due to merging)
  std::string m_moduleID;
  /// vector of number of cells to be merged along theta for each layer
  std::vector<int> m_mergedCellsTheta;
  /// vector of number of modules to be merged for each layer
  std::vector<int> m_mergedModules;

  /// to be seen if we need these
  /// and how to retrieve them at initialization step from the geometry
  /// number of modules (or, equivalently, the deltaPhi between adjacent modules)
  int m_nModules;
  /// innermost radius of electrodes (not needed)
  /// std::vector<double> m_rLayers;
  /// inclination in phi of electrodes (not needed)
  //// double m_alpha;
};
}
}
#endif /* DETSEGMENTATION_FCCSWGRIDMODULETHETAMERGED_H */
