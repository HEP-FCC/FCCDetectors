# SimpleSensitiveLayeredCylinder.cpp Description 
This file is describing the **SimpleSensitiveLayeredCylinder.cpp**
It's a factory for a shape from multiple cylinders, These cylinders can be configured to be sensitive, and by default, they are non-sensitive if the 'sensitive' keyword is not explicitly defined.

## Expected xml structure:

```xml
<detector type="SimpleSensitiveLayeredCylinder_o1_v00" ...>
  <dimensions rmin="..." rmax="..." dz="..." z_offset="..."> 
  <sensitive type="SimpleTrackerSD"/>
  <layer rmin="..." rmax="..." dz="..." z_offset="..." material="...">

  ...
  <layer rmin="..." rmax="..." dz="..." z_offset="..." material="..." sensitive="true">
</detector>

```

 - **rmin**: This attribute represents the inner radius of the layer. It defines the minimum radial distance from the center of the detector to the inner edge of the layer.

 - **rmax**: This attribute represents the outer radius of the layer. It defines the maximum radial distance from the center of the detector to the outer edge of the layer.

 - **dz**: This attribute stands for "delta-z" and represents the half-length of the layer along the z-axis. It defines how far the layer extends along the z-direction from its center.

 - **z_offset**: This attribute specifies the offset of the layer along the z-axis with respect to the center of the detector. It defines the vertical position of the layer within the detector's volume.(Take attention that this z_offset is summed to the dimension.z_offset"envelop offset")
 
 - **material**: This attribute indicates the material of the layer. It refers to the substance from which the layer is made. The material is typically defined in the DD4hep configuration and must match a material defined in the configuration.

 - **sensitive**: This attribute is a boolean value that indicates whether the layer is a sensitive layer or not. In DD4hep, a sensitive layer is one that interacts with particles and records hits for simulation purposes. If set to "true", it means that this layer is sensitive and will contribute to particle tracking and interaction simulations.
