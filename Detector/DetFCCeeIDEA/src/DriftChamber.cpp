/*****************************************************************************\
* DD4hep geometry code for the central drift chamber of the IDEA detector     *
* Author: Lorenzo Capriotti                                                   *
\*****************************************************************************/
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "DD4hep/detail/DetectorInterna.h"
#include "TClass.h"
#include "TMath.h"
#include "XML/Utilities.h"
#include <DD4hep/DetFactoryHelper.h>
#include <XML/Layering.h>
#include <iostream>

using namespace std;
using namespace dd4hep;

struct wire

{

  dd4hep::Volume layer;
  string type;
  int num;
  double radius;
  double theta;
  double thetaoffset;
  double stereo;
  double halfalpha;
  double thickness;
  double halflength;
  dd4hep::Volume volume;

};

namespace {

  struct CDCHBuild : public dd4hep::xml::tools::VolumeBuilder {
    std::vector<dd4hep::DetElement> deSuperLayer, deLayer, deSWire;

    CDCHBuild( dd4hep::Detector& description, xml_elt_t e, dd4hep::SensitiveDetector sens );

    double diff_of_squares(double a, double b);
    void PlaceWires( struct wire &w, double outwrap, double halflength, int copyNunOffset, int SL, int iring, int wirenum);
    void build_layer(DetElement parent, Volume parentVol);

  };

  // ******************************************************
  // Initializing constructor
  // ******************************************************

  CDCHBuild::CDCHBuild( dd4hep::Detector& dsc, xml_elt_t e, dd4hep::SensitiveDetector sens )
      : dd4hep::xml::tools::VolumeBuilder( dsc, e, sens ) {}

  double CDCHBuild::diff_of_squares( double a, double b) {
    
    double diff = pow(a,2) - pow(b,2);
    return diff;

  }


  void CDCHBuild::PlaceWires( struct wire &w, double outwrap, double halflength, int copyNunOffset=0, int SL=999, int iring=999, int wirenum=-1) {

    dd4hep::RotationZYX rot( 0., 0., w.stereo );
    dd4hep::RotationX     rot_stereo( w.stereo );
    dd4hep::Position    pos( w.radius, 0., 0. );
    dd4hep::Translation3D transl( w.radius, 0., 0. );

    dd4hep::Transform3D T(transl*rot_stereo);

    string wirewrapname = "lvWireWrap_SL";
    wirewrapname += std::to_string(SL);
    wirewrapname += "_ring";
    wirewrapname += std::to_string(iring);
    wirewrapname += "_type";
    wirewrapname += w.type;
    wirewrapname += "_stereo";
    wirewrapname += std::to_string(w.stereo);

    cout << "wirewrapname: " << wirewrapname << endl;

    string wirename = "lvWire_SL";
    wirename += std::to_string(SL);
    wirename += "_ring";
    wirename += std::to_string(iring);
    wirename += "_type";
    wirename += w.type;
    wirename += "_stereo";
    wirename += std::to_string(w.stereo);


    dd4hep::Tube WrapTube(w.thickness, w.thickness+0.5*outwrap, halflength);
    dd4hep::Volume lvWireWrapVol(wirewrapname, WrapTube, description.material( "G4_Au") );

    dd4hep::Tube TotalWire(0.0, w.thickness+0.5*outwrap, halflength);
    dd4hep::Volume lvWireVol(wirename, TotalWire, description.material("Air"));

    lvWireVol.placeVolume( w.volume, dd4hep::Position( 0.0, 0.0, 0.0) );
    lvWireVol.placeVolume( lvWireWrapVol, dd4hep::Position( 0.0, 0.0, 0.0) );

//    registerVolume(lvWireWrapVol.name(), lvWireWrapVol);
//    registerVolume(lvWireVol.name(), lvWireVol);

    for (int n=0; n<w.num; n++) {
        dd4hep::RotationZ iRot(w.thetaoffset + w.theta*n);
        if (n%1==0) w.layer.placeVolume( lvWireVol, dd4hep::Transform3D( iRot*T ) );
    }
  }



  void CDCHBuild::build_layer(DetElement parent, Volume parentVol) {

    // ******************************************************
    // Loading parameters
    // ******************************************************

    double halfalpha = 0.5 * dd4hep::_toDouble("CDCH:alpha");
    double inner_radius = dd4hep::_toDouble( "CDCH:r0" );
    double outer_radius = dd4hep::_toDouble( "CDCH:rOut" );
    double halflength = dd4hep::_toDouble( "CDCH:zHalfLength" );
    double CarbonInnerWallThick = dd4hep::_toDouble("CDCH:CarbonInnerWallThick");
    double CopperInnerWallThick	= dd4hep::_toDouble("CDCH:CopperInnerWallThick");
    double GasInnerWallThick	= dd4hep::_toDouble("CDCH:GasInnerWallThick");
    double Carbon1OuterWallThick = dd4hep::_toDouble("CDCH:Carbon1OuterWallThick");
    double Carbon2OuterWallThick = dd4hep::_toDouble("CDCH:Carbon2OuterWallThick");
    double CopperOuterWallThick = dd4hep::_toDouble("CDCH:CopperOuterWallThick");
    double FoamOuterWallThick    = dd4hep::_toDouble("CDCH:FoamOuterWallThick");
    double GasEndcapWallThick = dd4hep::_toDouble("CDCH:GasEndcapWallThick");
    double CopperEndcapWallThick = dd4hep::_toDouble("CDCH:CopperEndcapWallThick");
    double KaptonEndcapWallThick = dd4hep::_toDouble("CDCH:KaptonEndcapWallThick");
    double CarbonEndcapWallThick = dd4hep::_toDouble("CDCH:CarbonEndcapWallThick");
    double FWireShellThickIn = dd4hep::_toDouble( "CDCH:FWireShellThickIn" );
    double FWireShellThickOut = dd4hep::_toDouble( "CDCH:FWireShellThickOut" );
    double SWireShellThickIn = dd4hep::_toDouble( "CDCH:SWireShellThickIn" );
    double SWireShellThickOut = dd4hep::_toDouble( "CDCH:SWireShellThickOut" );
    double CntFWireShellThickIn = dd4hep::_toDouble( "CDCH:CntFWireShellThickIn" );
    double CntFWireShellThickOut = dd4hep::_toDouble( "CDCH:CntFWireShellThickOut" );
    double InGWireShellThickIn = dd4hep::_toDouble( "CDCH:InGWireShellThickIn" );
    double InGWireShellThickOut = dd4hep::_toDouble( "CDCH:InGWireShellThickOut" );
    double OutGWireShellThickIn = dd4hep::_toDouble( "CDCH:InGWireShellThickIn" );
    double OutGWireShellThickOut = dd4hep::_toDouble( "CDCH:InGWireShellThickOut" );
    double secure = dd4hep::_toDouble( "CDCH:secure" );
    double capGasLayer = dd4hep::_toDouble( "CDCH:capGasLayer" );
    double extShiftFW = dd4hep::_toDouble( "CDCH:extShiftFW" );
    double cellDimension = dd4hep::_toDouble( "CDCH:cellDimension" );
    double inGuardRad = dd4hep::_toDouble( "CDCH:inGuardRad" );
    double outGuardRad = dd4hep::_toDouble( "CDCH:outGuardRad" );
    int nSDeltaWire = dd4hep::_toInt("CDCH:nSDeltaWire");
    int nSWire = dd4hep::_toInt("CDCH:nSWire");
    int nInGWire = dd4hep::_toInt("CDCH:nInGWire");
    int nOutGWire = dd4hep::_toInt("CDCH:nOutGWire");
    int nStoFWireRatio = dd4hep::_toInt("CDCH:nStoFWireRatio");
    int nVerticalFWire = dd4hep::_toInt("CDCH:nVerticalFWire");
    int nSuperLayer = dd4hep::_toInt("CDCH:nSuperLayer");
    int	nRing = dd4hep::_toInt("CDCH:nRing");
    int	nFieldWireShells = dd4hep::_toInt("CDCH:nFieldWireShells");

    double epsilon = 0.0;
    double theta_ring = 0.0;
    double theta_ring1 = 0.0;
    int nFWire = 0;
    int nFWire1 = 0;
    int num_wire = 0;
    int nHorizontalFWire = nStoFWireRatio - nVerticalFWire;
    int sign_epsilon = -1;
    double phi = 0.0;
    double scaleFactor = 0.0;
    double dropFactor = 0.0;
    double epsilonFactor = 0.0;
    double delta_radius_ring = cellDimension;
    double senseWireRing_radius_0 = 0.0;
    double iradius = 0.0;
    double idelta_radius = 0.0;

    double envelop_Inner_thickness = CarbonInnerWallThick + CopperInnerWallThick + GasInnerWallThick;
    double envelop_Outer_thickness = Carbon1OuterWallThick + Carbon2OuterWallThick + CopperOuterWallThick + FoamOuterWallThick;
    double FWireDiameter = FWireShellThickIn + FWireShellThickOut;
    double FWradii = 0.5*FWireDiameter;
    double CntFWireDiameter = CntFWireShellThickIn + CntFWireShellThickOut;
    double CntFWradii = 0.5*CntFWireDiameter;
    double SWireDiameter = SWireShellThickIn + SWireShellThickOut;
    double SWradii = 0.5*SWireDiameter;
    double inGWireDiameter = InGWireShellThickIn + InGWireShellThickOut;
    double inGWradii = 0.5*inGWireDiameter;
    double fakeLayerInIWthick = -0.0001+GasInnerWallThick;
    double inner_radius_0 = inner_radius+envelop_Inner_thickness-fakeLayerInIWthick;

    double radius_ring_0 = inner_radius + envelop_Inner_thickness + FWradii + secure + capGasLayer;
    double radius_ringOut_0  = radius_ring_0-FWradii-secure;

    double drop = 0.0;
    double radius_ring = 0.0;
    double radius_ringIn_0 = 0.0;
    double radius_ringIn = 0.0;
    double radius_ringOut = 0.0;
    double epsilonIn = 0.0;
    double epsilonOut = 0.0;
    double ringangle = 0.0;
    double cellBase = 0.0;
    double inscribedRadius = 0.0;
    double circumscribedRadius = 0.0;
    double zlength = 0.0;
    double cellStaggering = 0.0;
    double epsilonInGwRing = 0.0;
    double epsilonOutGwRing = 0.0;

    //------------------------------------------------------------------------
    // The enlarge parameter is used to see the wires in the rendering
    //------------------------------------------------------------------------

    double enlarge = 60.;

    //------------------------------------------------------------------------
    // Build the inner, outer and endcap walls first
    //------------------------------------------------------------------------

    dd4hep::Tube Endcap_Gas(inner_radius, outer_radius, 0.5*GasEndcapWallThick);
    dd4hep::Tube Endcap_Copper(inner_radius, outer_radius, 0.5*CopperEndcapWallThick);
    dd4hep::Tube Endcap_Kapton(inner_radius, outer_radius, 0.5*KaptonEndcapWallThick);
    dd4hep::Tube Endcap_Carbon(inner_radius, outer_radius, 0.5*CarbonEndcapWallThick);

    dd4hep::Volume lvEndcapWallGas = dd4hep::Volume("lvEndcapWallGasVol", Endcap_Gas, description.material( "GasHe_90Isob_10" ));
    dd4hep::Volume lvEndcapWallCopper = dd4hep::Volume("lvEndcapWallCopperVol", Endcap_Copper, description.material( "G4_Cu" ));
    dd4hep::Volume lvEndcapWallKapton = dd4hep::Volume("lvEndcapWallKaptonVol", Endcap_Kapton, description.material( "Kapton" ));
    dd4hep::Volume lvEndcapWallCarbon = dd4hep::Volume("lvEndcapWallCarbonVol", Endcap_Carbon, description.material( "CarbonFiber" ));


    dd4hep::Tube InnerWall_Carbon(inner_radius, inner_radius+CarbonInnerWallThick, halflength);
    dd4hep::Tube InnerWall_Copper(inner_radius+CarbonInnerWallThick, inner_radius+CarbonInnerWallThick+CopperInnerWallThick, halflength);
    dd4hep::Tube InnerWall_Gas(inner_radius+CarbonInnerWallThick+CopperInnerWallThick, inner_radius+envelop_Inner_thickness, halflength);

    dd4hep::Volume lvInnerWallCarbon = dd4hep::Volume("lvInnerWallCarbonVol", InnerWall_Carbon, description.material( "CarbonFiber" ));
    dd4hep::Volume lvInnerWallCopper = dd4hep::Volume("lvInnerWallCopperVol", InnerWall_Copper, description.material( "G4_Cu" ));
    dd4hep::Volume lvInnerWallGas = dd4hep::Volume("lvInnerWallGasVol", InnerWall_Gas, description.material( "GasHe_90Isob_10" ));


    dd4hep::Tube OuterWall_Copper(outer_radius-envelop_Outer_thickness, outer_radius-Carbon1OuterWallThick-Carbon2OuterWallThick-FoamOuterWallThick, halflength);
    dd4hep::Tube OuterWall_Carbon1(outer_radius-Carbon1OuterWallThick-Carbon2OuterWallThick-FoamOuterWallThick, outer_radius-Carbon2OuterWallThick-FoamOuterWallThick, halflength);
    dd4hep::Tube OuterWall_Foam(outer_radius-Carbon2OuterWallThick-FoamOuterWallThick, outer_radius-Carbon2OuterWallThick, halflength);
    dd4hep::Tube OuterWall_Carbon2(outer_radius-Carbon2OuterWallThick, outer_radius, halflength);

    dd4hep::Volume lvOuterWallCarbon1 = dd4hep::Volume("lvOuterWallCarbon1Vol", OuterWall_Carbon1, description.material( "CarbonFiber" ));
    dd4hep::Volume lvOuterWallCarbon2 = dd4hep::Volume("lvOuterWallCarbon2Vol", OuterWall_Carbon2, description.material( "CarbonFiber" ));
    dd4hep::Volume lvOuterWallCopper = dd4hep::Volume("lvOuterWallCopperVol", OuterWall_Copper, description.material( "G4_Cu" ));
    dd4hep::Volume lvOuterWallFoam = dd4hep::Volume("lvOuterWallFoamVol", OuterWall_Foam, description.material( "GasHe_90Isob_10" ));

    //------------------------------------------------------------------------
    // Now we are ready to loop over the SuperLayers and fill the gas volume!
    //------------------------------------------------------------------------


    std::vector<dd4hep::Volume> lvLayerVol;
    std::vector<dd4hep::Hyperboloid> HypeLayer;
    std::vector<dd4hep::Volume> lvFwireVol, lvSwireVol, lvGwireVol, lvWrapVol;

    string lvName1, lvName2, lvName3, HypeName1, HypeName2, HypeName3, wirecol, gascol;
    string lvFwireName, lvSwireName;

    struct wire w;
    //nSuperLayer = 1;

    for (int SL = 0; SL<nSuperLayer; ++SL) {

      num_wire = nSWire + SL*nSDeltaWire;
      phi = 2.*TMath::Pi() / num_wire;
      nFWire = nHorizontalFWire * num_wire;
      theta_ring = 2.*TMath::Pi() / nFWire;
      nFWire1 = nFWire/2;
      theta_ring1 = 2.0*theta_ring;
      scaleFactor = (1.0 + TMath::Pi() / num_wire)/(1.0 - TMath::Pi() / num_wire);
      dropFactor = (1.0/cos(halfalpha)-1.0);
      epsilonFactor = sin(halfalpha)/halflength;
      ringangle = -0.5*phi;

      if (SL%3==0) gascol = "vCDCH:Gas1";
      else if ((SL+1)%3==0) gascol = "vCDCH:Gas2";
      else if ((SL+2)%3==0) gascol = "vCDCH:Gas3";

      if (SL%3==0) wirecol = "vCDCH:Wire1";
      else if ((SL+1)%3==0) wirecol = "vCDCH:Wire2";
      else if ((SL+2)%3==0) wirecol = "vCDCH:Wire3";


      if (SL==0) {

        double stereoOut0 = atan(radius_ringOut_0 * (1.0*dropFactor*epsilonFactor));

        dd4hep::Hyperboloid HypeLayer0(inner_radius_0, 0.0, radius_ringOut_0-secure, stereoOut0, halflength);
        lvLayerVol.push_back( dd4hep::Volume("lvLayerInit", HypeLayer0, description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, "vCDCH:Pb" );

        epsilonInGwRing  = atan(inGuardRad*(1.0+dropFactor)*epsilonFactor);
        zlength = halflength;
        zlength-=sin(epsilonInGwRing)*inGWradii;
        zlength/=cos(epsilonInGwRing);

        w.layer = lvLayerVol.back();
        w.type = "G";
        w.num = nInGWire/2;
        w.radius = inGuardRad-inGWradii;
        w.theta = theta_ring1;
        w.thetaoffset = ringangle;
        w.stereo = epsilonInGwRing;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*InGWireShellThickIn*enlarge; //half the inner thickness as radius of tube
        w.halflength = zlength;

        dd4hep::Tube Gwire(0.0, w.thickness, halflength);
        lvGwireVol.push_back(dd4hep::Volume("Gwire_inner", Gwire, description.material( "G4_Al") ));
        lvGwireVol.back().setVisAttributes( description, wirecol );

        w.volume = lvGwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, -1);

        w.radius = inGuardRad+inGWradii+extShiftFW;
        w.thetaoffset = ringangle+theta_ring;
        w.stereo = -1.0*epsilonInGwRing;
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, nInGWire/2, SL, -1);


        drop = radius_ring_0*dropFactor;
        radius_ring = radius_ring_0+drop;
        epsilon = atan(radius_ring*epsilonFactor);
        radius_ringIn_0  = radius_ring_0-FWradii-2.0*secure;
        radius_ringIn    = radius_ringIn_0+drop;
        radius_ringOut_0 = radius_ring_0+FWradii;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonIn        = atan(sqrt(pow(radius_ringIn,2) - pow(radius_ringIn_0,2)) / halflength);
        epsilonOut       = atan(sqrt(pow(radius_ringOut,2)- pow(radius_ringOut_0,2))/halflength);
    

        dd4hep::Hyperboloid HypeLayer1(radius_ringIn_0, epsilonIn, radius_ringOut_0, epsilonOut, halflength);
        lvLayerVol.push_back( dd4hep::Volume("lvLayer_0", HypeLayer1, description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, "vCDCH:Plastic" );


        zlength = halflength;
	zlength-=sin(epsilon)*FWradii;
	zlength/=cos(epsilon);

        w.layer = lvLayerVol.back();
        w.type = "F";
        w.num = nFWire1;
        w.radius = radius_ringIn_0-FWradii-extShiftFW;
        w.theta = theta_ring1;
        w.thetaoffset = ringangle+cellStaggering-theta_ring;
        w.stereo = sign_epsilon*epsilon;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*FWireShellThickIn*enlarge;
        w.halflength = zlength;

        lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_init");

        dd4hep::Tube Fwire(0.0, w.thickness, halflength);
        lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material( "G4_Al") ));
        lvFwireVol.back().setVisAttributes( description, wirecol );

        w.volume = lvFwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, -1);

	radius_ring_0+=FWradii;

      } else {

        delta_radius_ring = 2.*TMath::Pi()*radius_ringOut_0/(num_wire-TMath::Pi());

      }

      //------------------------------------------------------------------------
      // Starting the layer ("ring") loop. nRing=8
      //------------------------------------------------------------------------

      for (int iring=0; iring< nRing; iring++ ){

        //------------------------------------------------------------------------
        // First define some useful names for the future
        //------------------------------------------------------------------------


        HypeName1 = dd4hep::_toString(SL, "HypeLayer1_%d") + dd4hep::_toString(iring, "_%d");
        HypeName2 = dd4hep::_toString(SL, "HypeLayer2_%d") + dd4hep::_toString(iring, "_%d");
        HypeName3 = dd4hep::_toString(SL, "HypeLayer3_%d") + dd4hep::_toString(iring, "_%d");
        lvName1 = dd4hep::_toString(SL, "lvLayer1_%d") + dd4hep::_toString(iring, "_%d");
        lvName2 = dd4hep::_toString(SL, "lvLayer2_%d") + dd4hep::_toString(iring, "_%d");
        lvName3 = dd4hep::_toString(SL, "lvLayer3_%d") + dd4hep::_toString(iring, "_%d");

        //------------------------------------------------------------------------
        // Next, fill the geometry parameters of the layer. Each layer lies
        // on top of the following one, so new ringIn = old ringOut 
        //------------------------------------------------------------------------


        inscribedRadius        = 0.5*delta_radius_ring;
        circumscribedRadius    = inscribedRadius*sqrt(2.0);
        senseWireRing_radius_0 = radius_ring_0+inscribedRadius;
        sign_epsilon *=-1;

        radius_ringIn_0 = radius_ringOut_0;
        radius_ringIn = radius_ringOut;
        epsilonIn        = epsilonOut;

        radius_ringOut_0 = radius_ringIn_0+FWireDiameter+2.0*secure;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonOut       = atan(sqrt(diff_of_squares(radius_ringOut, radius_ringOut_0))/halflength);

        zlength = halflength;

        //------------------------------------------------------------------------
        // Build the hyperboloid shape and the volume of the layer. This is the
        // base layer of the cell.
        //------------------------------------------------------------------------

        HypeLayer.push_back( dd4hep::Hyperboloid(radius_ringIn_0,epsilonIn, radius_ringOut_0, epsilonOut, zlength ) );
        lvLayerVol.push_back( dd4hep::Volume(lvName1, HypeLayer.back(), description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, gascol );

        //------------------------------------------------------------------------
        // Reduce zlength to avoid volume extrusions and check the staggering
        //------------------------------------------------------------------------

        zlength-=sin(epsilon)*FWradii;
        zlength/=cos(epsilon);

        if (iring%2==1) cellStaggering=theta_ring;
        else cellStaggering=0.0;

        //------------------------------------------------------------------------
        // Fill the field wire struct with all the relevant information to be 
        // passed to the PlaceWires function. This is the base of the cell.
        //------------------------------------------------------------------------

        w.layer = lvLayerVol.back();
        w.type = "F";
        w.num = nFWire1;
       	w.radius = radius_ringIn_0+FWradii+extShiftFW;
       	w.theta = theta_ring1;
       	w.thetaoffset = ringangle+cellStaggering;
       	w.stereo = sign_epsilon*epsilon;
       	w.halfalpha = halfalpha;
       	w.thickness = 0.5*FWireShellThickIn*enlarge;
       	w.halflength = zlength;

        //------------------------------------------------------------------------
        // Define the field wire name and build the field wire volume
	//------------------------------------------------------------------------

        lvFwireName = dd4hep::_toString(SL, "lvFwire_%d") + dd4hep::_toString(iring, "_%d");

        dd4hep::Tube Fwire(0.0, w.thickness, halflength);
        lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material( "G4_Al") ));
	lvFwireVol.back().setVisAttributes( description, wirecol );

        //------------------------------------------------------------------------
        // Add the field wire volume to the struct and call the magic function
        //------------------------------------------------------------------------

	w.volume = lvFwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, iring);

        //------------------------------------------------------------------------
        // Next, fill the geometry parameters of the central layer.
        //------------------------------------------------------------------------

        iradius          = radius_ring_0;
        radius_ring_0    += delta_radius_ring;
        drop             = radius_ring_0*dropFactor;

        radius_ringIn_0  = radius_ringOut_0;
        radius_ringIn    = radius_ringOut;
        epsilonIn        = epsilonOut;
        radius_ringOut_0 = radius_ring_0-FWireDiameter-2.0*secure;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonOut       = atan(sqrt(diff_of_squares(radius_ringOut, radius_ringOut_0))/halflength);
        zlength = halflength;
     
        //------------------------------------------------------------------------
        // Build the hyperboloid shape and the volume of the layer. This is the
        // central layer of the cell.
        //------------------------------------------------------------------------

        HypeLayer.push_back( dd4hep::Hyperboloid(radius_ringIn_0,epsilonIn, radius_ringOut_0, epsilonOut,zlength));
        lvLayerVol.push_back( dd4hep::Volume(lvName2, HypeLayer.back(), description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, gascol );

        //------------------------------------------------------------------------
        // Reduce zlength to avoid volume extrusions
        //------------------------------------------------------------------------

	zlength-=sin(epsilon)*CntFWradii;
        zlength/=cos(epsilon);

        //------------------------------------------------------------------------
        // Fill the sense wire struct with all the relevant information to be
        // passed to the PlaceWires function. 
        //------------------------------------------------------------------------

        w.layer = lvLayerVol.back();
        w.type = "S";
        w.num = num_wire;
        w.radius = senseWireRing_radius_0;
        w.theta = phi;
        w.thetaoffset = cellStaggering;
        w.stereo = sign_epsilon*epsilon;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*SWireShellThickIn*enlarge;
        w.halflength = zlength;

        //------------------------------------------------------------------------
        // Define the sense wire name and build the sense wire volume
        //------------------------------------------------------------------------

        lvSwireName = dd4hep::_toString(SL, "lvSwire_%d") + dd4hep::_toString(iring, "_%d");

        dd4hep::Tube Swire(0.0, w.thickness, halflength);
        lvSwireVol.push_back(dd4hep::Volume(lvSwireName, Swire, description.material( "G4_W") ));
        lvSwireVol.back().setVisAttributes( description, wirecol );

        //------------------------------------------------------------------------
        // Add the sense wire volume to the struct and call the magic function
        //------------------------------------------------------------------------

        w.volume = lvSwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, iring);

        //------------------------------------------------------------------------
        // Tune the radius and epsilon of the central field wires
        //------------------------------------------------------------------------

        idelta_radius = delta_radius_ring * 0.5;
        iradius+=idelta_radius;
        epsilon = atan(iradius*(1.0+dropFactor)*epsilonFactor);

        //------------------------------------------------------------------------
        // Fill the central field wire struct with all the relevant information 
        // and call the magic function.
        //------------------------------------------------------------------------

        w.layer = lvLayerVol.back();
        w.type = "F";
        w.num = num_wire;
        w.radius = iradius;
        w.theta = phi;
        w.thetaoffset = ringangle+cellStaggering;
        w.stereo = sign_epsilon*epsilon;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*FWireShellThickIn*enlarge;
        w.halflength = zlength;
        w.volume = lvFwireVol.back();

        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, iring);

        //------------------------------------------------------------------------
        // Next, fill the geometry parameters of the upper layer.
        //------------------------------------------------------------------------

        radius_ringIn_0  = radius_ringOut_0;
        radius_ringIn    = radius_ringOut;
        epsilonIn        = epsilonOut;
        radius_ringOut_0 = radius_ringIn_0+FWireDiameter+2.0*secure;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonOut       = atan(sqrt(diff_of_squares(radius_ringOut, radius_ringOut_0))/halflength);
        zlength = halflength;

        //------------------------------------------------------------------------
        // Build the hyperboloid shape and the volume of the layer. This is the
        // central layer of the cell.
        //------------------------------------------------------------------------

        HypeLayer.push_back( dd4hep::Hyperboloid(radius_ringIn_0,epsilonIn, radius_ringOut_0, epsilonOut,zlength));
        lvLayerVol.push_back( dd4hep::Volume(lvName3, HypeLayer.back(), description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, gascol );

        //------------------------------------------------------------------------
        // Reduce zlength to avoid volume extrusions
        //------------------------------------------------------------------------

        zlength-=sin(epsilon)*FWradii;
        zlength/=cos(epsilon);

        //------------------------------------------------------------------------
        // Fill the field wire struct with all the relevant information and 
        // call the magic function. This is the top of the cell.
        //------------------------------------------------------------------------

        w.layer = lvLayerVol.back();
        w.type = "F";
        w.num = nFWire1;
        w.radius = radius_ringIn_0-FWradii-extShiftFW;
        w.theta = theta_ring1;
        w.thetaoffset = ringangle+cellStaggering;
        w.stereo = sign_epsilon*epsilon;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*FWireShellThickIn*enlarge;
        w.halflength = zlength;
        w.volume = lvFwireVol.back();

        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, iring);

        //------------------------------------------------------------------------
        // Scale the delta radius of the ring for next iteration
        //------------------------------------------------------------------------

        delta_radius_ring *= scaleFactor;

      }

      if ( SL== (nSuperLayer-1) ) {

        radius_ringIn_0  = radius_ringOut_0;
        radius_ringIn    = radius_ringOut;
        epsilonIn        = epsilonOut;
        radius_ringOut_0 = radius_ring_0+FWireDiameter+2.0*secure;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonOut       = atan(sqrt(diff_of_squares(radius_ringOut, radius_ringOut_0))/halflength);

        dd4hep::Hyperboloid HypeLayerOut(radius_ringIn_0, epsilonIn, radius_ringOut_0, epsilonOut, halflength);
        lvLayerVol.push_back( dd4hep::Volume("lvLayerOut", HypeLayerOut, description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, "vCDCH:Plastic" );

        zlength = halflength;
        zlength-=sin(epsilon)*FWradii;
        zlength/=cos(epsilon);

        w.layer = lvLayerVol.back();
        w.type = "F";
        w.num = nFWire1;
        w.radius = radius_ringIn_0+FWradii+extShiftFW;
        w.theta = theta_ring1;
        w.thetaoffset = ringangle+cellStaggering+theta_ring;
        w.stereo = -1.*epsilon;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*FWireShellThickIn*enlarge;
        w.halflength = zlength;

        lvFwireName = dd4hep::_toString(SL, "lvFwire_%d_out");

        dd4hep::Tube Fwire(0.0, w.thickness, halflength);
        lvFwireVol.push_back(dd4hep::Volume(lvFwireName, Fwire, description.material( "G4_Al") ));
        lvFwireVol.back().setVisAttributes( description, wirecol );

        w.volume = lvFwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, -1);
 

        //------------------------------------------------------------------------
        // Start placing the outer layer of guard wires
        //------------------------------------------------------------------------

        radius_ringIn_0  = radius_ringOut_0;
        radius_ringIn    = radius_ringOut;
        epsilonIn        = epsilonOut;
        radius_ringOut_0 = radius_ring_0+FWireDiameter+2.0*secure;
        radius_ringOut   = radius_ringOut_0+drop;
        epsilonOut	 = atan(sqrt(diff_of_squares(radius_ringOut, radius_ringOut_0))/halflength);

        dd4hep::Hyperboloid HypeLayerOutG(radius_ringIn_0, epsilonOut, outer_radius-envelop_Outer_thickness-0.0001, 0.0, halflength);
        lvLayerVol.push_back( dd4hep::Volume("lvLayerOutG", HypeLayerOutG, description.material( "GasHe_90Isob_10" ) ) );
        lvLayerVol.back().setVisAttributes( description, "vCDCH:Pb" );

        epsilonOutGwRing  = atan(outGuardRad*(1.0+dropFactor)*epsilonFactor);
        zlength = halflength;
        zlength-=sin(epsilonOutGwRing)*inGWradii;
        zlength/=cos(epsilonOutGwRing);

        w.layer = lvLayerVol.back();
        w.type = "G";
        w.num = nOutGWire/2;
        w.radius = outGuardRad-inGWradii;
        w.theta = theta_ring1;
        w.thetaoffset = ringangle;
        w.stereo = epsilonOutGwRing;
        w.halfalpha = halfalpha;
        w.thickness = 0.5*OutGWireShellThickIn*enlarge;
        w.halflength = zlength;

        dd4hep::Tube Gwire(0.0, w.thickness, halflength);
        lvGwireVol.push_back(dd4hep::Volume("Gwire_outer", Gwire, description.material( "G4_Al") ));
        lvGwireVol.back().setVisAttributes( description, wirecol );

        w.volume = lvGwireVol.back();
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, 0, SL, -1);

        w.radius = outGuardRad+inGWradii+extShiftFW;
        w.thetaoffset = ringangle+theta_ring;
        w.stereo = -1.0*epsilonOutGwRing;
        CDCHBuild::PlaceWires(w, FWireShellThickOut, halflength, nOutGWire/2, SL, -1);

      }

    }  

    dd4hep::PlacedVolume pv;
    dd4hep::DetElement   CDCHDetector( parent, "Ecal_DP", parent.id() );

    Int_t sizeLayer = lvLayerVol.size();

    for (Int_t i=0; i<sizeLayer; i++) {
        registerVolume( lvLayerVol.at(i).name(), lvLayerVol.at(i) );
        cout << "Placing Volume: " << lvLayerVol.at(i).name() << endl;
        pv = parentVol.placeVolume( volume( lvLayerVol.at(i).name() ) );
        CDCHDetector.setPlacement( pv );
    }

    double PosEndcapGas = halflength + 0.5*GasEndcapWallThick;
    double PosEndcapCopper = halflength + GasEndcapWallThick + 0.5*CopperEndcapWallThick;
    double PosEndcapKapton = halflength	+ GasEndcapWallThick + CopperEndcapWallThick + 0.5*KaptonEndcapWallThick;
    double PosEndcapCarbon = halflength	+ GasEndcapWallThick + CopperEndcapWallThick + KaptonEndcapWallThick + 0.5*CarbonEndcapWallThick;

    parentVol.placeVolume(lvInnerWallCarbon);
    parentVol.placeVolume(lvInnerWallCopper);
    parentVol.placeVolume(lvInnerWallGas);
    parentVol.placeVolume(lvOuterWallCarbon1);
    parentVol.placeVolume(lvOuterWallCarbon2);
    parentVol.placeVolume(lvOuterWallCopper);
    parentVol.placeVolume(lvOuterWallFoam);
    parentVol.placeVolume(lvEndcapWallGas, dd4hep::Position(0., 0., PosEndcapGas));
    parentVol.placeVolume(lvEndcapWallCopper, dd4hep::Position(0., 0., PosEndcapCopper));
    parentVol.placeVolume(lvEndcapWallKapton, dd4hep::Position(0., 0., PosEndcapKapton));
    parentVol.placeVolume(lvEndcapWallCarbon, dd4hep::Position(0., 0., PosEndcapCarbon));
    parentVol.placeVolume(lvEndcapWallGas, dd4hep::Position(0., 0., -PosEndcapGas));
    parentVol.placeVolume(lvEndcapWallCopper, dd4hep::Position(0., 0., -PosEndcapCopper)); 
    parentVol.placeVolume(lvEndcapWallKapton, dd4hep::Position(0., 0., -PosEndcapKapton)); 
    parentVol.placeVolume(lvEndcapWallCarbon, dd4hep::Position(0., 0., -PosEndcapCarbon));  

  }


} // namespace

static dd4hep::Ref_t create_element( dd4hep::Detector& description, xml_h e, dd4hep::SensitiveDetector sens_det ) {

  xml_det_t x_det = e;
  CDCHBuild builder( description, x_det, sens_det );
  string    det_name = x_det.nameStr();

  dd4hep::printout( dd4hep::DEBUG, "CreateCDCH", "Detector name: %s with ID: %s", det_name, x_det.id() );

  DetElement  CDCH_det = builder.detector; // ( det_name, x_det.id() );
  dd4hep::Box CDCH_box( "5000/2", "5000/2", "5000/2" );

  Volume envelope( "lvCDCH", CDCH_box, description.air() );
  envelope.setVisAttributes( description, "vCDCH:Air" );
  PlacedVolume pv;

  dd4hep::printout( dd4hep::DEBUG, "CreateCDCH", "MotherVolume is: %s", envelope.name() );
  sens_det.setType( "tracker" );

  builder.buildVolumes( e );
  builder.placeDaughters( CDCH_det, envelope, e );

  // ******************************************************
  // Build CDCH cable
  // ******************************************************

  builder.build_layer(CDCH_det, envelope);

  // ******************************************************
  // Build CDCH cell and beam plug
  // ******************************************************

//  builder.build_cell();
//  builder.build_beamplug();

  // ******************************************************
  // Assemble CDCH
  // ******************************************************

//  builder.build_CDCH( Ecal_det, envelope );

  // ******************************************************
  // Place the CDCH in the world
  // ******************************************************

  pv = builder.placeDetector( envelope );
  pv.addPhysVolID( "system", x_det.id() );

  return CDCH_det;
}

DECLARE_DETELEMENT( IDEA_CDCH_v1_0, create_element )
