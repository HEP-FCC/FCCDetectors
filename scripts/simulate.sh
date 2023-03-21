### Use ddsim to make a simple particle gun test
### like in the FCC tutorial https://hep-fcc.github.io/fcc-tutorials/full-detector-simulations/Geometry/Geometry.html#modify-an-existing-xml-file

ddsim --compactFile Detector/DetFCCeeIDEA/compact/FCCee_DectMaster.xml \
      --enableGun \
      --gun.distribution uniform \
      --gun.energy "10*GeV" \
      --gun.particle mu- \
      --numberOfEvents 10000 \
      --outputFile Step1_edm4hep.root
