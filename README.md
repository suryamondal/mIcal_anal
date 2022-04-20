# mIcal reconstruction : Kalman

This code is Geant4 simulation for mIcal with 10 layers of RPCs. Total 10 RPCs are placed in the center of each layer.

The code is tested for `root6.20.04`, `clhep2404` and `geant4.10.04.p03`.

The field file is `B_mical_hist.root`, the file name is hard-coded in `src/FieldPropagator.cc`. The field map could be scaled down. Please look for `fieldxin->Scale(` in the same file.


*Warning:*
- The code has memory leakage. It may overflow memory. I do not recommend it for running more than 1M events.
- One might need to comment out `G4DataQuestionaire` portions in `src/micalPhysicsList.cc` for later geant4 versions.

Consult `batch_execute` to avoid the scenario.

The name of the geometry file `geo_mical_world.gdml` is hard-coded in `src/anal_ical.cc`.

Please copy the file from `mIcal_mc` if running for the first time.

Source in sim01: `source env.sh`

Requirement only once:
```
rm src/Hitsdict.cc
rm src/HitPosdict.cc
rm src/HitPosdict_rdict.pcm
rm src/Hitsdict_rdict.pcm
cd include
rootcint ../src/HitPosdict.cc -c HitPos.h
rootcint ../src/Hitsdict.cc -c Hits.h
cd ..
```

Compile: `mkdir -p build && cd build && cmake3 .. && make -j 4 && cd ..`

Run: `./anal_ical test.log 0 0.1 1000` where `test.log` contains the list of input files and for the other arguments, please constult with `src/anal_ical.cc`.
