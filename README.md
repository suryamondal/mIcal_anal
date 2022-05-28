# mIcal reconstruction : Kalman

This code is for Kalman based reconstruction of digi files generated by `mIcal_mc`.

The code is tested for `root6.20.04`, `clhep2404` and `geant4.10.04.p03`.

The field file is `B_mical_hist.root`, the file name is hard-coded in `src/FieldPropagator.cc`. The field map could be scaled down. Please look for `fieldxin->Scale(` in the same file.


*Warning:*
- The code has memory leakage. It may overflow memory. I do not recommend it for running more than 1M events.
- Consult `batch_execute` to avoid the scenario.

The name of the geometry file `geo_mical_world.gdml` is hard-coded in `src/anal_ical.cc`.

Please copy the file from `mIcal_mc` if running for the first time.

Source in sim01: `source env.sh`

Requirement only once:
```
rm src/Hitsdict.cc
rm src/HitPosdict.cc
rm src/HitPosdict_rdict.pcm
rm src/Hitsdict_rdict.pcm
cd inc
rootcint ../src/HitPosdict.cc -c HitPos.h
rootcint ../src/Hitsdict.cc -c Hits.h
cd ..
```

Compile: `mkdir -p build && cd build && cmake3 .. && make -j 4 && cd ..`

Create necessary directories in needed: `simdata`, `recodata`, `digidata`, etc.

Run: `./anal_ical test.log 1 0.1 1 1000` where `test.log` contains the list of input files and for the other arguments, please constult with `src/anal_ical.cc`.

InputOutput:
```
0: SIM  -> DIGI
1: DIGI -> RECO
2: DATA -> RECO
```

Format of `test.log`: <filename><no of events><start event no>

*Warning:*
- For `SIM to DIGI`, please check the `mIcal_mc`. The colleted file and `Cal0SD` has to be updated. 
- For `DATA to RECO`, please check the input tree if needed. 
  
  ## Tracking flow
```mermaid
  graph TD;
  id0[anal_ical]-->id1[InoRecoAlg.ReadEvent]-->id2[InoRecoAlg.PerformTrackReconstruction];
  id2-->id3[InoTrackFinder.RunTheFinder]-->id4[InoTrackFitAlg.RunAlg];
  id4-->id5[InitialFramework_new]-->id6[RunTheFitter_new];
  id6-->id8[GetInitialCovarianceMatrix]-->id12[StoreFilteredData]-->id9[GoBackwords_new true]-->id10[ResetCovarianceMatrix];
  id10-->id11[RemoveTrkHitsInShw]-->id14[StoreFilteredData]-->id13[GoForwards_new false]-->id15[ShowerSwim]-->id16[ResetCovarianceMatrix];
```

  
## Backup
  ### In RunAlg
  fFinderTrack is looped over ptrackCollection->InoTrack_list.
  - InitialFramework_new -> TrkClustsData and InitTrkClustData are set here from fFinderTrack
  - RunTheFitter_new ->
    - GetInitialCovarianceMatrix(true)
    - StoreFilteredData(MaxPlane)
    - GoBackwards_new(false)
    - ResetCovarianceMatrix()
    - FilteredData[i].clear()
    - StoreFilteredData(MinPlane)
    - GoForwards_new(false)
    - ResetCovarianceMatrix()
    - Iteration L932
      - GoBackwards_new(true)
      - ResetCovarianceMatrix()
      - StoreFilteredData(MinPlane)
      - GoForwards_new(false)
      - ResetCovarianceMatrix() L1037
      - Swim(StateVector, Prediction, MaxPlane, loc_zend, GoForward) and then break
    - FillGapsInTrack()
    - SetTrackProperties(Prediction)
