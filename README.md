
# Done:
* Comment continues to not skip the candidate with the muon matches.

# Things to add trigger efficiency measurement.
* Add a new trigger (so far maybe not necesary)
* Run for the whole EGamma dataset.
* Check out all the cuts
  * For the dimuons
  ```
    preVtxSelection        = cms.string(' && '.join([
        'abs(userCand("mu1").bestTrack.dz - userCand("mu2").bestTrack.dz) <= 0.4 ',
        'mass() > 2',
        'mass() < 4',
        'userFloat("muons12_deltaR") > 0.01',
        #'mass() > 0.0',
        #'mass() < 5.0',
        #'userFloat("muons12_deltaR") > 0.03',
        ])
    ),
    postVtxSelection   = cms.string(
        'userFloat("sv_prob") > 1.e-5 '
        '&& userFloat("sv_success") > 0. '
        '&& pt > 3 '
        #        '&& userFloat("sv_chi2") < 998 ' 
    ),
  ```
  * For the three muons system:
  ```
    preVtxSelection       = cms.string(''),
    postVtxSelection      = cms.string(' && '.join([
        'userFloat("fitted_cos_theta_2D") >= 0',
        'mass < 10.',
        'userInt("sv_OK") == 1',
        'userFloat("sv_prob") > 1e-8',
        ])
    ),
  ```

# Output file info
- With the trigger flag selection
```shell
-rw-r--r--. 1 garamire zh 9.4M May 31 19:15 RJPsi_mc_10614.root
```
- Without the trigger flags selection
```shell
-rw-r--r--. 1 garamire zh 9.5M May 31 19:43 RJPsi_mc_10614.root
```

- For the first file of the EGamma dataset:

- Without the trigger flags selection
```shell
-rw-r--r--. 1 garamire zh  77M May 31 20:41 RJPsi_data_10614.root
```


# nanoAOD producer for RJPsi analysis

Add your fork of the repository as remote and relevant packages:

```shell
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src/
cmsenv
git cms-init
git cms-merge-topic -u friti:TransientTracks
git cms-merge-topic -u friti:KinParticleVtxFitter
git clone git@github.com:friti/RJPsiAnalysis.git  ./PhysicsTools
git cms-merge-topic -u friti:GenParticlesPrecision
git cms-addpkg PhysicsTools/NanoAOD
scram b
```
