# mucoll

Studies for the Muon Collider. Currently set up for MAIA simulations. If you're here from a poster/presentation, read [this](https://mcd-wiki.web.cern.ch/software/howto/maia/) to get up to date on how to generally run MAIA sims. I'll just point out a couple of things that are specific to my setup.

## Gen
Inside `/jobs` is the machinery to generate `.slcio` files to feed into sim and reco.

## Sim/Reco 
Run sh `.enterApptainer.sh` to enter the image required for sim and reco. Running `source setup_shower_profiles.sh` will encode all the required dependencies. From there you can run `ddsim` and `k4reco` as normal.

## Analysis
All of the analysis code is stored in `/analysisCode`. It is a mix of Python and Julia.
