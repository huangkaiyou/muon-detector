The code is to analyze the data from muon detector, which is composed of scintillator, PMT, and Oscilloscope.
In the setup, we use two layer of setup.
When muon hit the setup from atmosphere, there is first peak in both channel
When the muon decay into electron, the signal would only exist in second(lower) channel
the raw data are filtered into new_lifetimedatatxt with required condition for first peak.
Main_txt_scipy_stable_sigma code is uesd to select the possible muon event with second peak judgement.
Lifetime_dis code plot the result into lifetime distribution, which get the tau from fit parameter.
