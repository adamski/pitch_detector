# Pitch Detector
JUCE module for pitch estimation

PitchYIN class based on the YIN implementation found in the aubio library

PitchMPM class adapted from the McLeod Pitch Method implementation in https://github.com/sevagh/pitch-detection

The updated version of the PitchMPM class now uses FFT for the Normalised Squared Difference Function using the AudioFFT library (via the module wrapper at https://github.com/adamski/audio_fft). The previous time-based version is now in the `time-based` branch. 

### TODO
- Seperate time-based method into another class that can be used as an alternative to the FFT based method
- Add FFT based YIN implementation (not a priority, MPM works well for my needs - PR's welcome)
- Create base (virtual) `Pitch` class and add implementations as subclasses.  
- Add other methods, e.g. Wavelet? 
- Remove JUCE dependency from implementations so that they can be used on embedded platforms, e.g. Arduino/Teensy. Will also need 'pluggable' FFT methods for this to work. 
