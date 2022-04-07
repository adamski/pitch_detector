/*
 
 BEGIN_JUCE_MODULE_DECLARATION
  ID:               pitch_detector
  vendor:           adamski
  version:          0.2.0
  name:             Pitch Detector
  description:      Pitch estimation methods
  website:          http://www.github.com/adamski/pitch_detector
  license:          MIT
  dependencies:     juce_core, juce_audio_basics, audio_fft
  OSXFrameworks:    
  iOSFrameworks:    
 END_JUCE_MODULE_DECLARATION
 */

#pragma once
#include <juce_core/juce_core.h>
#include <juce_audio_basics/juce_audio_basics.h>
#include <audio_fft/audio_fft.h>
#include <float.h>
#include <complex>

namespace adamski {
    
using namespace juce;

#include "source/PitchMPM.h"
#include "source/PitchYIN.h"

}
