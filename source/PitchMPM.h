

#if JUCE_IOS || JUCE_MAC // TODO: Move this to Projucer project?
#define AUDIOFFT_APPLE_ACCELERATE 1
#endif

#define CUTOFF 0.93 //0.97 is default
#define SMALL_CUTOFF 0.5
#define LOWER_PITCH_CUTOFF 80 //hz

/**
 * TODO: Provide switch between time-based and FFT based methods
 */

class PitchMPM
{

public:

    PitchMPM(size_t bufferSize) : PitchMPM(44100, bufferSize) {}

    PitchMPM(int sampleRate, size_t bufferSize) : bufferSize (bufferSize),
                                                  sampleRate (sampleRate),
                                                  fftSize (2 * bufferSize), // Needs to be a power of 2!
                                                  real (audiofft::AudioFFT::ComplexSize(fftSize)),
                                                  imag (audiofft::AudioFFT::ComplexSize(fftSize)),
                                                  output (fftSize)


    {
        //nsdf.insertMultiple(0, 0.0, bufferSize);

    }
    
    
    ~PitchMPM()
    {
        
        //nsdf.clear();
        maxPositions.clear();
        ampEstimates.clear();
        periodEstimates.clear();
        
    }
    
    float getPitch(const float *audioBuffer)
    {
        
        maxPositions.clearQuick();
        periodEstimates.clearQuick();
        ampEstimates.clearQuick();
        
        //nsdfTimeDomain(audioBuffer);
        //nsdfFrequencyDomain(audioBuffer);
        
        if (audioBuffer == nullptr)
        {
            DBG ("audioBuffer NULL");
            return 0.0f;
        }
        //nsdf = Array<float> (nsdfFrequencyDomain(audioBuffer).data());
        std::vector<float> _nsdf = nsdfFrequencyDomain(audioBuffer);
        std::vector<int> max_positions = peak_picking(_nsdf);
        std::vector<std::pair<float, float>> estimates;

        //peakPicking();
        
        float highestAmplitude = -FLT_MAX;
        
        for (auto tau : max_positions)
        {
            highestAmplitude = jmax(highestAmplitude, _nsdf[tau]);
            
            if (_nsdf[tau] > SMALL_CUTOFF)
            {
                auto x = parabolic_interpolation(_nsdf, tau);
                estimates.push_back(x);
                highestAmplitude = std::max(highestAmplitude, std::get<1>(x));
            }
        }

        if (estimates.empty()) return -1;

        float actualCutoff = CUTOFF * highestAmplitude;
        float period = 0;

        for (auto estimate : estimates)
        {
            if (std::get<1>(estimate) >= actualCutoff)
            {
                period = std::get<0>(estimate);
                break;
            }
        }

        float pitchEstimate = (sampleRate / period);
        return (pitchEstimate > LOWER_PITCH_CUTOFF) ? pitchEstimate : -1;
    }
    
    void setSampleRate (int newSampleRate)
    {
        sampleRate = newSampleRate;
    }

    void setBufferSize (int newBufferSize)
    {
        bufferSize = newBufferSize;
        input.resize (bufferSize);
        fftSize = 2 * bufferSize;
        real.resize (audiofft::AudioFFT::ComplexSize(fftSize));
        imag.resize (audiofft::AudioFFT::ComplexSize(fftSize));
        output.resize (fftSize);
    }

private:
    size_t bufferSize;

    audiofft::AudioFFT fft;
    size_t fftSize;
    std::vector<float> input;
    std::vector<float> real;
    std::vector<float> imag;
    std::vector<float> output;

    float sampleRate;
    
    float turningPointX, turningPointY;
    //Array<float> nsdf;
    
    Array<int> maxPositions;
    Array<float> periodEstimates;
    Array<float> ampEstimates;
    
    /*
    void parabolicInterpolation(int tau)
    {
        float nsdfa = nsdf.getUnchecked (tau - 1);
        float nsdfb = nsdf.getUnchecked (tau);
        float nsdfc = nsdf.getUnchecked (tau + 1);
        float bValue = tau;
        float bottom = nsdfc + nsdfa - 2 * nsdfb;
        if (bottom == 0.0)
        {
            turningPointX = bValue;
            turningPointY = nsdfb;
        }
        else
        {
            float delta = nsdfa - nsdfc;
            turningPointX = bValue + delta / (2 * bottom);
            turningPointY = nsdfb - delta * delta / (8 * bottom);
        }
    }
    
    void peakPicking()
    {
        
        int pos = 0;
        int curMaxPos = 0;
        float* nsdfPtr = nsdf.getRawDataPointer();
        
        while (pos < (bufferSize - 1) / 3 && nsdfPtr[pos] > 0) {
            pos++;
        }
        
        while (pos < bufferSize - 1 && nsdfPtr[pos] <= 0.0) {
            pos++;
        }

        if (pos == 0) {
            pos = 1;
        }
        
        while (pos < bufferSize - 1) {
            if (nsdfPtr[pos] > nsdfPtr[pos - 1] && nsdfPtr[pos] >= nsdfPtr[pos + 1]) {
                if (curMaxPos == 0) {
                    curMaxPos = pos;
                } else if (nsdfPtr[pos] > nsdfPtr[curMaxPos]) {
                    curMaxPos = pos;
                }
            }
            pos++;
            if (pos < bufferSize - 1 && nsdfPtr[pos] <= 0) {
                if (curMaxPos > 0) {
                    maxPositions.add (curMaxPos);
                    curMaxPos = 0;
                }
                while (pos < bufferSize - 1 && nsdfPtr[pos] <= 0.0f) {
                    pos++;
                }
            }
        }
        if (curMaxPos > 0)
        {
            maxPositions.add (curMaxPos);
        }
    }
     */

    inline std::pair<float, float> parabolic_interpolation(std::vector<float> array, float x)
    {
        int x_adjusted;

        if (x < 1) {
            x_adjusted = (array[x] <= array[x+1]) ? x : x+1;
        } else if (x > signed(array.size())-1) {
            x_adjusted = (array[x] <= array[x-1]) ? x : x-1;
        } else {
            float den = array[x+1] + array[x-1] - 2 * array[x];
            float delta = array[x-1] - array[x+1];
            return (!den) ? std::make_pair(x, array[x]) : std::make_pair(x + delta / (2 * den), array[x] - delta*delta/(8*den));
        }
        return std::make_pair(x_adjusted, array[x_adjusted]);
    }
    
    static std::vector<int> peak_picking(std::vector<float> nsdf)
    {
        std::vector<int> max_positions{};
        int pos = 0;
        int curMaxPos = 0;
        ssize_t size = nsdf.size();
        
        while (pos < (size - 1) / 3 && nsdf[pos] > 0) pos++;
        while (pos < size - 1 && nsdf[pos] <= 0.0) pos++;
        
        if (pos == 0) pos = 1;
        
        while (pos < size - 1) {
            if (nsdf[pos] > nsdf[pos - 1] && nsdf[pos] >= nsdf[pos + 1]) {
                if (curMaxPos == 0) {
                    curMaxPos = pos;
                } else if (nsdf[pos] > nsdf[curMaxPos]) {
                    curMaxPos = pos;
                }
            }
            pos++;
            if (pos < size - 1 && nsdf[pos] <= 0) {
                if (curMaxPos > 0) {
                    max_positions.push_back(curMaxPos);
                    curMaxPos = 0;
                }
                while (pos < size - 1 && nsdf[pos] <= 0.0) {
                    pos++;
                }
            }
        }
        if (curMaxPos > 0) {
            max_positions.push_back(curMaxPos);
        }
        return max_positions;
    }
    
    /*
    void nsdfTimeDomain(const float *audioBuffer)
    {
        int tau;
        for (tau = 0; tau < bufferSize; tau++) {
            float acf = 0;
            float divisorM = 0;
            for (int i = 0; i < bufferSize - tau; i++) {
                acf += audioBuffer[i] * audioBuffer[i + tau];
                divisorM += audioBuffer[i] * audioBuffer[i] + audioBuffer[i + tau] * audioBuffer[i + tau];
            }
            nsdf.setUnchecked(tau, 2 * acf / divisorM);
        }
    }
    */

    // FFT based methods
    std::vector<float> nsdfFrequencyDomain (const float *audioBuffer)
    {
        //std::vector<std::complex<float>> acf(size2);
        //std::vector<float> acf_real{};
    
        real.resize (fftSize);
        imag.resize (fftSize);

        if (audioBuffer == nullptr)
            DBG ("audioBuffer NULL: nsdfFrequencyDomain");
        
        std::vector<float> acf (autoCorrelation (audioBuffer));

        /*
        for (auto it = acf.begin() + size2/2; it != acf.end(); ++it)
            acf_real.push_back((*it) / acf[size2 / 2]);
        */
        
        /** This code is for interleaved data, above is not
        for (auto it = acf.begin() + size2/2; it != acf.end(); ++it)
            acf_real.push_back((*it)/acf[size2/2]);
            //acf_real.push_back((*it).real()/acf[size2/2].real());
         ****************************************************/

//        for (int i = size2/2; i < acf.size(); ++i)
//            nsdf.setUnchecked(i, acf[i] / acf[size2 / 2]);
        

        //return acf_real;
        return acf;
    }

    std::vector<float> autoCorrelation(const float *audioBuffer)
    {
        if (audioBuffer == nullptr)
            DBG ("audioBuffer NULL: autoCorrelation");
        
        //AudioSampleBuffer paddedAudioBuffer (audioBuffer, 1, fftSize);
        std::vector<float> input (audioBuffer, audioBuffer + bufferSize);
        input.resize(fftSize, 0.0f);
        
        if (audioBuffer == nullptr)
            DBG ("audioBuffer NULL: autoCorrelation post resize");
        
        if (input.data() == nullptr)
            DBG ("input.data() NULL: autoCorrelation post resize");

        fft.init(fftSize);
        fft.fft(input.data(), real.data(), imag.data());
        //fft.fft(audioBuffer, real.data(), imag.data());

        // Complex Conjugate
        for (int i = 0; i < fftSize; ++i)
        {
            /**
             * std::complex method
             */
            std::complex<float> complex(real[i], imag[i]);
            complex = complex * std::conj(complex); // No need to scale as AudioFFT does this already
            real[i] = complex.real();
            imag[i] = complex.imag();

            /**
             * calculate via real[i] * real[i] + imag[i] * imag[i].
             * And if you really mean complex conjugation, just negate imag[i]
             */

            //imag[i] *= -1;
            //real[i] = real[i] * real[i]; // + imag[i] * imag[i];
        }

        fft.ifft(output.data(), real.data(), imag.data());
        return output;
    }



};
