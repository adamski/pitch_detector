

#if JUCE_IOS || JUCE_MAC // TODO: Move this to Projucer project?
#define AUDIOFFT_APPLE_ACCELERATE 1
#endif

#define CUTOFF 0.93f //0.97 is default
#define SMALL_CUTOFF 0.5f
#define LOWER_PITCH_CUTOFF 80.f //hz

/**
 * TODO: Provide switch between time-based and FFT based methods
 */

class PitchMPM
{
public:
    PitchMPM (size_t detectionBufferSize) : PitchMPM (44100, detectionBufferSize) {}

    PitchMPM (int detectionSampleRate, size_t detectionBufferSize) : bufferSize (detectionBufferSize),
                                                                     sampleRate ((float) detectionSampleRate),
                                                                     fftSize (2 * bufferSize), // Needs to be a power of 2!
                                                                     real (audiofft::AudioFFT::ComplexSize (fftSize)),
                                                                     imag (audiofft::AudioFFT::ComplexSize (fftSize)),
                                                                     output (fftSize)

    {
        //nsdf.insertMultiple(0, 0.0, bufferSize);
    }

    ~PitchMPM()
    {
        //nsdf.clear();
//        maxPositions.clear();
//        ampEstimates.clear();
//        periodEstimates.clear();
    }

    float getPitch (const float* audioBuffer)
    {
//        periodEstimates.clearQuick();
//        ampEstimates.clearQuick();

        if (audioBuffer == nullptr)
        {
            DBG ("audioBuffer NULL");
            return 0.0f;
        }

        std::vector<float> nsdf = nsdfFrequencyDomain (audioBuffer);
        std::vector<int> maxPositions = peakPicking (nsdf);
        std::vector<std::pair<float, float>> estimates;

        float highestAmplitude = -FLT_MAX;

        for (auto tau : maxPositions)
        {
            auto tauPosition = size_t (tau);
            highestAmplitude = std::max (highestAmplitude, nsdf[tauPosition]);

            if (nsdf[tauPosition] > SMALL_CUTOFF)
            {
                auto x = parabolicInterpolation (nsdf, float (tauPosition));
                estimates.push_back (x);
                highestAmplitude = std::max (highestAmplitude, std::get<1> (x));
            }
        }

        if (estimates.empty())
            return -1;

        float actualCutoff = CUTOFF * highestAmplitude;
        float period = 0;

        for (auto estimate : estimates)
        {
            if (std::get<1> (estimate) >= actualCutoff)
            {
                period = std::get<0> (estimate);
                break;
            }
        }

        float pitchEstimate = (sampleRate / period);
        return (pitchEstimate > LOWER_PITCH_CUTOFF) ? pitchEstimate : -1;
    }

    void setSampleRate (int newSampleRate)
    {
        sampleRate = float (newSampleRate);
    }

    void setBufferSize (int newBufferSize)
    {
        bufferSize = size_t (newBufferSize);
        input.resize (bufferSize);
        fftSize = 2 * bufferSize;
        real.resize (audiofft::AudioFFT::ComplexSize (fftSize));
        imag.resize (audiofft::AudioFFT::ComplexSize (fftSize));
        output.resize (fftSize);
    }

private:
    size_t bufferSize;
    float sampleRate;

    audiofft::AudioFFT fft;
    size_t fftSize;
    std::vector<float> input;
    std::vector<float> real;
    std::vector<float> imag;
    std::vector<float> output;

    juce::Array<float> periodEstimates;
    juce::Array<float> ampEstimates;

    inline std::pair<float, float> parabolicInterpolation (std::vector<float> array, float x)
    {
        int xAdjusted;
        auto xPosition = size_t (x);

        if (x < 1)
        {
            xAdjusted = int ((array[xPosition] <= array[xPosition + 1]) ? x : x + 1);
        }
        else if (x > signed (array.size()) - 1)
        {
            xAdjusted = int ((array[xPosition] <= array[xPosition - 1]) ? x : x - 1);
        }
        else
        {
            float den = array[xPosition + 1] + array[xPosition - 1] - 2 * array[xPosition];
            float delta = array[xPosition - 1] - array[xPosition + 1];
            return (den == 0.0f) ? std::make_pair (x, array[xPosition]) : std::make_pair (x + delta / (2 * den), array[xPosition] - delta * delta / (8 * den));
        }
        return std::make_pair (float (xAdjusted), array[size_t (xAdjusted)]);
    }

    static std::vector<int> peakPicking (std::vector<float> nsdf)
    {
        std::vector<int> max_positions {};
        int pos = 0;
        int curMaxPos = 0;
        juce::ssize_t size = juce::ssize_t (nsdf.size());

        while (pos < (size - 1) / 3 && nsdf[size_t (pos)] > 0)
            pos++;
        while (pos < size - 1 && nsdf[size_t (pos)] <= 0.0)
            pos++;

        if (pos == 0)
            pos = 1;

        while (pos < size - 1)
        {
            if (nsdf[size_t (pos)] > nsdf[size_t (pos - 1)] && nsdf[size_t (pos)] >= nsdf[size_t (pos + 1)])
            {
                if (curMaxPos == 0)
                {
                    curMaxPos = pos;
                }
                else if (nsdf[size_t (pos)] > nsdf[size_t (curMaxPos)])
                {
                    curMaxPos = pos;
                }
            }
            pos++;
            if (pos < size - 1 && nsdf[size_t (pos)] <= 0)
            {
                if (curMaxPos > 0)
                {
                    max_positions.push_back (curMaxPos);
                    curMaxPos = 0;
                }
                while (pos < size - 1 && nsdf[size_t (pos)] <= 0.0)
                {
                    pos++;
                }
            }
        }
        if (curMaxPos > 0)
        {
            max_positions.push_back (curMaxPos);
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
    std::vector<float> nsdfFrequencyDomain (const float* audioBuffer)
    {
        //std::vector<std::complex<float>> acf(size2);
        //std::vector<float> acf_real{};

        real.resize (fftSize);
        imag.resize (fftSize);

        if (audioBuffer == nullptr)
        {
            DBG ("audioBuffer NULL: nsdfFrequencyDomain");
        }

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

    std::vector<float> autoCorrelation (const float* audioBuffer)
    {
        if (audioBuffer == nullptr)
        {
            DBG ("audioBuffer NULL: autoCorrelation");
        }

        //AudioSampleBuffer paddedAudioBuffer (audioBuffer, 1, fftSize);
        std::vector<float> inputBuf (audioBuffer, audioBuffer + bufferSize);
        inputBuf.resize (fftSize, 0.0f);

        if (audioBuffer == nullptr)
        {
            DBG ("audioBuffer NULL: autoCorrelation post resize");
        }

        if (inputBuf.data() == nullptr)
        {
            DBG ("inputBuf.data() NULL: autoCorrelation post resize");
        }

        fft.init (fftSize);
        fft.fft (inputBuf.data(), real.data(), imag.data());
        //fft.fft(audioBuffer, real.data(), imag.data());

        // Complex Conjugate
        for (int i = 0; i < int (fftSize); ++i)
        {
            /**
             * std::complex method
             */
            std::complex<float> complex (real[size_t (i)], imag[size_t (i)]);
            complex = complex * std::conj (complex); // No need to scale as AudioFFT does this already
            real[size_t (i)] = complex.real();
            imag[size_t (i)] = complex.imag();

            /**
             * calculate via real[i] * real[i] + imag[i] * imag[i].
             * And if you really mean complex conjugation, just negate imag[i]
             */

            //imag[i] *= -1;
            //real[i] = real[i] * real[i]; // + imag[i] * imag[i];
        }

        fft.ifft (output.data(), real.data(), imag.data());
        return output;
    }
};
