#include "JuceHeader.h"

class PitchYIN
{

public:

    PitchYIN (unsigned int bufferSize)
    {
        yin = AudioSampleBuffer (1, bufferSize);
        tolerence = 0.15;
    }

    /** Output the difference function */
    void difference(AudioSampleBuffer input) 
    {
        float tmp;
        float *yinData = yin.getWritePointer(0);
        const float *inputData = input.getReadPointer(0);

        FloatVectorOperations::fill(yinData, 0.0, yin.getNumSamples());

        for (int tau = 1; tau < yin.getNumSamples(); tau++) 
        {
            for (int j = 0; j < yin.getNumSamples(); j++) 
            {
                tmp = inputData[j] - inputData[j + tau];
                yinData[tau] += (tmp * tmp);
            }
        }
        
    }

    /** cumulative mean normalized difference function */
    void cumulativeMean ()
    {
        float *yinData = yin.getWritePointer(0);
        float tmp = 0.;
        yinData[0] = 1.;
        //AUBIO_DBG("%f\t",yinData[0]);
        for (int tau = 1; tau < yin.getNumSamples(); tau++)
        {
            tmp += yinData[tau];
            yinData[tau] *= tau / tmp;
            //AUBIO_DBG("%f\t",yinData[tau]);
        }
        //AUBIO_DBG("\n");
    }

    int getPitch (Array<float> buffer) 
    {
        int tau = 0;
        float *yinData = yin.getWritePointer(0);
        do {
            if (yinData[tau] < 0.1) 
            {
                while (yinData[tau + 1] < yinData[tau]) 
                {
                    tau++;
                }
                return tau;
            }
            tau++;
        } while (tau < yin.getNumSamples());
        //AUBIO_DBG("No pitch found");
        return 0;
        
    }

    float calculatePitch (AudioSampleBuffer input) 
    {
        int period;
        float tmp = 0.0, tmp2 = 0.0;
        float *yinData = yin.getWritePointer(0);
        float *inputData = input.getReadPointer(0);

        yinData[0] = 1.;
        for (int tau = 1; tau < yin.getNumSamples(); tau++) {
            yinData[tau] = 0.;
            for (int j = 0; j < yin.getNumSamples(); j++) {
                tmp = inputData[j] - inputData[j + tau];
                yinData[tau] += (tmp * tmp);
            }
            tmp2 += yinData[tau];
            if (tmp2 != 0) {
                yinData[tau] *= tau / tmp2;
            } else {
                yinData[tau] = 1.;
            }
            period = tau - 3;
            if (tau > 4 && (yinData[period] < tolerence) &&
                    (yinData[period] < yinData[period + 1])) {
                return quadraticPeakPosition (yin, period);
            }
        }
        return quadraticPeakPosition (yin, minElement (yin));
    }

    float getPitchInHz (AudioSampleBuffer input)
    {
        float pitch = 0.;
        //slideBlock (input);
        //pitch = calculate (buf);
        pitch = calculate (input);
        
        if (pitch > 0) {
            pitch = sampleRate / (pitch + 0.);
        } else {
            pitch = 0.;
        }
        currentPitch = pitch;
        return pitch;
    }


private:
    AudioSampleBuffer yin, buf;
    unsigned int tolerence, confidence;
    unsigned int sampleRate;
    float currentPitch;

    /** adapter to stack ibuf new samples at the end of buf, and trim `buf` to `bufsize` */
    void slideBlock (AudioSampleBuffer ibuf)
    {
        float *bufData = buf.getWritePointer(0);
        float *ibufData = ibuf.getReadPointer(0);

        unsigned int j = 0, overlapSize = 0;
        overlapSize = buf.getNumSamples() - ibuf.getNumSamples();
        for (j = 0; j < overlapSize; j++) 
        {
            bufData[j] = bufData[j + ibuf.getNumSamples()];
        }
        for (j = 0; j < ibuf.getNumSamples(); j++) 
        {
            bufData[j + overlapSize] = ibufData[j];
        }
    }

    // Below functions should go in seperate utilities class

    float quadraticPeakPosition (AudioSampleBuffer* x, unsigned int pos) {
        float s0, s1, s2; 
        float *xData = x.getReadPointer(0);
        unsigned int x0, x2;
        if (pos == 0 || pos == x.getNumSamples() - 1) return pos;
        x0 = (pos < 1) ? pos : pos - 1;
        x2 = (pos + 1 < x.getNumSamples()) ? pos + 1 : pos;
        if (x0 == pos) return (xData[pos] <= xData[x2]) ? pos : x2;
        if (x2 == pos) return (xData[pos] <= xData[x0]) ? pos : x0;
        s0 = xData[x0];
        s1 = xData[pos];
        s2 = xData[x2];
        return pos + 0.5 * (s0 - s2 ) / (s0 - 2.* s1 + s2);
    }

    unsigned int minElement (AudioSampleBuffer * s)
    {
        float *sData = s->getReadPointer(0);
    #ifndef JUCE_USE_VDSP_FRAMEWORK
        unsigned int j, pos = 0;
        float tmp = sData[0];
        for (j = 0; j < s.getNumSamples(); j++) {
            pos = (tmp < sData[j]) ? pos : j;
            tmp = (tmp < sData[j]) ? tmp : sData[j];
        }
    #else
        float tmp = 0.0;
        unsigned int pos = 0;
        #if !DOUBLE_SAMPLES
        vDSP_minvi(sData, 1, &tmp, (vDSP_Length *)&pos, s.getNumSamples());
        #else
        vDSP_minviD(sData, 1, &tmp, (vDSP_Length *)&pos, s.getNumSamples());
        #endif
    #endif
        return pos;
    }    
}
