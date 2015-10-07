class PitchYIN 
{
private:
    Array<float> yin;
    unsigned int tolerence, confidence;

    PitchYIN ()
    {
        tolerence = 0.15;
    }

    /** Output the difference function */
    void difference(Array<float> input) 
    {
        float tmp;
        // fill with zeros
        for (int tau = 0; tau < yin.size(); tau++) 
        {
            yin.setUnchecked(tau, 0.0);
        }

        for (int tau = 1; tau < yin.size(); tau++) 
        {
            for (int j = 0; j < yin.size(); j++) 
            {
                tmp = input.getUnchecked (j) - input.getUnchecked (j + tau);
                yin.setUnchecked (yin.getUnchecked (tau) + (tmp * tmp));
            }
        }
        
    }

    /** cumulative mean normalized difference function */
    void cumulativeMean (Array<float> yin)
        {
            float tmp = 0.;
            yin.setUnchecked(0) = 1.0;
            //AUBIO_DBG("%f\t",yin->data[0]);
            for (int tau = 1; tau < yin.size(); tau++) 
            {
                tmp += yin.getUnchecked (tau);
                yin.setUnchecked (tau, yin.getUnchecked (tau) * tau / tmp);
                //AUBIO_DBG("%f\t",yin->data[tau]);
            }
            //AUBIO_DBG("\n");
        }

    int getPitch (Array<float> buffer) 
    {
        do {
            if (yin.getUnchecked (tau) < 0.1) {
                while (yin.getUnchecked (tau + 1) < yin.getUnchecked (tau)) {
                    tau++;
                }
                return tau;
            }
            tau++;
        } while (tau < yin.size());
        //AUBIO_DBG("No pitch found");
        return 0;
        
    }

    void do (Array<float> input, Array<float> output) 
    {
        int period;
        float tmp = 0.0, tmp2 = 0.0;
        yin.setUnchecked(0) = 1.0;
        for (int tau = 1; tau < yin.size(); tau++) {
            yin.setUnchecked (tau, 0.0);
            for (j = 0; j < yin.size(); j++) {
                tmp = input.getUnchecked(j) - input.getUnchecked(j + tau);
                yin.setUnchecked (yin.getUnchecked (tau) + (tmp * tmp));            }
            tmp2 += yin.getUnchecked (tau);
            if (tmp2 != 0) {
                yin.setUnchecked (tau, yin.getUnchecked (tau) * tau / tmp2);
            } else {
                yin.setUnchecked (tau, 1.0);
            }
            period = tau - 3;
            if (tau > 4 && (yin.getUnchecked (period) < tolerence) &&
                    (yin.getUnchecked (period) < yin.getUnchecked (period)+ 1)) {
                out.setUnchecked(0) = fvec_quadratic_peak_pos (yin, period);
                return;
            }
        }
        out.setUnchecked(0) = fvec_quadratic_peak_pos (yin, fvec_min_elem (yin));
    }

    float fvec_quadratic_peak_pos (Array<float> * x, unsigned int pos) {
        float s0, s1, s2; 
        unsigned int x0, x2;
        if (pos == 0 || pos == x->length - 1) return pos;
        x0 = (pos < 1) ? pos : pos - 1;
        x2 = (pos + 1 < x->length) ? pos + 1 : pos;
        if (x0 == pos) return (x->data[pos] <= x->data[x2]) ? pos : x2;
        if (x2 == pos) return (x->data[pos] <= x->data[x0]) ? pos : x0;
        s0 = x->data[x0];
        s1 = x->data[pos];
        s2 = x->data[x2];
        return pos + 0.5 * (s0 - s2 ) / (s0 - 2.* s1 + s2);
    }

    uint_t fvec_min_elem (fvec_t * s)
    {
#ifndef HAVE_ACCELERATE
        uint_t j, pos = 0.;
        smpl_t tmp = s->data[0];
        for (j = 0; j < s->length; j++) {
            pos = (tmp < s->data[j]) ? pos : j;
            tmp = (tmp < s->data[j]) ? tmp : s->data[j];
        }
#else
        smpl_t tmp = 0.;
        uint_t pos = 0.;
#if !HAVE_AUBIO_DOUBLE
        vDSP_minvi(s->data, 1, &tmp, (vDSP_Length *)&pos, s->length);
#else
        vDSP_minviD(s->data, 1, &tmp, (vDSP_Length *)&pos, s->length);
#endif
#endif
        return pos;
    }    

}
