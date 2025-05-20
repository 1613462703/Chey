#pragma once
#include <stdio.h>
#include <complex>

//#define ARMA_DONT_USE_CXX11
#define ARMA_DONT_USE_CXX11_MUTEX
#include <armadillo>


#ifndef Chebyshev_Filter_H
#define Chebyshev_Filter_H

#ifdef __cplusplus
extern "C" {  // only need to export C interface if
              // used by C++ source code
#endif

    namespace CH_FI
    {

        class Chebyshev_Filter

        {

        private:

            //get analog, pre - warped frequencies
            void freq_pre_wrapped(int, double, double);

            //convert to low-pass prototype estimate
            void Wn_f1_Wn_f2(int, double, double);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<std::complex<double>> poly(std::vector<std::complex<double>>, int);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<double> poly(arma::mat, int);

            //Calculate the coefficients of the polynomial (based on Matlab code)
            std::vector<double> poly(arma::cx_vec, int);

            //Sort complex numbers into complex conjugate pairs, starting from with the number with the lowest real part.
            //If the real part is the same for all the number, order according to the absolute value of the highest imaginary part.
            //Within a pair, the element with negative imaginary part comes first. 
            std::vector<std::complex<double>> cplxpair(std::vector<std::complex<double>>);

            //Get N - th order Chebyshev type-I lowpass analog prototype 
            void cheb1ap(int, double);

            //Transform to state-space
            void zp2ss();

            //Bilinear transformation to find discrete equivalent
            void bilinear(arma::mat, arma::mat, arma::mat, arma::mat, double, int, int);

            //Transform to zero - pole - gain and polynomial forms
            void zero_pole_gain(arma::mat, int, int, double, double, double, int);
            

        public:

            //Estimate the coeffients of a band-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2bp(int, double, double, double);

            //Estimate the coeffients of a band-stop filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2bs(int, double, double, double);

            //Estimate the coeffients of a low-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2lp(int, double, double);

            //Estimate the coeffients of a high-pass filter and return a 2 rows x N coefficients matrix. Row 1 = Numerator; Row 2 = Denumerator
            std::vector<std::vector<double> > lp2hp(int, double, double);

            //Check the stability of the filter. Returns "true" is the filter is stable, false if it is unstable 
            bool check_stability_iir(std::vector<std::vector<double> >);

            //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
            void Filter_Data(double pre_filt_signal[], int size, double filt_signal[], std::vector<std::vector<double>> coeff_filt);

            std::vector<std::vector<double>> Calculate_Chebyshev(double f1, double f2, double sf, double Rp, int Oreder_filter, int typefilt);

            void Cheby1(double Input_Array[], int Number_of_Elements, double Sampling_Frequency, double Cutoff_Frequency, double Ripple_db, int Order, double Output_Array[], int type_filt);

            void Cheby1band(double Input_Array[], int Number_of_Elements, double Sampling_Frequency, double Lower_Cutoff_Freq, double Upper_Cutoff_Freq, double Ripple_db, int Order, double Output_Array[], int type_filt);

            class StreamInstance {
            public:
                StreamInstance(const std::vector<std::vector<double>>& coeff);
                double step(double x);
                void reset();
            private:
                std::vector<double> b, a;
                std::vector<double> w;
                int order;
            };
        };

    }

#endif

#ifdef __cplusplus

}

#endif

