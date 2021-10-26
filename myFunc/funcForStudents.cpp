#include "funcForStudents.h"
#include "complex.h"
#include <iostream>


int min(int a, int b) {
return (a <= b) ? a : b;
}
int max(int a, int b) {
return (a >= b) ? a : b;
}

void xcorr(std::vector<complex> tr1, std::vector<complex> tr2, std::vector<complex> corp, 
int shift, int shift_zero,int window, int demean, int normalize, int ndat1,
int ndat2, int ndat1d, int ndat2d)
 {
	int a, a2, b, b2, bmin, bmax, flag = 0, ind1, ind2, ind3, ind4;
	complex sum, sum1, sum2, cmax;
	std::vector<complex> tra1(ndat1);
	for (int i=0;i<ndat1;i++)
	    tra1[i] =0;
	std::vector<complex> tra2(ndat2);
	for (int i=0;i<ndat2;i++)
	    tra2[i] =0;
	    // Set standard values 
    if (window == 0) {
        window = min(ndat1, ndat2);
    }
    if (ndat1d == 0) {
        ndat1d = window;
    }
    if (ndat2d == 0) {
        ndat2d = window;
    }

    ind1 = max(0, (ndat1 - window) / 2);
    ind2 = min(ndat1, (ndat1 + window) / 2);
    ind3 = max(0, (ndat2 - window) / 2);
    ind4 = min(ndat2, (ndat2 + window) / 2);*/
    printf(" %d %d %d %d \n",ind1,ind2,ind3,ind4);

    // Demean data (Zero offset) 
    if (demean > 0) {
        sum = 0;
        for (a = ind1; a < ind2; a++) {
            sum += tr1[a];
        }
        sum /= ndat1d;
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a] - (float) sum;
        }
        if (sum == 0.0)
            flag = 1;
        sum = 0;
        for (a = ind3; a < ind4; a++) {
            sum += tr2[a];
        }
        sum /= ndat2d;
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a] - (float) sum;
        }
        if (sum == 0.0)
            flag += 1;
    } else {
        for (a = 0; a < ndat1; a++) {
            tra1[a] = tr1[a];
        }
        for (a = 0; a < ndat2; a++) {
            tra2[a] = tr2[a];
        }
    }
	
	cout <<flag;
   
    if (flag == 0) {
        a = 0;
        a2 = -shift_zero - shift;
        if (ndat1 != ndat2) {
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2 + (ndat2 - ndat1) / 2);
                bmax = min(ndat2, -a2 + (ndat1 + ndat2) / 2);
                b2 = bmin;
                b = b2 + (ndat1 - ndat2) / 2;
                corp[a] = 0;
                //printf("%d - %d %d - %d %d\n", a2, b + a2, b + a2 + bmax - bmin, b2, bmax);
                if (bmin >= bmax) {
                    continue;
                }
                for (; b2 < bmax; b++, b2++) {
                    corp[a] += tra1[b + a2] * tra2[b2];
                }
            }
        } else { // same as above only ndat2 = ndat1, we need one variable less in the second loop
            for (; a < (2 * shift + 1); a++, a2++) {
                bmin = max(0, -a2);
                bmax = min(ndat1, -a2 + ndat1);
                corp[a] = 0;
                if (bmin >= bmax) {
                    continue;
                }
                for (b = bmin; b < bmax; b++) {
                    corp[a] += tra1[b + a2] * tra2[b];
                }
            }
        }

        // normalize xcorr function 
        if (normalize > 0) {
            sum1 = sum2 = 0.0;
            for (a = ind1; a < ind2; a++) {
                sum1 += (*(tra1 + a))*(*(tra1 + a));
            }
            for (a = ind3; a < ind4; a++) {
                sum2 += (*(tra2 + a))*(*(tra2 + a));
            }
            sum1 = sqrt(sum1);
            sum2 = sqrt(sum2);
            cmax = 1 / (sum1 * sum2);
            for (a = 0; a < (2 * shift + 1); a++) {
                corp[a] *= cmax;
            }
        }
        
    } else {
        for (a = 0; a < (2 * shift + 1); a++) {
            corp[a] = 0;
        }
    }
    
}



std::vector<complex> modQAM(std::vector<int> in, int k)
{
	const int I_2[2] = {0, 0};
	const int Q_2[2] = {-1, 1};
	const int I_4[4] = {-1, -1, 1, 1};
	const int Q_4[4] = {-1, 1, -1, 1};
	const int I_16[16] = {-3, -3, -3, -3, -1, -1, -1, -1, 3, 3, 3, 3, 1, 1, 1, 1};
	const int Q_16[16] = {-3, -1, 3, 1, -3, -1, 3, 1, -3, -1, 3, 1, -3, -1, 3, 1};
	const int I_64[64] = {-7, -7, -7, -7, -7, -7, -7, -7, -5, -5, -5, -5, -5, -5, -5,
		-5, -1, -1, -1, -1, -1, -1, -1, -1, -3, -3, -3, -3, -3, -3, -3, -3, 7, 7, 7,
		7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3};
	const int Q_64[64] = {-7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7,
		-5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1,
		3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3, -7, -5, -1, -3, 7, 5, 1, 3};
	std::vector<int> in_dec(int(in.size()/k));
	for (unsigned int i = 0; i < in_dec.size(); i++) {
		in_dec[i] = 0;
		for (int j = k-1; j > -1; j--){
			in_dec[i] = in_dec[i] << 1;
			in_dec[i] += in[k*i+j];
		}
	}
	std::vector<complex> mod(in_dec.size());
	switch (k)
	{
	case 1:
		for (unsigned int i = 0; i < in_dec.size(); i++) {
			mod[i] = (I_2[in_dec[i]] + complex::j* Q_2[in_dec[i]]);
		}
		break;
	case 2:
		for (unsigned int i = 0; i < in_dec.size(); i++) {
			mod[i] = (I_4[in_dec[i]]+ complex::j* Q_4[in_dec[i]]);
            //printf("%d %4.2f + %4.2fj\n",in_dec[i],mod[i].re(),mod[i].im());

		}
		break;
	case 4:
		for (unsigned int i = 0; i < in_dec.size(); i++) {
			mod[i] = (I_16[in_dec[i]]+ complex::j*Q_16[in_dec[i]]);
		}
		break;
	case 6:
		for (unsigned int i = 0; i < in_dec.size(); i++) {
			mod[i] = (I_64[in_dec[i]]+ complex::j* Q_64[in_dec[i]]);
		}
		break;
	default:
		break;
	}
	return mod;
}


 /*
void createPack(double* outPack_i, double* outPack_q,int Nsymb, double* dataIn_i, 
                double* dataIn_q, double* pilot_i, double* pilot_q, double* window, 
                int Msc, int NFFT, int Tpr, int Tpf, double* spec_koef_i, double* spec_koef_q,
                 int flag, double* phi, int Nsymb_aft_pilot)
 {
	int Msc2 = int(Msc/2);
	int complexSz = sizeof(complex);
	complex* pilot = new complex[NFFT];
	complex* spec_koef = new complex[NFFT];
	for (int i = 0; i < NFFT; i++)
		{
			pilot[i]		=  pilot_i[i] + complex::j * pilot_q[i];
			spec_koef[i]	=  spec_koef_i[i] + complex::j * spec_koef_q[i];
		}
	complex* dataIn = new complex[Msc*Nsymb];
	for (int i = 0; i < NFFT; i++)
		{
			dataIn[i]		=  dataIn_i[i] + complex::j * dataIn_q[i];
		}
	complex* Spilot = new complex[NFFT+Tpr+Tpf];
	rfft.Inverse(pilot, Spilot+Tpr, NFFT);
	memcpy(Spilot, Spilot+NFFT, complexSz*Tpr);
	memcpy(Spilot+Tpr+NFFT, Spilot+Tpr, complexSz*Tpf);
	complex* Sint = new complex[NFFT];
	complex* Stmp = new complex[NFFT + Tpr + Tpf];
	complex* Sout = new complex[(NFFT + Tpr + Tpf) * Nsymb];
	for (int i = 0; i< Nsymb; i++)
	{
		memset(Sint, 0, complexSz*NFFT);
		memcpy(Sint, dataIn+i*Msc, complexSz*Msc2);
		memcpy(Sint+NFFT-Msc2, dataIn+i*Msc+Msc2, complexSz*Msc2);
		rfft.Inverse(Sint, Stmp+Tpr, NFFT);
		memcpy(Stmp, Stmp+NFFT, complexSz*Tpr);
		memcpy(Stmp+Tpr+NFFT, Stmp+Tpr, complexSz*Tpf);
		for (int j = 0; j < NFFT+Tpr+Tpf; j++)
		{
			Stmp[j] = Stmp[j] * window[j];
		}
		memcpy(Sout + i * (NFFT + Tpr + Tpf), Stmp, complexSz*(NFFT + Tpr + Tpf));
	}
	complex* outPack = new complex[(NFFT + Tpr + Tpf) * (Nsymb + int(Nsymb/Nsymb_aft_pilot) + 1)];
	switch (flag)
	{
	case 1:
		{
			complex* packPtr = outPack;
			memcpy(outPack, Spilot, complexSz*(NFFT+Tpr+Tpf));
			packPtr += NFFT+Tpr+Tpf;
			for (int i = 0; i < Nsymb/Nsymb_aft_pilot; i++)
			{
				memcpy(packPtr, Spilot, complexSz*(NFFT+Tpr+Tpf));
				packPtr += NFFT+Tpr+Tpf;
				memcpy(packPtr, Sout + i * (NFFT+Tpr+Tpf) * Nsymb_aft_pilot, complexSz*(NFFT+Tpr+Tpf) * Nsymb_aft_pilot);
				packPtr += (NFFT+Tpr+Tpf) * Nsymb_aft_pilot;
			}
			break;
		}
	case 0:
		{
			memcpy(outPack, Spilot, complexSz*(NFFT+Tpr+Tpf));
            break;
		}
	default:
		break;
	}
	for (int i = 0; i < (NFFT + Tpr + Tpf) * (Nsymb + int(Nsymb/Nsymb_aft_pilot) + 1); i++)
		{
			outPack_i[i] =  outPack[i].re();
			outPack_q[i] =  outPack[i].im();
		}
	delete [] outPack;
	delete [] dataIn;
	delete [] pilot;
	delete [] spec_koef;
	delete [] Spilot;
	delete [] Stmp;
	delete [] Sint;
	delete [] Sout;
}
*/


template<typename T> 
std::vector<T> conv(std::vector<T> const &f, std::vector<T> const &g) {
  int const nf = f.size();
  int const ng = g.size();
  int const n  = nf + ng - 1;
  std::vector<T> out(n, T());
  for(auto i(0); i < n; ++i) {
    int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
    int const jmx = (i <  nf - 1)? i            : nf - 1;
    for(auto j(jmn); j <= jmx; ++j) {
      out[i] += (f[j] * g[i - j]);
    }
  }
  return out; 
}

std::vector<complex>  myMod(int M,int Nsymb)
{
    int wordsize = log2(M);   
    vector<int> data(Nsymb*wordsize/4);
    for (int i=0;i<Nsymb*wordsize/4; i++)
    {
        srand(807+i); 
		data[i] = rand()%2;
    }
    std::vector<complex> dataIn(Nsymb/4);
	dataIn = modQAM(data, wordsize);
    
    std::vector<complex> dataOut(Nsymb);
    for (int i =0; i<Nsymb;i++)
    {
        dataOut[i] =(i%4==0)?dataIn[int(i/4)]:0;
    }
    
    /////////////////////////////
    const double dRRC[41] = {-0.00375178919426627,-0.00571593347774432,-0.00146669878754364,0.00640237868930219,0.0106116623233924,0.00496465154503783,-0.00914936129145000,-0.0213492558100365,-0.0187589459713313,0.00301124744697750,0.0326521730597412,0.0470574653755551,0.0265291558084809,-0.0275011161285712,-0.0851595126299198,-0.0993711657502577,-0.0321226123649396,0.118945853798775,0.310937757839833,0.471641896778595,0.534222039323268,0.471641896778595,0.310937757839833,0.118945853798775,-0.0321226123649396,-0.0993711657502577,-0.0851595126299198,-0.0275011161285712,0.0265291558084809,0.0470574653755551,0.0326521730597412,0.00301124744697750,-0.0187589459713313,-0.0213492558100365,-0.00914936129145000,0.00496465154503783,0.0106116623233924,0.00640237868930219,-0.00146669878754364,-0.00571593347774432,-0.00375178919426627};
    std::vector<complex> RRC(41);
    for (int i =0; i<41;i++)
    {
        RRC[i] = dRRC[i];
    }
    std::vector<complex> tmp(Nsymb+40);
    tmp = conv(dataOut,RRC);
    std::vector<complex> res(Nsymb);
    for (int i =0; i<Nsymb;i++)
    {
        res[i] = tmp[i+40];
    } 
    
    return res;
}


void Chirp(int len, double f0, double f1, double Fd, complex *Sout)
{
    int p =1;
    double t1 = len/Fd;
    double phi = 0;
    const double pi = acos(-1);
    double beta  = (f1-f0)*pow(t1,(-p));
    for (int time =0; time<len;time++)
    {
        Sout[time] = (cos(2*pi*phi/360)+complex::j*sin(2*pi*phi/360))*(cos(2*pi*(beta/(1+p)*pow(time,(1+p))+f0*time))-complex::j*cos(2*pi*(beta/(1+p)*pow(time,(1+p))+f0*time+90/360)));
    }
}
