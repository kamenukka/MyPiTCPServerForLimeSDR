#include <stdio.h>

#include <vector>


#include <iostream>

#include <algorithm>

#include <cmath>

#include <math.h>

#include "fft.h"

#include "complex.h"

CFFT rfft;

//const complex complex::i(0., 1.);

//const complex complex::j(0., 1.);









int argminim(const double *arr, size_t length) {

    // returns the minimum value of array

    size_t i;

    double minimum = arr[0];

	int argmin = 0;

    for (i = 1; i < length; ++i) {

        if (minimum > arr[i]) {

            minimum = arr[i];

			argmin = i;

        }

    }

    return argmin;

}







void slice(complex* arr,int size_a, int pa, int pb,complex* b) 

{

	

	for (int i = 0; i<pb-pa; i++)

	{



			b[i] = arr[i+pa];



	}



}





void slice(double* arr,int size_a, int pa, int pb,double* b) 

{



	for (int i = 0; i<pb-pa; i++)

	{



			b[i] = arr[i+pa];



	}



}







std::vector<complex> slice(std::vector<complex> arr, int pa, int pb) 

{

	std::vector<complex> b(pb-pa);

	for (int i = 0; i<pb-pa; i++)

	{



			b[i] = arr[i+pa];



	}



	return b;



}



std::vector<int> slice(std::vector<int> arr, int pa, int pb) 

{

	std::vector<int> b(pb-pa);

	for (int i = 0; i<pb-pa; i++)

	{



			b[i] = arr[i+pa];



	}



	return b;



}



std::vector<double> slice(std::vector<double> arr, int pa, int pb) 

{

	std::vector<double> b(pb-pa);

	for (int i = 0; i<pb-pa; i++)

	{



			b[i] = arr[i+pa];



	}



	return b;



}



//void wconv1(double* x,double* wfilter, int len_x,int len_f, double* y)

//{

//	int n_conv = len_x + len_f-1;

//	int i1;

//	double tmp;

//	for(int i =0 ; i<n_conv; i++)

//	{

//		i1 =i;

//		tmp = 0.0;

//		for (int j =0; j<len_f; j++)

//		{

//			if (i1>=0 && i1<len_x)

//			{

//				tmp = tmp + x[i1]*wfilter[j];

//			}

//

//			i1 = i1-1;

//			y[i] = tmp;

//		}

//	}

//

//}



std::vector<double> wconv1(std::vector<double> x,std::vector<double> wfilter, int len_x,int len_f)

{

	int n_conv = len_x + len_f-1;

	int i1;

	double tmp;

	std::vector<double> y(n_conv);

	for(int i =0 ; i<n_conv; i++)

	{

		i1 =i;

		tmp = 0.0;

		for (int j =0; j<len_f; j++)

		{

			if (i1>=0 && i1<len_x)

			{

				tmp = tmp + x[i1]*wfilter[j];

			}



			i1 = i1-1;

			y[i] = tmp;

		}

	}

	return y;

}



//

//

//void wconv2(double* x,double* wfilter, int len_x,int len_f, double* y)

//{

//	int n_conv = len_x + len_f-1;

//	int i1;

//	double tmp;

//	double * zz = new double[n_conv];

//	for(int i =0 ; i<len_x; i++)

//	{

//		i1 =i;

//		tmp = 0.0;

//		for (int j =0; j<len_f; j++)

//		{

//			if (i1>=0 && i1<len_x)

//			{

//				tmp = tmp + x[i1]*wfilter[j];

//			}

//

//			i1 = i1-1;

//			zz[i] = tmp;

//		}

//	}

//

//	for (int i =0;i<(len_x+len_f-1); i++)

//	{

//		y[i] = *(zz +len_f-1+i); 

//		//y[i] = zz[i];

//	}

//	//delete [] zz;

//}

//



std::vector<double> fliplr(std::vector<double> x,int len_x)

{

	std::vector<double> y(len_x);

	for(int i=0; i<len_x; i++)

	{

		y[len_x-1-i] = x[i];

	}

	return y;

}





void fliplr(double* x,int len_x)

{

	double * tmp = new double[len_x];



	for(int i=0; i<len_x; i++)

	{

		tmp[len_x-1-i] = x[i];

	}

	for(int i=0; i<len_x; i++)

	{

		x[i] = tmp[i];

	}

	delete [] tmp;

}



void fliplr(complex* x,int len_x)

{

	complex * tmp = new complex[len_x];



	for(int i=0; i<len_x; i++)

	{

		tmp[len_x-1-i] = x[i];

	}

	for(int i=0; i<len_x; i++)

	{

		x[i] = tmp[i].conjugate();

	}

	delete [] tmp;

}



std::vector<int> fliplr(std::vector<int> x,int len_x)

{

	std::vector<int> y(len_x);

	for(int i=0; i<len_x; i++)

	{

		y[len_x-1-i] = x[i];

	}

	return y;

}



std::vector<double> wextend(std::vector<double> x,int len_x, int lenEXT)

{   // specialno ne budu dopisivat' cod dlya lenEXT>len_x, ibo pri divx64 32 podnesushih

	//I = [lf:-1:1 , 1:lx , lx:-1:lx-lf+1];

	std::vector<double> y(2*lenEXT+len_x);

	std::vector<int> I(2*lenEXT+len_x);

	

	//	for (int j=0; j<(len_x);j++)

	//{

	//	std::cout<<"x["<<j+1<<"] = "<<x[j]<<std::endl;



	//}

	//	std::cout<<"____________"<<std::endl;

	for (int j=0;j<lenEXT;j++)

	{

		I[j] = lenEXT-1-j;

		I[2*lenEXT+len_x-1-j] = len_x - lenEXT + j ;

	}



	for (int j=0;j<len_x;j++)

	{

		I[lenEXT+j] = j;

	}



	for (int j=0; j<(2*lenEXT+len_x);j++)

	{

		y[j]=x[I[j]];

		//std::cout<<"I["<<j+1<<"] = "<<I[j]<<std::endl;



	}

	//delete [] I;

	return y;

}









//void wavedec(double* x,int n,int name_db, double* c,double* l)

void dwt(std::vector<double> x,int len_x,int name_db, std::vector<double>  * a,std::vector<double>  * d)

{



	double Lo_D_[16] = {-0.000117476784002282,	0.000675449405998557,	-0.000391740372995977,	-0.00487035299301066,

		0.00874609404701566,	0.0139810279170155,	-0.0440882539310647,	-0.0173693010020221,	0.128747426620186,

		0.000472484573997973,	-0.284015542962428,	-0.0158291052560239,	0.585354683654869,	0.675630736298013,

		0.312871590914466,	0.0544158422430816};

	double Hi_D_[16] = {-0.0544158422430816,	0.312871590914466,	-0.675630736298013,	0.585354683654869,	0.0158291052560239,

		-0.284015542962428,	-0.000472484573997973,	0.128747426620186,	0.0173693010020221,	-0.0440882539310647,

		-0.0139810279170155,	0.00874609404701566,	0.00487035299301066,	-0.000391740372995977,	-0.000675449405998557,

		-0.000117476784002282};

	std::vector<double> Lo_D(16);

	std::vector<double> Hi_D(16);

	for (int i=0;i<16;i++)

	{

		Lo_D[i] = Lo_D_[i];

		Hi_D[i] = Hi_D_[i]; 

	}

	

	int first = 1;

	int lenEXT = 16-1;

	int last = len_x+16-1;

	int leny = 2*lenEXT + len_x;

	std::vector<double> y(leny);

	y = wextend(x,len_x,lenEXT);



	//for (int i=0;i<leny;i++)

	//{

	//	std::cout<<"y["<<i+1<<"] = "<<y[i]<<std::endl;

	//}



	int lenz = leny+16+1;

	std::vector<double> z(lenz);

	z = wconv1(y,Lo_D,leny,16);

	int lend = int((last - first+1+1)/2);

	std::vector<double> a_tmp(lend);

	std::vector<double> d_tmp(lend);

	int k = first;

	//std::cout<<"__________"<<std::endl;

	for (int j = 0;j<lend;j++)

	{

		a_tmp[j] = z[k+15];

		//a->push_back(z[k+15]);

		k = k+2;

	}



	z = wconv1(y,Hi_D,leny,16);

	k = first;

	//std::cout<<"__________"<<std::endl;

	for (int j = 0;j<lend;j++)

	{

		d_tmp[j] = z[k+15];

		//d->push_back(z[k+15]);

		//std::cout<<"d["<<j+1<<"] = "<<d[j]<<std::endl;

		k = k+2;

	}

	a->swap(a_tmp);

	d->swap(d_tmp);





	//for (int i = 0;i<lend;i++)

	//{

	//	std::cout<<a[i]<<"   "<<d[i]<<" "<<std::endl;

	//}

	/*delete [] x_;

	delete [] y;

	delete [] z;*/

}



std::vector<double> dyadup(std::vector<double> x, int len_x)

{

	int leny =len_x*2-1;

	std::vector<double> y(leny);

	int k =0;

	for (int i =0;i<len_x-1;i++)

	{

		y[k] = x[i];

		y[k+1] =0;

		k +=2;

	}

	y[leny-1] = x[len_x-1];

	return y;

}









//double * wkeep1(double *x, int len_x,int len)

//{

//	double d = (len_x -len)/2;

//	int first = int(d);

//	int kostyl = 0;

//	int kostyl2 = 0;

//	if ((len&1)^(len_x&1))  kostyl = 1;

//	int last = len_x - kostyl -int(d);

//	double * y  = new double[len];

//	slice(x,len_x,first,last,y);

//	return y;

//}



std::vector<double> wkeep1( std::vector<double> x, int len_x,int len)

{

	double d = (len_x -len)/2;

	

	int first = int(d);

	int kostyl = 0;

	int kostyl2 = 0;

	if ((len&1)^(len_x&1))  kostyl = 1;

	int last = len_x - kostyl -int(d);

	std::vector<double> y(last-first);

	y = slice(x,first,last);



	return y;



}



std::vector<double> upsconv1(std::vector<double> x, std::vector<double> f, int len_x,int len_f, int s)

{

	int lx = 2*len_x;

	int lf = len_f;

	std::vector<double> d(lx-1);

	//d = dyadup(x,len_x);

	d = dyadup(x,len_x);



	std::vector<double> z(lx-1+len_f-1);

	

	z = wconv1(d,f,2*len_x-1,len_f);

	std::vector<double> y(s);

	y = wkeep1(z,(lx-1+len_f-1),s);

	

	return y;

}



double cumsum1(std::vector<double> x, int len_x,int num)

{

	double res = 0;

	for (int i=0;i<num;i++)

	{

		res += x[i];

	}

	return res;

}



int cumsum1(std::vector<int> x, int len_x,int num)

{

	int res = 0;

	for (int i=0;i<num;i++)

	{

		res += x[i];

	}

	return res;

}



std::vector<double> cumsum(std::vector<double> x, int len_x)

{

	//double * res = new double[len_x];

	std::vector<double> res(len_x);

	for (int i=0;i<len_x;i++)

	{

		res[i] = cumsum1(x,len_x,i+1);

	}

	return res;

}



std::vector<int> cumsum(std::vector<int> x, int len_x)

{

	//double * res = new double[len_x];

	std::vector<int> res(len_x);

	for (int i=0;i<len_x;i++)

	{

		res[i] = cumsum1(x,len_x,i+1);

	}

	return res;

}



std::vector<double> detcoef(std::vector<double> coefs, std::vector<int> longs, int len_coef,int len_longs, int levels,int * len_tmp_ )

{

	int nmax = len_longs-2;

	std::vector<int> first_tmp(len_longs-2);







	first_tmp = cumsum(longs,len_longs-2);

	

	std::vector<int> first(len_longs-2);

	first = fliplr(first_tmp,len_longs-2);



	for (int i = 0; i<len_longs-2;i++)

	{

		first[i] +=1;

	}

	std::vector<int> longs_tmp(len_longs-2);

	std::vector<int> longs_(len_longs-2);



	longs_tmp =  slice(longs,1,len_longs-1);

	longs_ = fliplr(longs_tmp,len_longs-2);



	std::vector<int> last(len_longs-2);

	for (int i = 0;i<len_longs-2;i++)

	{

		last[i] = first[i] + longs[len_longs-2-i]-1;



			//std::cout<<first[i] <<"  "<< last[i]<<std::endl;



	}

	int len_tmp = int(last[levels-1]-first[levels-1]);









	std::vector<double> tmp(len_tmp);

	tmp = slice(coefs,int(first[levels-1])-1,int(last[levels-1]));

	len_tmp_[0] = len_tmp;



	return tmp;



	

}







std::vector<double> idwt(std::vector<double> a,std::vector<double> d, int len_a,int name_db, int l)

{



	double Lo_D_[16] = {0.0544158422430816, 0.312871590914466, 0.675630736298013, 0.585354683654869, -0.0158291052560239, -0.284015542962428, 0.000472484573997973, 0.128747426620186, -0.0173693010020221, -0.0440882539310647, 0.0139810279170155, 0.00874609404701566, -0.00487035299301066, -0.000391740372995977, 0.000675449405998557, -0.000117476784002282};

	double Hi_D_[16] ={-0.000117476784002282, -0.000675449405998557, -0.000391740372995977, 0.00487035299301066, 0.00874609404701566, -0.0139810279170155, -0.0440882539310647, 0.0173693010020221, 0.128747426620186, -0.000472484573997973, -0.284015542962428, 0.0158291052560239, 0.585354683654869, -0.675630736298013, 0.312871590914466, -0.0544158422430816};

	std::vector<double> Lo_D(16);

	std::vector<double> Hi_D(16);

	for (int i=0;i<16;i++)

	{

		Lo_D[i] = Lo_D_[i];

		Hi_D[i] = Hi_D_[i]; 

	}





	std::vector<double> x_Lo(l);

	std::vector<double> x_Hi(l);

	std::vector<double> x(l);

	x_Lo = upsconv1(a,Lo_D,len_a,16,l);

	x_Hi = upsconv1(d,Hi_D,len_a,16,l);

	for (int i=0; i<l;i++)

	{

		x[i] = x_Lo[i] + x_Hi[i];

	}



	return x;



}





std::vector<double> appcoef(std::vector<double> c, std::vector<int> l, int lenc, int lenl, int * lena)

{

	int rmax = lenl;

	int nmax = rmax-2;

	std::vector<double> a_(int(l[0]*lenl));

	a_ = slice(c,0,int(l[0]));

	int imax = rmax +1;

	std::vector<double> d(lenc);

	int * lend = new int[1];

	int lena_ = int(l[0]);

	int lena_2 = int(l[0]);

	std::vector<double> a_tmp(int(l[0]*lenl));





	for (int p = nmax;p>0;p--)

	{

		d = detcoef(c,l,lenc,lenl,p,lend);





		a_tmp = idwt(a_,d,lena_,8,l[imax-p-1]);

		lena_ = l[imax-p-1];

		

		for (int i=0;i<lena_2;i++)

		{

			a_[i] = a_tmp[i];

		}

		for (int i=0;i<lena_-lena_2;i++)

		{

			a_.push_back(a_tmp[i+lena_2]);

		}

		lena_2 = lena_;

	}

	std::vector<double> a(lena_);

		for (int i=0;i<lena_;i++)

	{

		a[i] = a_[i];

	}

	lena[0] = lena_;

	return a;

}



std::vector<double> wthresh(std::vector<double> x, int len_x, double t)

{

	std::vector<double> y(len_x);

	double tmp;

	int mn = 1;

	for (int i=0;i<len_x;i++)

	{

	mn = 1;

	tmp = (abs(x[i])-t + abs(abs(x[i])-t))/2;



	if (x[i]<0) 

	{

		mn = -1;

	}



	y[i] = mn * tmp;



	}

	return y;



}



//int cumsum1(int* x, int len_x,int num)

//{

//	int res = 0;

//	for (int i=0;i<num;i++)

//	{

//		res += x[i];

//	}

//	return res;

//}



std::vector<double> hstack(std::vector<double> x,std::vector<double> y, int lenx, int leny)

{

	int len =lenx + leny;

	std::vector<double> z(len);

	for (int i=0;i<lenx;i++)

	{

		z[i] = x[i];

	}

	for (int i=0;i<leny;i++)

	{

		z[i+lenx] = y[i];

	}

	return z;

}







void wavedec(std::vector<double> x,int len_x,int n, std::vector<double> *c,std::vector<int> *l, int* lenc)

{

	std::vector<double> x_(len_x);

	x_.swap(x);



	std::vector<int> l_tmp(n+2);

	for (int i=0;i<n+2;i++)

	{

		l_tmp[i] = 0;

	}

	l_tmp[n+1]=len_x;

	int len_x_ = len_x;

	int lena;



	int tmp = 0;



	len_x_ = len_x;

	//lenc_ = cumsum1(lena_arr,n,n);

	

	std::vector<double>  c_(len_x*(n+2));

	lena = (len_x_ +14)/2+1;

	std::vector<double> a(lena);

	std::vector<double> d(lena);

	int lenc__ = 0;

	for (int k=0;k<n;k++)

	{

		lena = int((len_x_ +14+1)/2);





		dwt(x_,len_x_,8,&a,&d);



		//std::cout<<"________"<<std::endl;

		for (int i=0;i<lena;i++)

		{

			x_[i]=a[i];

			//std::cout<<"x["<<i+1<<"] =   "<<x_[i]<<std::endl;



		}

		l_tmp[n-k]  = lena;

	

		c_ = hstack(d,c_,lena,lenc__);

		//hstack(d,c_,lena,lenc__, c_);



		lenc__ += lena;

		len_x_ = lena;



		

	}



	c_ = hstack(x_,c_,len_x_,lenc__);



	lenc[0] = len_x_+lenc__;





	l_tmp[0] =int(lena);



	l->swap(l_tmp);

	c->swap(c_);







}





std::vector<double> wden(std::vector<double> x,int len_x, int n)

{

	std::vector<double> c(len_x*n+1);

	std::vector<int> l(n+2);

	int * lenc = new int[1];



	wavedec(x,len_x,n,&c,&l,lenc);



	std::vector<int> s(n);

	for (int i=0;i<n;i++)

	{

		s[i] =1;

	}



	std::vector<int> first_tmp(n);

	first_tmp = cumsum(l,n);

	std::vector<int> first(n);

	first = fliplr(first_tmp,n);



	for (int i = 0; i<n;i++)

	{

		first[i] +=1;

	}



	std::vector<int> ld(n);

	std::vector<int> ld_(n);

	ld_ = slice(l,1,n+1);

	ld = fliplr(ld_,n);





	std::vector<int> last(n);

	for (int i = 0; i<n;i++)

	{

		last[i] = first[i] + ld[i] -1;

		//std::cout<<"first = "<<first[i]<<"     last = "<<last[i]<<std::endl;

		

	}

	double thr = sqrt(2*log(lenc[0]));

	int flk1 =0;

	int flk2 = 0;

	std::vector<double> cxd(lenc[0]);



	for (int i=0; i<lenc[0];i++)

	{

		cxd[i] = c[i];

	}

	std::vector<double> cxd_tmp(lenc[0]);

	std::vector<double> cxd_tmp2(lenc[0]);

	for (int k =0;k<n;k++)

	{

		flk1 = int(first[k]);

		flk2 = int(last[k]);

		thr = sqrt(2*log(lenc[0])) * s[k];

		cxd_tmp = slice(c,flk1-1,flk2);

		cxd_tmp2 = wthresh(cxd_tmp,flk2-flk1+1,thr);

		for (int i = 0;i<flk2-flk1+1;i++)

		{

			cxd[flk1-1+i] = cxd_tmp2[i];

			

		}

	}



	//std::cout<<"do appcoef rabotaet"<<std::endl;

	std::vector<double> xd(len_x);

	int * lenx = new int[1]; 

	xd = appcoef(cxd,l,lenc[0],n+2,lenx);

	//std::cout<<"appcoef rabotaet"<<std::endl;





	return xd;

}



void kostyleviywden(double * x, int len,int n, double * y)

{

	std::vector<double> xx(len);

	std::vector<double> yy(len);



	for (int i=0; i<len;i++)

	{

		xx[i] = x[i];

	}

	yy = wden(xx,len,n);

	for (int i=0; i<len;i++)

	{

		y[i] = yy[i];

	}



}





/* */double minim(const double *arr, size_t length) {

    // returns the minimum value of array

    size_t i;

    double minimum = arr[0];

    for (i = 1; i < length; ++i) {

        if (minimum > arr[i]) {

            minimum = arr[i];

        }

    }

    return minimum;

}



 void soft_detector(double *mod_rx_i, double *mod_rx_q, int len, int wordsize, double *bit_metric) 

{ 

double* Im0;

double* Qm0;

double* Im1;

double* Qm1;



if (wordsize == 2)

{

	for (int ji=0; ji<len; ji++)

	{

		mod_rx_i[ji] = mod_rx_i[ji] /sqrt(2);

		mod_rx_q[ji] = mod_rx_q[ji] /sqrt(2);

	}

}

else

	{

		if (wordsize == 4)

		{

			for (int ji=0; ji<len; ji++)

			{

				mod_rx_i[ji] = mod_rx_i[ji] /sqrt(10);

				mod_rx_q[ji] = mod_rx_q[ji] /sqrt(10);

			}

		}

		else

		{

			if (wordsize == 6)

			{

				for (int ji=0; ji<len; ji++)

				{

					mod_rx_i[ji] = mod_rx_i[ji] /sqrt(42);

					mod_rx_q[ji] = mod_rx_q[ji] /sqrt(42);

				}

			}

		}

}



double Im0_1[2][2] =  {{0, 0}, {0, 0}};

double Qm0_1[2][2] = {{-1, - 1}, {-1, - 1}};

double Im1_1[2][2] = {{0, 0}, {0, 0}};

double Qm1_1[2][2] = {{1, 1}, {1, 1}};





double Im0_2[2][2] =  {{-0.70711, 0.70711},

                            {-0.70711, - 0.70711}};

double Qm0_2[2][2] = {{-0.70711, - 0.70711},

                            {-0.70711, 0.70711}};

double Im1_2[2][2] ={{-0.70711, 0.70711},

                            {0.70711, 0.70711}};

double Qm1_2[2][2] = {{0.70711, 0.70711},

                            {-0.70711, 0.70711}};



//-0.94868



double Im0_4[4][8] = {{-0.94868, 0.94868, - 0.31623, 0.31623, - 0.94868, 0.94868, - 0.31623, 0.31623},

                            {-0.94868, 0.94868, - 0.31623, 0.31623, - 0.94868, 0.94868, - 0.31623, 0.31623},

                            {-0.94868, 0.94868, - 0.94868, 0.94868, - 0.94868, 0.94868, - 0.94868, 0.94868},

                            {-0.94868, - 0.31623, - 0.94868, - 0.31623, - 0.94868, - 0.31623, - 0.94868, - 0.31623}};

double Qm0_4[4][8] = {{-0.94868, - 0.94868, - 0.94868, - 0.94868, 0.94868, 0.94868, 0.94868, 0.94868},

                            {-0.94868, - 0.94868, - 0.94868, - 0.94868, - 0.31623, - 0.31623, - 0.31623, - 0.31623},

                            {-0.94868, - 0.94868, 0.94868, 0.94868, - 0.31623, - 0.31623, 0.31623, 0.31623},

                            {-0.94868, - 0.94868, 0.94868, 0.94868, - 0.31623, - 0.31623, 0.31623, 0.31623}};

double Im1_4[4][8] = {{-0.94868, 0.94868, - 0.31623, 0.31623, - 0.94868, 0.94868, - 0.31623, 0.31623},

                            {-0.94868, 0.94868, - 0.31623, 0.31623, - 0.94868, 0.94868, - 0.31623, 0.31623},

                            {-0.31623, 0.31623, - 0.31623, 0.31623, - 0.31623, 0.31623, - 0.31623, 0.31623},

                            {0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623, 0.94868, 0.31623}};

double Qm1_4[4][8] = {{-0.31623, - 0.31623, - 0.31623, - 0.31623, 0.31623, 0.31623, 0.31623, 0.31623},

                            {0.94868, 0.94868, 0.94868, 0.94868, 0.31623, 0.31623, 0.31623, 0.31623},

                            {-0.94868, - 0.94868, 0.94868, 0.94868, - 0.31623, - 0.31623, 0.31623, 0.31623},

                            {-0.94868, - 0.94868, 0.94868, 0.94868, - 0.31623, - 0.31623, 0.31623, 0.31623}};







double Im0_6[6][32] = {{-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801,

                             -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543,

                             -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291},



                            {-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801,

                             -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543,

                             -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291},



                            {-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801,

                             -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543,

                             -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291},



                            {-1.0801, 1.0801, -0.1543, 0.1543, -1.0801, 1.0801, -0.1543, 0.1543, -1.0801, 1.0801,

                             -0.1543, 0.1543, -1.0801, 1.0801, -0.1543, 0.1543, -1.0801, 1.0801, -0.1543, 0.1543,

                             -1.0801, 1.0801, -0.1543, 0.1543, -1.0801, 1.0801, -0.1543, 0.1543, -1.0801, 1.0801,

                             -0.1543, 0.1543},



                            {-1.0801, 1.0801, -0.77152, 0.77152, -1.0801, 1.0801, -0.77152, 0.77152, -1.0801, 1.0801,

                             -0.77152, 0.77152, -1.0801, 1.0801, -0.77152, 0.77152, -1.0801, 1.0801, -0.77152, 0.77152,

                             -1.0801, 1.0801, -0.77152, 0.77152, -1.0801, 1.0801, -0.77152, 0.77152, -1.0801, 1.0801,

                             -0.77152, 0.77152},



                            {-1.0801, -0.1543, -0.77152, -0.46291, -1.0801, -0.1543, -0.77152, -0.46291, -1.0801,

                             -0.1543, -0.77152, -0.46291, -1.0801, -0.1543, -0.77152, -0.46291, -1.0801, -0.1543,

                             -0.77152, -0.46291, -1.0801, -0.1543, -0.77152, -0.46291, -1.0801, -0.1543, -0.77152,

                             -0.46291, -1.0801, -0.1543, -0.77152, -0.46291}

                            };

double Qm0_6[6][32] = {

                            {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801,

                            1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543,

                            0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543},



                            {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801,

                             1.0801, 1.0801, 1.0801, 1.0801, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152,

                             -0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152},



                            {-1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -1.0801, -0.1543, -0.1543, -0.1543,

                             -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152,

                             -0.77152, -0.77152, -0.77152, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291,

                             -0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291}

                            };

double Im1_6[6][32] = {

                            {-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543,

                             0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291},



                            {-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543,

                             0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291},



                            {-1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543,

                             0.1543, -0.77152, 0.77152, -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152,

                             -0.46291, 0.46291, -1.0801, 1.0801, -0.1543, 0.1543, -0.77152, 0.77152, -0.46291, 0.46291},



                            {-0.77152, 0.77152, -0.46291, 0.46291, -0.77152, 0.77152, -0.46291, 0.46291, -0.77152, 0.77152,

                             -0.46291, 0.46291, -0.77152, 0.77152, -0.46291, 0.46291, -0.77152, 0.77152, -0.46291, 0.46291,

                             -0.77152, 0.77152, -0.46291, 0.46291, -0.77152, 0.77152, -0.46291, 0.46291, -0.77152, 0.77152,

                             -0.46291, 0.46291},



                            {-0.1543, 0.1543, -0.46291, 0.46291, -0.1543, 0.1543, -0.46291, 0.46291, -0.1543, 0.1543, -0.46291,

                             0.46291, -0.1543, 0.1543, -0.46291, 0.46291, -0.1543, 0.1543, -0.46291, 0.46291, -0.1543, 0.1543,

                             -0.46291, 0.46291, -0.1543, 0.1543, -0.46291, 0.46291, -0.1543, 0.1543, -0.46291, 0.46291},



                            {1.0801, 0.1543, 0.77152, 0.46291, 1.0801, 0.1543, 0.77152, 0.46291, 1.0801, 0.1543, 0.77152, 0.46291,

                             1.0801, 0.1543, 0.77152, 0.46291, 1.0801, 0.1543, 0.77152, 0.46291, 1.0801, 0.1543, 0.77152, 0.46291,

                             1.0801, 0.1543, 0.77152, 0.46291, 1.0801, 0.1543, 0.77152, 0.46291}

                            };

double Qm1_6[6][32] = {

                            {-0.77152, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152,

                             0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291,

                             -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, -0.1543, 0.1543, 0.1543, 0.1543, 0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291, -0.46291,

                             -0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 1.0801, 0.1543, 0.1543, 0.1543, 0.1543, 0.1543,

                             0.1543, 0.1543, 0.1543, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152, 0.77152,

                             0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291},



                            {-1.0801, -1.0801, -1.0801, -1.0801, 1.0801, 1.0801, 1.0801, 1.0801, -0.1543, -0.1543, -0.1543, -0.1543,

                             0.1543, 0.1543, 0.1543, 0.1543, -0.77152, -0.77152, -0.77152, -0.77152, 0.77152, 0.77152, 0.77152,

                             0.77152, -0.46291, -0.46291, -0.46291, -0.46291, 0.46291, 0.46291, 0.46291, 0.46291}

                            };









int count;

if (wordsize == 1)

{

Im0 = &Im0_1[0][0];

Qm0 = &Qm0_1[0][0];

Im1 = &Im1_1[0][0];

Qm1 = &Qm1_1[0][0];

count = 2;



}

else

{

	if (wordsize == 2)

{

Im0 = &Im0_2[0][0];

Qm0 = &Qm0_2[0][0];

Im1 = &Im1_2[0][0];

Qm1 = &Qm1_2[0][0];

count = 2;

}

else

{

if (wordsize == 4)

{

Im0 = &Im0_4[0][0];

Qm0 = &Qm0_4[0][0];

Im1 = &Im1_4[0][0];

Qm1 = &Qm1_4[0][0];

count = 8;

}

else

{

if (wordsize == 6)

{

Im0 = &Im0_6[0][0];

Qm0 = &Qm0_6[0][0];

Im1 = &Im1_6[0][0];

Qm1 = &Qm1_6[0][0];

count = 32;

}

}

	}

}





double* eq1 = new double[count];

double* eq2 = new double[count];



double m0;

int ij=0;



        for (int sds=0; sds <len; sds++){

            for (int sdb = 0 ; sdb< wordsize; sdb ++ ){

				for (int k = 0; k< count; k ++ ){

					*(eq1+k) = abs(mod_rx_i[sds] - *(Im0+sdb*count+k)) + abs(mod_rx_q[sds] - *(Qm0+sdb*count+k));

					*(eq2+k) = abs(mod_rx_i[sds] - *(Im1+sdb*count+k)) + abs(mod_rx_q[sds] - *(Qm1+sdb*count+k));

				}

				m0 = minim(eq1,count) - minim(eq2,count);

				//std::cout << min(eq1,count)<<"          "<<min(eq2,count)<<std::endl;



                *(bit_metric+ij) = m0;

				ij += 1;

		}

		}



		delete [] eq1;

		delete [] eq2;



		}







/* */double maxim(const double *arr, size_t length) {

    // returns the minimum value of array

    size_t i;

	

    double maximum = arr[0];

    for (i = 1; i < length; ++i) {

	

        if (maximum < arr[i]) {

            maximum = arr[i];

        }

    }

    return maximum;

}



/* */double maxi(double a, double b) 

{   double maximum;



	if (a>b)

	{

		maximum =a;

	}

	else

	{

		maximum = b;

	}



    return maximum;

}





/**/ void logmapo(double* sys, double* parity, double* L_u, int L_total, int terminate, int K, double* L_all){



        double Inftymat = -10000;

		double* m11 = new double[L_total];

		double* m10 = new double[L_total];

		int cntr = 16;





		for(int j = 0; j < L_total; j++)

			{

				m11[j] = (sys[j] + L_u[j] + parity[j]) / 2;

				m10[j] = (sys[j] + L_u[j] - parity[j]) / 2;

				//std::cout<<"m11[j] = "<<m11[j]<<",  m10[j] = "<<m10[j]<<std::endl;

			}





		int L_2  = int(pow(2,(K-1)));

		double **Alpha = new double*[L_total+1] ;

		double **Beta = new double*[L_total+1] ;

		for (int i=0; i<L_total+1;i++)

			{

				Alpha[i] = new double[L_2];

				Beta[i] = new double[L_2];

			}





		for(int i = 0; i < L_total+1; i++)

			{

			for(int j = 0; j < L_2; j++)

				{

				if (i==0){

					Alpha[i][ j] = Inftymat;

					}

				else

					{

					Alpha[i][ j] = 0;

				}

				Beta[i][ j] = 0;

				}

			} 

		Alpha[0][0] = 0;

 



        if (terminate == 1)

			{

				Beta[L_total][0] = 0;

				for(int j = 1; j < L_2; j++)

				{

					Beta [L_total][ j] = Inftymat;

				}

			}

		int counter=0;



		double** enumerator = new double*[cntr];

		for (int i=0; i<cntr;i++)

			{

				enumerator[i] = new double[L_total];

			}



		double** denominator = new double*[cntr];

		for (int i=0; i<cntr;i++)

			{

				denominator[i] = new double[L_total];

			}



		for(int i = 0; i < cntr; i++)

			{

			for(int j = 0; j < L_total; j++)

				{

					enumerator[i][ j]  = 0;

					denominator[i][ j] = 0;

				}

			}





        if (K == 3)

			{

			counter = 4;

            for (int k=0; k<L_total; k++)

			{

                Alpha[k + 1][ 0] = maxi(Alpha[k][ 0] - m11[k], Alpha[k][ 2] + m11[k]);

                Alpha[k + 1][ 1] = maxi(Alpha[k][ 0] + m11[k], Alpha[k][ 2] - m11[k]);

                Alpha[k + 1][ 2] = maxi(Alpha[k][ 1] + m10[k], Alpha[k][ 3] - m10[k]);

                Alpha[k + 1][3] = maxi(Alpha[k][ 1] - m10[k], Alpha[k][3] + m10[k]);

			}



             // âû÷èñëåíèå ìåòðèê ñîñòî¤íèé Beta



			for (int k=L_total-1; k > -1 ; k--)

			{

				Beta[k][0] = maxi(Beta[k + 1][ 0] - m11[k], Beta[k + 1][ 1] + m11[k]);

                Beta[k][1] = maxi(Beta[k + 1][ 2] + m10[k], Beta[k + 1][ 3] - m10[k]);

                Beta[k][2] = maxi(Beta[k + 1][ 0] + m11[k], Beta[k + 1][ 1] - m11[k]);

                Beta[k][3] = maxi(Beta[k + 1][ 2] - m10[k], Beta[k + 1][ 3] + m10[k]);





			}

			//std::cout<<"  alpha"<<std::endl;

			//			for(int i = 0; i < L_total+1; i++)

			//		{

			//			for(int j = 0; j < L_2; j++)

			//			{



			//				std::cout<<Alpha[i][ j]<<"  ";



			//		}

			//			std::cout<<std::endl;

			//	} 

			//			std::cout<<"  beta"<<std::endl;



			//									for(int i = 0; i < L_total+1; i++)

			//		{

			//			for(int j = 0; j < L_2; j++)

			//			{



			//				std::cout<<Beta[i][ j]<<"  ";



			//		}

			//			std::cout<<std::endl;

			//	} 







			for (int k=0; k<L_total;k++)

			{

            enumerator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 1] + m11[k];

            enumerator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 2] + m10[k];

            enumerator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 0] + m11[k];

            enumerator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 3] + m10[k];



            denominator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 0] - m11[k];

            denominator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 3] - m10[k];

            denominator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 1] - m11[k];

            denominator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 2] - m10[k];

			}







	//		std::cout<<"  enum"<<std::endl;

	//		for(int i = 0; i < 4; i++)

	//	{

	//		for(int j = 0; j < L_total; j++)

	//		{



	//			std::cout<<enumerator[i][ j]<<"  ";



	//	}

	//		std::cout<<std::endl;

	//} 

	//					std::cout<<"  denom"<<std::endl;

	//		for(int i = 0; i < 4; i++)

	//	{

	//		for(int j = 0; j < L_total; j++)

	//		{



	//			std::cout<<denominator[i][ j]<<"  ";



	//	}

	//		std::cout<<std::endl;

	//} 



		}

        else

		{



            if (K == 4)

			{ counter =8;

                // âû÷èñëåíèå ìåòðèê ñîñòî¤íèé Alpha

                 for (int k=0; k<L_total; k++)

				{

                    Alpha[k + 1][ 0] = maxi(Alpha[k][ 0] - m11[k], Alpha[k][ 1] + m11[k]);

                    Alpha[k + 1][ 1] = maxi(Alpha[k][ 2] - m10[k], Alpha[k][ 3] + m10[k]);

                    Alpha[k + 1][ 2] = maxi(Alpha[k][ 4] + m10[k], Alpha[k][ 5] - m10[k]);

                    Alpha[k + 1][ 3] = maxi(Alpha[k][ 6] + m11[k], Alpha[k][ 7] - m11[k]);

                    Alpha[k + 1][ 4] = maxi(Alpha[k][ 0] + m11[k], Alpha[k][ 1] - m11[k]);

                    Alpha[k + 1][ 5] = maxi(Alpha[k][ 2] + m10[k], Alpha[k][ 3] - m10[k]);

                    Alpha[k + 1][ 6] = maxi(Alpha[k][ 4] - m10[k], Alpha[k][ 5] + m10[k]);

                    Alpha[k + 1][ 7] = maxi(Alpha[k][ 6] - m11[k], Alpha[k][ 7] + m11[k]);



					//std::cout<<"Alpha["<<k + 1<<"][ 0]"<<Alpha[k + 1][ 0]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 1]"<<Alpha[k + 1][ 1]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 2]"<<Alpha[k + 1][ 2]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 3]"<<Alpha[k + 1][ 3]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 4]"<<Alpha[k + 1][ 4]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 5]"<<Alpha[k + 1][ 5]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 6]"<<Alpha[k + 1][ 6]<<std::endl;

					//std::cout<<"Alpha["<<k + 1<<"][ 7]"<<Alpha[k + 1][ 7]<<std::endl;



				 }

                // âû÷èñëåíèå ìåòðèê ñîñòî¤íèé Beta

				for (int k=L_total-1; k > -1 ; k--)

				{

					Beta[k][ 0] = maxi(Beta[k + 1][ 0] - m11[k], Beta[k + 1][ 4] + m11[k]);

                    Beta[k][ 1] = maxi(Beta[k + 1][ 0] + m11[k], Beta[k + 1][ 4] - m11[k]);

                    Beta[k][ 2] = maxi(Beta[k + 1][ 1] - m10[k], Beta[k + 1][ 5] + m10[k]);

                    Beta[k][ 3] = maxi(Beta[k + 1][ 1] + m10[k], Beta[k + 1][ 5] - m10[k]);

                    Beta[k][ 4] = maxi(Beta[k + 1][ 2] + m10[k], Beta[k + 1][ 6] - m10[k]);

                    Beta[k][ 5] = maxi(Beta[k + 1][ 2] - m10[k], Beta[k + 1][ 6] + m10[k]);

                    Beta[k][ 6] = maxi(Beta[k + 1][ 3] + m11[k], Beta[k + 1][ 7] - m11[k]);

                    Beta[k][ 7] = maxi(Beta[k + 1][ 3] - m11[k], Beta[k + 1][ 7] + m11[k]);

				}















				for (int k=0; k<L_total;k++)

				{

					enumerator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 4] + m11[k];

					enumerator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 0] + m11[k];

					enumerator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 5] + m10[k];

					enumerator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 1] + m10[k];

					enumerator[4][k] = Alpha[k][ 4] + Beta[k + 1][ 2] + m10[k];

					enumerator[5][k] = Alpha[k][ 5] + Beta[k + 1][ 6] + m10[k];

					enumerator[6][k] = Alpha[k][ 6] + Beta[k + 1][ 3] + m11[k];

					enumerator[7][k] = Alpha[k][ 7] + Beta[k + 1][ 7] + m11[k];



					denominator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 0] - m11[k];

					denominator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 4] - m11[k];

					denominator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 1] - m10[k];

					denominator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 5] - m10[k];

					denominator[4][k] = Alpha[k][ 4] + Beta[k + 1][ 6] - m10[k];

					denominator[5][k] = Alpha[k][ 5] + Beta[k + 1][ 2] - m10[k];

					denominator[6][k] = Alpha[k][ 6] + Beta[k + 1][ 7] - m11[k];

					denominator[7][k] = Alpha[k][ 7] + Beta[k + 1][ 3] - m11[k];

				}





			}

			else

			{

                if (K == 5)

				{ 

					counter = 16;

                    // âû÷èñëåíèå ìåòðèê ñîñòî¤íèé Alpha

                 for (int k=0; k<L_total; k++)

				{

                        Alpha[k + 1][ 0] = maxi(Alpha[k][ 0] - m11[k], Alpha[k][ 1] + m10[k]);

                        Alpha[k + 1][ 1] = maxi(Alpha[k][ 3] - m10[k], Alpha[k][ 2] + m11[k]);

                        Alpha[k + 1][ 2] = maxi(Alpha[k][ 5] - m10[k], Alpha[k][ 4] + m11[k]);

                        Alpha[k + 1][ 3] = maxi(Alpha[k][ 6] - m11[k], Alpha[k][ 7] + m10[k]);

                        Alpha[k + 1][ 4] = maxi(Alpha[k][ 9] - m10[k], Alpha[k][ 8] + m11[k]);

                        Alpha[k + 1][ 5] = maxi(Alpha[k][ 10] - m11[k], Alpha[k][ 11] + m10[k]);

                        Alpha[k + 1][ 6] = maxi(Alpha[k][ 12] - m11[k], Alpha[k][ 13] + m10[k]);

                        Alpha[k + 1][ 7] = maxi(Alpha[k][ 15] - m10[k], Alpha[k][ 14] + m11[k]);

                        Alpha[k + 1][ 8] = maxi(Alpha[k][ 1] - m10[k], Alpha[k][ 0] + m11[k]);

                        Alpha[k + 1][ 9] = maxi(Alpha[k][ 2] - m11[k], Alpha[k][ 3] + m10[k]);

                        Alpha[k + 1][ 10] = maxi(Alpha[k][ 4] - m11[k], Alpha[k][ 5] + m10[k]);

                        Alpha[k + 1][ 11] = maxi(Alpha[k][ 7] - m10[k], Alpha[k][ 6] + m11[k]);

                        Alpha[k + 1][ 12] = maxi(Alpha[k][ 8] - m11[k], Alpha[k][ 9] + m10[k]);

                        Alpha[k + 1][ 13] = maxi(Alpha[k][ 11] - m10[k], Alpha[k][ 10] + m11[k]);

                        Alpha[k + 1][ 14] = maxi(Alpha[k][ 13] - m10[k], Alpha[k][ 12] + m11[k]);

                        Alpha[k + 1][ 15] = maxi(Alpha[k][ 14] - m11[k], Alpha[k][ 15] + m10[k]);

				 }

                    // âû÷èñëåíèå ìåòðèê ñîñòî¤íèé Beta

				for (int k=L_total-1; k > -1 ; k--)

				{

                        Beta[k][ 0] = maxi(Beta[k + 1][ 0] - m11[k], Beta[k + 1][ 8] + m11[k]);

                        Beta[k][ 1] = maxi(Beta[k + 1][ 8] - m10[k], Beta[k + 1][ 0] + m10[k]);

                        Beta[k][ 2] = maxi(Beta[k + 1][ 9] - m11[k], Beta[k + 1][ 1] + m11[k]);

                        Beta[k][ 3] = maxi(Beta[k + 1][ 1] - m10[k], Beta[k + 1][ 9] + m10[k]);

                        Beta[k][ 4] = maxi(Beta[k + 1][ 10] - m11[k], Beta[k + 1][ 2] + m11[k]);

                        Beta[k][ 5] = maxi(Beta[k + 1][ 2] - m10[k], Beta[k + 1][ 10] + m10[k]);

                        Beta[k][ 6] = maxi(Beta[k + 1][ 3] - m11[k], Beta[k + 1][ 11] + m11[k]);

                        Beta[k][ 7] = maxi(Beta[k + 1][ 11] - m10[k], Beta[k + 1][ 3] + m10[k]);

                        Beta[k][ 8] = maxi(Beta[k + 1][ 12] - m11[k], Beta[k + 1][ 4] + m11[k]);

                        Beta[k][ 9] = maxi(Beta[k + 1][ 4] - m10[k], Beta[k + 1][ 12] + m10[k]);

                        Beta[k][ 10] = maxi(Beta[k + 1][ 5] - m11[k], Beta[k + 1][ 13] + m11[k]);

                        Beta[k][ 11] = maxi(Beta[k + 1][ 13] - m10[k], Beta[k + 1][ 5] + m10[k]);

                        Beta[k][ 12] = maxi(Beta[k + 1][ 6] - m11[k], Beta[k + 1][ 14] + m11[k]);

                        Beta[k][ 13] = maxi(Beta[k + 1][ 14] - m10[k], Beta[k + 1][ 6] + m10[k]);

                        Beta[k][ 14] = maxi(Beta[k + 1][ 15] - m11[k], Beta[k + 1][ 7] + m11[k]);

                        Beta[k][ 15] = maxi(Beta[k + 1][ 7] - m10[k], Beta[k + 1][ 15] + m10[k]);

				}









				for (int k=0; k<L_total;k++)

				{

                    enumerator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 8] + m11[k];

                    enumerator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 0] + m10[k];

                    enumerator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 1] + m11[k];

                    enumerator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 9] + m10[k];

                    enumerator[4][k] = Alpha[k][ 4] + Beta[k + 1][ 2] + m11[k];

                    enumerator[5][k] = Alpha[k][ 5] + Beta[k + 1][ 10] + m10[k];

                    enumerator[6][k] = Alpha[k][ 6] + Beta[k + 1][ 11] + m11[k];

                    enumerator[7][k] = Alpha[k][ 7] + Beta[k + 1][ 3] + m10[k];

                    enumerator[8][k] = Alpha[k][ 8] + Beta[k + 1][ 4] + m11[k];

                    enumerator[9][k] = Alpha[k][ 9] + Beta[k + 1][ 12] + m10[k];

                    enumerator[10][k] = Alpha[k][ 10] + Beta[k + 1][ 13] + m11[k];

                    enumerator[11][k] = Alpha[k][ 11] + Beta[k + 1][ 5] + m10[k];

                    enumerator[12][k] = Alpha[k][ 12] + Beta[k + 1][ 14] + m11[k];

                    enumerator[13][k] = Alpha[k][ 13] + Beta[k + 1][ 6] + m10[k];

                    enumerator[14][k] = Alpha[k][ 14] + Beta[k + 1][ 7] + m11[k];

                    enumerator[15][k] = Alpha[k][ 15] + Beta[k + 1][ 15] + m10[k];



                    denominator[0][k] = Alpha[k][ 0] + Beta[k + 1][ 0] - m11[k];

                    denominator[1][k] = Alpha[k][ 1] + Beta[k + 1][ 8] - m10[k];

                    denominator[2][k] = Alpha[k][ 2] + Beta[k + 1][ 9] - m11[k];

                    denominator[3][k] = Alpha[k][ 3] + Beta[k + 1][ 1] - m10[k];

                    denominator[4][k] = Alpha[k][ 4] + Beta[k + 1][ 10] - m11[k];

                    denominator[5][k] = Alpha[k][ 5] + Beta[k + 1][ 2] - m10[k];

                    denominator[6][k] = Alpha[k][ 6] + Beta[k + 1][ 3] - m11[k];

                    denominator[7][k] = Alpha[k][ 7] + Beta[k + 1][ 11] - m10[k];

                    denominator[8][k] = Alpha[k][ 8] + Beta[k + 1][ 12] - m11[k];

                    denominator[9][k] = Alpha[k][ 9] + Beta[k + 1][ 4] - m10[k];

                    denominator[10][k] = Alpha[k][ 10] + Beta[k + 1][ 5] - m11[k];

                    denominator[11][k] = Alpha[k][ 11] + Beta[k + 1][ 13] - m10[k];

                    denominator[12][k] = Alpha[k][ 12] + Beta[k + 1][ 6] - m11[k];

                    denominator[13][k] = Alpha[k][ 13] + Beta[k + 1][ 14] - m10[k];

                    denominator[14][k] = Alpha[k][ 14] + Beta[k + 1][ 15] - m11[k];

                    denominator[15][k] = Alpha[k][ 15] + Beta[k + 1][ 7] - m10[k];

				}



				}

				}

				}

				

		double* ar_enum	 = new double[counter];

		double* ar_denom = new double[counter];

		//std::cout<< counter<<std::endl;

		for (int i = 0; i< L_total; i++)

		{

			for (int j = 0; j< counter; j++)

			{

				*(ar_enum+ j) = enumerator[j][i];

				*(ar_denom + j) = denominator[j][i];



				//std::cout<< enumerator[j][i]<<"    "<< denominator[j][i]<<std::endl;



			}





			//std::cout<< maxim(ar_enum, counter)<<"   max   "<< maxim(ar_denom, counter)<<std::endl;



			*(L_all+i) = maxim(ar_enum, counter) - maxim(ar_denom, counter);

			



		}





		for (int i=0; i<L_total+1;i++)

		{

			delete [] Alpha[i];

		}

		for (int i=0; i<L_total+1;i++)

		{

			delete [] Beta[i];

		}	





		for (int i=0; i<cntr;i++)

		{

			delete [] enumerator[i];

			delete [] denominator[i];



		}

		delete [] Alpha;

		delete [] Beta;

		delete [] enumerator;

		delete [] denominator;

		delete [] ar_enum;

		delete [] ar_denom;

		delete [] m11;

		delete [] m10;

}



/* */void intrlv(double* x,int* alpha, int len, double* y)

{		for (int j=0; j<len;j++)

		{

            y[alpha[j] - 1] = x[j];

		}

}



/* */void intrlv2(double* x,int* alpha, int len, double* y)

{

	for (int j=0; j<len;j++)

	{

            y[j] = x[alpha[j] - 1];

	}

}



/* */void intrlv2(int* x,int* alpha, int len, int* y)

{

	for (int j=0; j<len;j++)

	{

            y[j] = x[alpha[j] - 1];

	}

}





/* */void demultiplex_new(double *r, int frame, double puncture, int code, int* puncture_pattern, double *sys_, double *parity1, double *parity2)

{



        //double *parity1 = new double[frame];

        //double *parity2 = new double[frame];

        //double *sys_ = new double[frame];

		int ptsize, pt1, pt2;

		int k=0;

        if (puncture == 1)

		{

			for  (int i=0; i<code;i++)

			{

				

				if ((i) % 3 == 0)

				{

				sys_[k]    = r[i];

				parity1[k] = r[i+1];

				parity2[k] = r[i+2];

				k +=1;

				}

			}

		}

        else

		{

            ptsize = puncture_pattern[0];

            pt1 = puncture_pattern[1];

            pt2 = puncture_pattern[2];



			for (int i=0; i<int(frame/ptsize);i++)

			{

				for (int k=0; k<ptsize; k++)

				{

					sys_[i * ptsize + k] = r[i * (ptsize + 2)+k];

				}

					parity1[i * (ptsize) + pt1] = r[i * (ptsize + 2) + ptsize + 0];

					parity2[i * (ptsize) + pt2] = r[i * (ptsize + 2) + ptsize + 1];

				

			}



		}

		//std::cout<<puncture_pattern[0]<<"  "<<puncture_pattern[1]<<"  "<<puncture_pattern[2]<<"  "<<std::endl;

}







 void decod(double* metrics_, int cod_length,int total_length,double puncture, int* puncture_patter, int* alpha, int* alpha_0, int niter, int terminate,int K,double* L_ext)

{

double* metrics = new double[cod_length];

int k=0;

for (int j=0; j<cod_length;j++)

{

	metrics[j] = metrics_[j];

}



double* r = new double[cod_length];

intrlv(metrics, alpha_0, cod_length,r);



double* sys1 = new double[total_length];

double* parity1 = new double[total_length];

double* parity2 = new double[total_length];



for (int i=0; i<total_length; i++)

{

	sys1[i] = 0;

	parity1[i] = 0;

	parity2[i] = 0;

}





demultiplex_new(r, total_length, puncture, cod_length, puncture_patter,sys1, parity1, parity2);



double* sys2 = new double[total_length];

intrlv2(sys1, alpha, total_length,sys2);



double* L_e = new double[total_length];

double* L_a = new double[total_length];

double* L_all = new double[total_length];

double* L_new = new double[total_length];

double* L_ext_ = new double[total_length];



for (int k=0;k<total_length;k++)

{

	L_e[k]=0;

	L_all[k]=0;

}



for (int iter=0;iter<niter; iter++)

{

     intrlv(L_e, alpha,total_length,L_a);

    logmapo(sys1, parity1,L_a, total_length, terminate, K, L_all);

	for (int k=0; k<total_length; k++)

	{

		L_e[k] = L_all[k] - L_a[k] - sys1[k];

	}

    intrlv2(L_e, alpha,total_length,L_a);

	logmapo(sys2, parity2,L_a, total_length, 0, K, L_all);

	for (int k=0; k<total_length; k++)

	{



		L_e[k] = L_all[k] - L_a[k] - sys2[k];

		L_new[k] = L_e[k] + L_a[k] + sys2[k];

	}

	intrlv(L_new, alpha, total_length, L_ext_);

}

for (int k=0;k<total_length;k++)

{

	*(L_ext+k)=L_ext_[k];

	

}

delete [] metrics;

delete []	r;

delete []	sys1;

delete []	sys2;

delete []	parity1;

delete []	parity2;

delete []	L_e;

delete []	L_a;

delete []	L_all;

delete []	L_new;

delete []	L_ext_;

}





void cicl(complex * arr,int len)

{

	complex tmp = arr[--len];

	while (len>0) arr[len--] = arr[len-1];

	arr[0]=tmp;

}





 void treatment(double* Snoisy_i, double* Snoisy_q, int lenSnoisy,int M, int Msc, int NFFT,

												 int Tpr, int Tpf, int Nsymb, int Nsymb_aft_pilot,

												 double* PilotRefi, double* PilotRefq, int flag_freq_div, 

												 int flag_correct,	 double *phi_scrmbl, int Num_otbr, int flag_angle_corr,

												 int flag_angle_mem, int flag_clip_afhi, int flag_noise_red, 

												 int n_wavelet, int clip_th_dB, int flag_fipspec,

												 int M_start, int flag_mean_afh,

												 double* ConDatai, double* ConDataq,

												 double* AFH_i, double* AFH_q, 

												 double* AFH_i_nr, double* AFH_q_nr)

{



	double pi = 3.141592653589793;

	complex *Snoisy = new complex[lenSnoisy];

	for (int i=0;i<lenSnoisy;i++)

	{

		Snoisy[i] =  Snoisy_i[i] + complex::j * Snoisy_q[i];

	}

	int T_tpr = NFFT+Tpr+Tpf;

	complex *ConData = new complex[int((Msc)*Nsymb/(flag_freq_div))];

	complex *SpilotIn = new complex[NFFT];

	complex *PilotIn = new complex[NFFT];

	complex *PilotInMsc = new complex[Msc];

	complex *PilotRef = new complex[NFFT];

	complex *AFHI = new complex[Msc];

	complex* tmp_Psb = new complex[int((Msc)/(flag_freq_div))];

	complex *Stmp1 = new complex[NFFT];

	complex *Sfft = new complex[NFFT];

	complex *Scon = new complex[Msc];

	complex *Scon_kostyl = new complex[Msc];

	complex *Scon1 = new complex[Msc];

	complex *Scon2 = new complex[int((Msc)/(flag_freq_div))];

	complex *Scon21 = new complex[int((Msc)/(flag_freq_div))];

	complex *SC = new complex[Msc];

	complex * AFH = new complex[Msc];

	complex * AFH_tmp = new complex[Msc];

	complex * AFH_0 = new complex[Msc+2];

	double *clip_arr = new double[NFFT];

	double *P_sb = new double[int((flag_freq_div))];

	double * AFH_zero_arr = new double[Msc];

	double * AFH_phi_arr = new double[Msc];

	double *AFH_0_re = new double[Msc+2];

	double *AFH_0_im = new double[Msc+2];

	double *AFH_0_re_w_tmp = new double[(Msc+2)];

	double *AFH_0_im_w_tmp = new double[(Msc+2)];

	complex meanAFH;

	double angle_afhi=0;

	double rmsAFH, clip_th;

	double sum_real=0;

	double sum_imag=0;

	double phi;

	double angle_mem=0;

	double alpha_mean=0;

	int * len_x_w = new int[];

	int cntr = 0;

	for (int k=0;k<Msc;k++)

	{

		PilotRef[k] = PilotRefi[k] + complex::j * PilotRefq[k];

	}	

	for (int ipilot = 0; ipilot<int(Nsymb/Nsymb_aft_pilot); ipilot++)

	{

			for (int k=0;k<NFFT;k++)

			{

				SpilotIn[k] = Snoisy[ipilot*T_tpr*(Nsymb_aft_pilot+1)+Tpr+k];

			}

			rfft.Forward(SpilotIn, PilotIn, NFFT);

			if (flag_fipspec == 1)

				{

					fliplr(PilotIn,NFFT);

					cicl(PilotIn,NFFT);

				}

			cntr=0;

			rmsAFH = 0;

			meanAFH = 0;



			int counterj=0;

			for (int k=0;k<Msc;k++)

			{

				PilotInMsc[k] = PilotIn[k+M_start+1] ;

			}



			for (int k=0;k<Msc;k++)

			{

				AFH_zero_arr[k] = 1;

				AFH_phi_arr[k] = 0;

				AFH[k] = PilotInMsc[k] / PilotRef[k] ;

				meanAFH += AFH[k];



				*(AFH_i + k) = AFH[k].re(); 

				*(AFH_q + k) = AFH[k].im();

			}

			meanAFH /= Msc;

			  //////////////////////////////////////////////////

			 ///////////////// velvet /////////////////////////

			//////////////////////////////////////////////////

			if (flag_noise_red*(1-flag_mean_afh) == 1)

			{

				for (int k=0;k<Msc;k++)

				{

					AFH_tmp[k] = AFH[k] ;

				}

				for (int k=0;k<Msc;k++)

				{

					//AFH[k] -= meanAFH;

					AFH_0[k+1] =AFH_tmp[k] /*-  meanAFH*/ ;

				}



				AFH_0[0] = AFH_tmp[0];





				AFH_0[Msc+1] = AFH_tmp[Msc-1];







				for (int k=0;k<Msc+2;k++)

				{

					AFH_0_re[k] =AFH_0[k].re();

					AFH_0_im[k] =AFH_0[k].im();

				}



				kostyleviywden(AFH_0_re,Msc+2,n_wavelet,AFH_0_re_w_tmp);



				kostyleviywden(AFH_0_im,Msc+2,n_wavelet,AFH_0_im_w_tmp);



				for (int k=0;k<Msc;k++)

				{

					AFH_tmp[k] =AFH_0_re_w_tmp[k+1] + complex::j * AFH_0_im_w_tmp[k+1];

				}



				for (int k=0;k<Msc;k++)

				{

					AFH[k] = AFH_tmp[k] ;

				}



			}



			for (int k=0;k<Msc;k++)

			{

				rmsAFH += sqrt(AFH[k].norm());

			}

			rmsAFH = rmsAFH/Msc;

			clip_th = pow(10,((20*log10(rmsAFH)+clip_th_dB)/20));

			int kkkkk =0;

			if (flag_clip_afhi == 1)

			{

				for (int k=0; k<Msc;k++)

				{

						if (sqrt(AFH[k].norm()) < clip_th)

						{

							//phi_afhi = AFHI[k].re()/AFHI[k].im();

							//angle_afhi = atan(phi_afhi);

							//AFHI[k] = clip_th * (cos(angle_afhi) - complex::j * sin(angle_afhi));

						AFH_zero_arr[k] = 0;	//AFH[k] = 0;

						AFH_phi_arr[k] = atan(AFH[k].re()/AFH[k].im());

						kkkkk += 1;

						}

				}

			}

			//printf("AFH %d \n",ipilot);



			for (int k=0; k<Msc;k++)

			{

				if (flag_mean_afh*flag_noise_red == 1)

				{

					AFH[k] = meanAFH; //rmsAFH*(cos(angle_afhi) + complex::j * sin(angle_afhi));

					AFHI[k] = 1/(AFH[k]); // * AFH_zero_arr[k]+ (1-AFH_zero_arr[k])*clip_th*(cos(AFH_phi_arr[k]) + complex::j * sin(AFH_phi_arr[k])));



				}

				else

				{

				AFHI[k] = 1/(AFH[k] * AFH_zero_arr[k]+ (1-AFH_zero_arr[k])*clip_th*(cos(AFH_phi_arr[k]) + complex::j * sin(AFH_phi_arr[k])));

				}

				*(AFH_i_nr + k) = AFH[k].re() * AFH_zero_arr[k] + (1-AFH_zero_arr[k])*clip_th;

				*(AFH_q_nr + k) = AFH[k].im() * AFH_zero_arr[k];

			}

			//printf("AFH done %d \n",ipilot);



			angle_mem=0;

			for (int i=0; i<int(Nsymb_aft_pilot);i++)

			{

				for (int k=0;k<NFFT;k++)

				{

					Stmp1[k] = Snoisy[ipilot*T_tpr*(Nsymb_aft_pilot+1)+ T_tpr+ i * T_tpr + Tpr + k];

			

				}



				rfft.Forward(Stmp1, Sfft, NFFT);

				if (flag_fipspec == 1)

				{

					fliplr(Sfft,NFFT);

					cicl(Sfft,NFFT);

				}

				//printf("fftshechka %d %d\n",i,ipilot);



				for (int k=0;k<Msc;k++)

				{

					Scon[k] = Sfft[k+M_start+1] ;

				}



				if (flag_correct == 1)

				{

					for (int k=0; k<Msc; k++)

					{

					Scon[k] *= AFHI[k];



					}

				}



				for (int k=0; k<Msc; k++)

				{

					Scon_kostyl[k] = Scon[k];

				}

		

				double tmp_abs_Psb = 0;



				for (int k=0; k<int((flag_freq_div)); k++)

				{

					slice(Scon_kostyl,Msc, int(k*Msc/(flag_freq_div)), int((k+1)*Msc/(flag_freq_div)), tmp_Psb) ;

					P_sb[k] = 0;

					for (int jkj = 0;jkj <(Msc)/(flag_freq_div);jkj++)

					{

						P_sb[k] +=  sqrt(tmp_Psb[jkj].norm());

						Scon2[jkj] = 0;

					}

				}

				//printf("gogtovimsja k razneseniju %d %d\n",i,ipilot);



				if (flag_freq_div > 1)

				{

					for (int k=0; k<Msc; k++)

					{

						Scon[k] *= (cos(phi_scrmbl[k]) - complex::j * sin(phi_scrmbl[k]));

					}

				}



				for (int k=0; k<Msc; k++)

				{

					Scon_kostyl[k] = Scon[+k];

				}



				int* W = new int [int((flag_freq_div))];

				int imin = 0;



					for (int k = 0; k<int((flag_freq_div));k++)

					{

						W[k] = 1;

					}

					for (int k=0;k<Num_otbr;k++)

					{

						imin = argminim(P_sb,int((flag_freq_div)));	

						W[imin] = 0;

						P_sb[imin] = 10000;



					}



					for (int k = 0; k<int((flag_freq_div));k++)

					{

						slice(Scon_kostyl,Msc, int(k*Msc/(flag_freq_div)), int((k+1)*Msc/(flag_freq_div)), tmp_Psb);

						for (int ijk = 0;ijk<int((Msc)/(flag_freq_div)); ijk++)

						{

							Scon2[ijk] += W[k] * tmp_Psb[ijk] / ((flag_freq_div)-Num_otbr);

						}



					}



					//printf("raznesli %d %d\n",i,ipilot);



					if (flag_angle_mem == 1)

				{

					for (int k=0; k<int((Msc)/(flag_freq_div)); k++)

					{

							Scon2[k] = Scon2[k] * (cos(angle_mem) - complex::j * sin(angle_mem));

					}  

				}



				if (flag_angle_corr == 1)

				{



					if (M == 2)

					{

						for (int ji = 0; ji<(Msc)/(flag_freq_div); ji++)

						{

							if ( Scon2[ji].im() < 0)

							{

								Scon21[ji] = Scon2[ji] * (-1 );

							}

							else

							{

								Scon21[ji] = Scon2[ji];

							}

							sum_real += Scon21[ji].re();

							sum_imag += Scon21[ji].im();

						}



						phi = sum_real/sum_imag;

						alpha_mean = - atan(phi);

						angle_mem += alpha_mean;

					}

					else

					{

						for (int ji = 0; ji<(Msc)/(flag_freq_div); ji++)

						{

							if ((Scon2[ji].re() > 0) && (Scon2[ji].im()< 0))

							{

								Scon21[ji] = Scon2[ji] * complex::j;

							}

							else

							{

								if ((Scon2[ji].re() < 0) && (Scon2[ji].im()< 0))

								{

									Scon21[ji] = Scon2[ji] * (-1);

								}

								else

								{

									if ((Scon2[ji].re() < 0) && (Scon2[ji].im()> 0))

										{

										Scon21[ji] = Scon2[ji] * (0 - complex::j);

										} 

									else

										{

											Scon21[ji] = Scon2[ji];

										}

								}

							}

							sum_real += Scon21[ji].re();

							sum_imag += Scon21[ji].im();

						}

						phi = sum_real/sum_imag;

						alpha_mean = pi / 4 - atan(phi);

						angle_mem += alpha_mean;

					}

					for (int k=0; k<(Msc)/(flag_freq_div);k++)

					{

					Scon2[k] = Scon2[k] * (cos(alpha_mean) - complex::j * sin(alpha_mean));

					}

				}

				//printf("vse corecii %d %d\n",i,ipilot);



				for (int k=0; k<int((Msc)/(flag_freq_div)); k++)

				{

					ConData[(ipilot*int(Nsymb_aft_pilot))*int((Msc)/(flag_freq_div))+int((Msc)/(flag_freq_div))*i  + k] = Scon2[k];

				}

				//printf("sunul v condatu %d %d\n",i,ipilot);



			}

	}



	for (int k=0; k<int((Msc)*Nsymb/(flag_freq_div)); k++)

	{

		*(ConDatai+k) = ConData[k].re();

		*(ConDataq+k) = ConData[k].im();

	}

	//printf("poluchi vse vrode %d %d\n",0);



	delete [] Snoisy;

	delete []  ConData;

	delete []  SpilotIn;

	delete []  PilotIn ;

	delete []  PilotInMsc ;

	delete []  PilotRef;

	delete []  AFHI;

	delete []  Stmp1;

	delete []  Sfft;

	delete []  Scon ;

	delete []  Scon_kostyl;

	delete []  Scon1 ;

	delete []  SC ;

	delete [] P_sb;

	delete [] tmp_Psb;

	delete [] Scon2;

	delete [] Scon21;

	delete []  AFH_tmp ;

	delete []  AFH_zero_arr ;

	delete [] AFH_phi_arr;

	delete []  AFH_0 ;

	delete [] AFH_0_re;

	delete [] AFH_0_im;

	delete [] AFH_0_re_w_tmp;

	delete [] AFH_0_im_w_tmp;

		//printf("udalili vse vrode %d %d\n",0);



}









std::vector<complex> modQAM(std::vector<int> in, int k){



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

			mod[i] = complex::complex(I_2[in_dec[i]], Q_2[in_dec[i]]);

		}

		break;

	case 2:

		for (unsigned int i = 0; i < in_dec.size(); i++) {

			mod[i] = complex::complex(I_4[in_dec[i]], Q_4[in_dec[i]]);

		}

		break;

	case 4:

		for (unsigned int i = 0; i < in_dec.size(); i++) {

			mod[i] = complex::complex(I_16[in_dec[i]], Q_16[in_dec[i]]);

		}

		break;

	case 6:

		for (unsigned int i = 0; i < in_dec.size(); i++) {

			mod[i] = complex::complex(I_64[in_dec[i]], Q_64[in_dec[i]]);

		}

		break;

	default:

		break;

	}



	return mod;

}



 void SNR_est(double * in, int lenin, int wordsize, double * ConData_i,double * ConData_q, double * SNR)

{

	std::vector<int> data(lenin);

	for (int i=0;i<lenin; i++)

		data[i] = int(in[i]);

	std::vector<complex> ConData(lenin/wordsize);

	for (int i=0;i<lenin/wordsize; i++)

		ConData[i] = ConData_i[i] + complex::j * ConData_q[i];





	std::vector<complex> data_mod(lenin/wordsize);

	data_mod = modQAM(data, wordsize);

	// time 2 estimate



  //      SNR = 10 * np.log10((sumIr + sumQr) / (sumI + sumQ))



	double sumI = 0.0;

	double sumQ = 0.0;

	double sumIr = 0.0;

	double sumQr = 0.0;



	for (int k = 0;k <lenin/wordsize;k++)

	{

		sumI += abs(data_mod[k].re() - ConData[k].re()) * abs(data_mod[k].re() - ConData[k].re());

		sumQ += abs(data_mod[k].im() - ConData[k].im()) * abs(data_mod[k].im() - ConData[k].im());

		sumIr += abs(ConData[k].re())*abs(ConData[k].re());

		sumQr += abs(ConData[k].im())*abs(ConData[k].im());

	}

	SNR[0] =  10 * log10((sumIr + sumQr) / (sumI + sumQ));



}



void rsc_encode_new(int* y, int* x_info, int* x, int niter, int info_length, int total_length, int* g, 

					double puncture, int cod_length, int K, int m, int* puncture_pattern,

					int* alpha, int* alpha_0, int terminate){



	

	int length_;

	if (terminate == 1)

	{

		length_ = info_length;

	} 

	else

	{

		length_ = total_length;

	}



	int state = 0; // initial coder state (1 in MATLAB)



	switch (K)

	{

	case 3:

		{

		const int state_array[2][4]		=	{{0, 2, 3, 1},

											{2, 0, 1, 3}};

		const int out_array[2][4]		=	{{0, 0, 1, 1},

											{1, 1, 0, 0}};

		const int terminate_array[4][2]	=	{{0, 0}, {1, 0}, 

											{0, 1}, {1, 1}};

        const int info_array[4][2]		=	{{0, 0}, {1, 0}, 

											{1, 1}, {0, 1}};

		for (int i = 0; i<length_; i++){

			y[i] = out_array[int(x[i])][state];

			state = state_array[int(x[i])][state];

		}

		

		if (terminate == 1)

		{

			for (int i = 0; i < m; i++)

			{

				y[i + length_]		= terminate_array[state][i];

				x_info[i + length_] = info_array[state][i];

			}

		}



		}

		break;

		

	case 4: 

		{

		const int state_array[2][8]		=	{{0, 4, 1, 5, 6, 2, 7, 3},

											{4, 0, 5, 1, 2, 6, 3, 7}};

		const int out_array[2][8]		=	{{0, 0, 1, 1, 1, 1, 0, 0},

											{1, 1, 0, 0, 0, 0, 1, 1}};

		const int terminate_array[8][3]	=	{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, 

											{0, 1, 1}, {1, 1, 1}, {1, 0, 1}, {0, 0, 1}};

        const int info_array[8][3]		=	{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, 

											{1, 0, 1}, {0, 0, 1}, {1, 1, 1}, {0, 1, 1}};





		for (int i = 0; i<length_; i++){

			y[i] = out_array[int(x[i])][state];

			state = state_array[int(x[i])][state];

		}



		if (terminate == 1)

		{

			for (int i = 0; i < m; i++)

			{

				y[i + length_]		= terminate_array[state][i];

				x_info[i + length_] = info_array[state][i];

			}

		}		



		}

		break;

		

	case 5: 

		{

		const int state_array[2][16]		=	{{0, 8, 9, 1, 10, 2, 3, 11, 12, 4, 5, 13, 6, 14, 15, 7},

												{8, 0, 1, 9, 2, 10, 11, 3, 4, 12, 13, 5, 14, 6, 7, 15}};

		const int out_array[2][16]			=	{{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1},

												{1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0}};

		const int terminate_array[16][4]	= 	{{0, 0, 0, 0}, {0, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0}, {0, 1, 0, 0}, {1, 1, 1, 0}, {1, 1, 1, 0},

												{0, 1, 1, 0}, {0, 1, 1, 0}, {0, 0, 1, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 0, 1, 0}};

		const int info_array[16][4] 		= 	{{0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0}, {1, 1, 1, 0}, {0, 1, 1, 0}, {0, 0, 1, 0}, {1, 0, 1, 0}, {1, 1, 1, 1}, {0, 1, 1, 1},

												{0, 0, 1, 1}, {1, 0, 1, 1}, {0, 0, 0, 1}, {1, 0, 0, 1}, {1, 1, 0, 1}, {0, 1, 0, 1}};

		

		for (int i = 0; i<length_; i++){

			y[i] = out_array[int(x[i])][state];

			state = state_array[int(x[i])][state];

		}



		if (terminate == 1)

		{

			for (int i = 0; i < m; i++)

			{

				y[i + length_]		= terminate_array[state][i];

				x_info[i + length_] = info_array[state][i];

			}

		}



		}

		break;

		

	}



}







 void encoderm_new(int* en_output, int* x, int niter, int info_length, int total_length, int* g, double puncture, int cod_length, 

				  int K, int m, int* puncture_pattern, int* alpha, int* alpha_0, int terminate){



	int L_total = total_length;



	int* y0 = new int[L_total];

	int* y1 = new int[L_total];

	int* y2 = new int[L_total];



	int* x_info = new int[L_total];

	int* x_intrl = new int[L_total];



	memcpy(x_info, x, sizeof(int)*info_length);

	rsc_encode_new(y1, x_info, x, niter, info_length, total_length, g, puncture, cod_length, 

		K, m, puncture_pattern, alpha, alpha_0, terminate);



	memcpy(y0, x_info, sizeof(int)*L_total);



	intrlv2(x_info, alpha, L_total, x_intrl);



	rsc_encode_new(y2, x_info, x_intrl, niter, info_length, total_length, g, puncture, cod_length, 

		K, m, puncture_pattern, alpha, alpha_0, 0);



	if (puncture == 1)

	{

		int j = 0;

		for (int i = 0; i < cod_length; i+=3)

		{

			en_output[i]	= y0[j];

			en_output[i+1]	= y1[j];

			en_output[i+2]	= y2[j];

			j++;

		}

	}

	else

	{

		int ptsize = puncture_pattern[0];

		int pt1 = puncture_pattern[1];

		int pt2 = puncture_pattern[2];



		for (int i = 0; i < L_total/ptsize; i++)

		{

			for (int j = 0; j < ptsize; j++)

			{

				en_output[i * (ptsize + 2)+j]			= y0[i * ptsize +j];

			}

			en_output[i * (ptsize + 2) + ptsize]		= y1[i * ptsize + pt1];

			en_output[i * (ptsize + 2) + ptsize + 1]	= y2[i * ptsize + pt2];

		}

	}



	delete [] y0;

	delete [] y1;

	delete [] y2;



	delete [] x_info;

	delete [] x_intrl;



}



 void createPack(double* outPack_i, double* outPack_q,int Nsymb, double* dataIn_i, double* dataIn_q, double* pilot_i, double* pilot_q, double* window, int Msc, int NFFT,

												  int Tpr, int Tpf, double* spec_koef_i, double* spec_koef_q, int flag, double* phi, int Nsymb_aft_pilot){



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

