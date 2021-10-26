//   fft.h - declaration of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

//   Include complex numbers header
#include "complex.h"

//   FORWARD FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
bool Forward(const complex *const Input, complex *const Output, const unsigned int N);

//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
//     Data - both input data and output
//     N    - length of input data
bool Forward(complex *const Data, const unsigned int N);

//   INVERSE FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
//     Scale  - if to scale result
bool Inverse(const complex *const Input, complex *const Output, const unsigned int N, const bool Scale = true);

//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
//     Data  - both input data and output
//     N     - length of both input data and result
//     Scale - if to scale result
bool Inverse(complex *const Data, const unsigned int N, const bool Scale = true);
bool test(int N);

//   Rearrange function and its inplace version
static void Rearrange(const complex *const Input, complex *const Output, const unsigned int N);
static void Rearrange(complex *const Data, const unsigned int N);

//   FFT implementation
static void Perform(complex *const Data, const unsigned int N, const bool Inverse = false);

//   Scaling of inverse FFT result
//static void Scale(complex *const Data, const unsigned int N);



bool Forward(const complex *const Input, complex *const Output, const unsigned int N)

{

	//   Check input parameters

	if (!Input || !Output || N < 1 || N & (N - 1))

		return false;

	//   Initialize data

	Rearrange(Input, Output, N);

	//   Call FFT implementation

	Perform(Output, N);

	//   Succeeded

	return true;

}



//   FORWARD FOURIER TRANSFORM, INPLACE VERSION

//     Data - both input data and output

//     N    - length of input data

bool Forward(complex *const Data, const unsigned int N)

{

	//   Check input parameters

	if (!Data || N < 1 || N & (N - 1))

		return false;

	//   Rearrange

	Rearrange(Data, N);

	//   Call FFT implementation

	Perform(Data, N);

	//   Succeeded

	return true;

}



//   INVERSE FOURIER TRANSFORM

//     Input  - input data

//     Output - transform result

//     N      - length of both input data and result

//     Scale  - if to scale result

bool Inverse(const complex *const Input, complex *const Output, const unsigned int N, const bool Scale /* = true */)

{

	//   Check input parameters

	if (!Input || !Output || N < 1 || N & (N - 1))

		return false;

	//   Initialize data

	Rearrange(Input, Output, N);

	//   Call FFT implementation

	Perform(Output, N, true);

	//   Scale if necessary

	//if (Scale)

		//Scale(Output, N);

	//   Succeeded

	return true;

}



//   INVERSE FOURIER TRANSFORM, INPLACE VERSION

//     Data  - both input data and output

//     N     - length of both input data and result

//     Scale - if to scale result

bool Inverse(complex *const Data, const unsigned int N, const bool Scale /* = true */)

{

	//   Check input parameters

	if (!Data || N < 1 || N & (N - 1))

		return false;

	//   Rearrange

	Rearrange(Data, N);

	//   Call FFT implementation

	Perform(Data, N, true);

	//   Scale if necessary

    //Scale(Data, N);

	//   Succeeded

	return true;

}



//   Rearrange function

void Rearrange(const complex *const Input, complex *const Output, const unsigned int N)

{

	//   Data entry position

	unsigned int Target = 0;

	//   Process all positions of input signal

	for (unsigned int Position = 0; Position < N; ++Position)

	{

		//  Set data entry

		Output[Target] = Input[Position];

		//   Bit mask

		unsigned int Mask = N;

		//   While bit is set

		while (Target & (Mask >>= 1))

			//   Drop bit

			Target &= ~Mask;

		//   The current bit is 0 - set it

		Target |= Mask;

	}

}



//   Inplace version of rearrange function

void Rearrange(complex *const Data, const unsigned int N)

{

	//   Swap position

	unsigned int Target = 0;

	//   Process all positions of input signal

	for (unsigned int Position = 0; Position < N; ++Position)

	{

		//   Only for not yet swapped entries

		if (Target > Position)

		{

			//   Swap entries

			const complex Temp(Data[Target]);

			Data[Target] = Data[Position];

			Data[Position] = Temp;

		}

		//   Bit mask

		unsigned int Mask = N;

		//   While bit is set

		while (Target & (Mask >>= 1))

			//   Drop bit

			Target &= ~Mask;

		//   The current bit is 0 - set it

		Target |= Mask;

	}

}



//   FFT implementation

void Perform(complex *const Data, const unsigned int N, const bool Inverse /* = false */)

{

	const double pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;

	//   Iteration through dyads, quadruples, octads and so on...

	for (unsigned int Step = 1; Step < N; Step <<= 1)

	{

		//   Jump to the next entry of the same transform factor

		const unsigned int Jump = Step << 1;

		//   Angle increment

		const double delta = pi / double(Step);

		//   Auxiliary sin(delta / 2)

		const double Sine = sin(delta * .5);

		//   Multiplier for trigonometric recurrence

		const complex Multiplier(-2. * Sine * Sine, sin(delta));

		//   Start value for transform factor, fi = 0

		complex Factor(1.);

		//   Iteration through groups of different transform factor

		for (unsigned int Group = 0; Group < Step; ++Group)

		{

			//   Iteration within group 

			for (unsigned int Pair = Group; Pair < N; Pair += Jump)

			{

				//   Match position

				const unsigned int Match = Pair + Step;

				//   Second term of two-point transform

				const complex Product(Factor * Data[Match]);

				//   Transform for fi + pi

				Data[Match] = Data[Pair] - Product;

				//   Transform for fi

				Data[Pair] += Product;

			}

			//   Successive transform factor via trigonometric recurrence

			Factor = Multiplier * Factor + Factor;

		}

	}

}



//   Scaling of inverse FFT result
/*
void Scale(complex *const Data, const unsigned int N)

{

	const double Factor = 1. / double(N);

	//   Scale all data entries

	for (unsigned int Position = 0; Position < N; ++Position)

		Data[Position] *= Factor;

}
*/


bool test(int N)

{

            //printf("Test %d  samples\n", N);

            return true;

}
