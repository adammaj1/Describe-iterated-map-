
/*

	c++ code by marcm200
	

	determine numerically the cycles of a polynomial Julia set
	https://fractalforums.org/image-threads/25/julia-and-mandelbrot-sets-w-or-wo-lyapunov-sequences/2696/new;topicseen#new
	
	
	
	What method/algorithm  you use to compute period ?

1. For the numerical periodicity detection I first find the critical points via classical Newton using a set of equidistant starting points at a square the size 3 times the Lagrangian estimate for the location of the roots of f'(z) (usually a total of 1024 points). 

2. Then I construct the orbit of those using double precision (or sometimes float128), but without rounding control, for usually 25 000 iterations (escaping is tested using an escape radius either manually set or computed via the Douady estimate). 


3. At the end of the orbit computation, I check backwards whether the last iterate has already been visited before (squared Euclidean distance of two points below 10^-15) to get the cycle length.
	
	The C++ code is in the zip-file below. Principal functions are findeNullstellen(...) which searches for roots with the described rectangular starting set, newton(..) which performs the actual Newton iterations from one starting complex number and AnalyseFunction(...) which analyzes the critical points. The code is based on two structs, Complex for complex number arithmetics and Polynom that stores a complex polynomial in coefficient form. The function itself is located in an external file named _infunction.txt.
	


One thing I forgot to mention. The code uses the format specifier %lld in printf and sscanf, that doesn't work on all compilers, some call it %LLd I guess or something similar. So you might need to change all those 64-bit references back to %i and int32_t (for the polynomials I usually compute, 64-bit is not necessary, but I started to use them some time ago). Or just use floating point numbers in the _infunction.txt file (so al lines starting with "r"). That should take care of the warnings and the incorrectly computed values in the output you posted.
=======================
As for the license: From my side, there isn't one. If you want to use the code, by all means do so, no reference needed as it is mostly basic stuff (the actual work was the (modified) universal set of Newton starting points, but that's due to Sutherland, so I suggest keeping that reference in the comments).
============================


https://stackoverflow.com/questions/13590735/printf-long-long-int-in-c-with-gcc

---------- translation ======================================
german	(abbrev)	english
-----------------------------
ausgabe		output
grade			degree
nenner			denominator
Ableitung		derivative
Anzahl	( anz)		number
ergebnisse (erg) 	result	
====================================================

g++ main.cpp -lm -Wall
./a.out > rabbit_out.txt

	
	
*/

#include "stdio.h"
#include "stdint.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include <inttypes.h> 


const int32_t MAXDEGREE=1024; // highest polynomial degree
const int32_t MAXZYKLUS=32; // max number of cycles found
const int32_t MAXZYKLUSLEN=512; // max length of a cycle
const int32_t MAXCP=1024; // max number of critical points


/* (There is no special reason to have two values in my code for orbit len, so you can combine them and maybe use maxit throughout. 
near parabolic ones need a longer iteration. I recommend you change the following two lines in the code:

int32_t maxit=25000;

to some other value - 50000 is working here.

If you want to use a higher maxit than 2^16 you also need to change the memory allocation constant:

const int32_t MAXFIXORBITLEN=(1 << 16); // max length of a total orbit to some higher value than maxit.


ike this from a (much) slower computer to here and then I just left it that way).

I normally do not go beyond 50000, because I usually also compute the
image - and more than 50000 per pixel takes too long for those.


*/
const int32_t MAXFIXORBITLEN=50000; // (1 << 16); // max length of a total orbit = 2^16
int32_t maxit=MAXFIXORBITLEN;





int64_t DENOMINATOR=(int64_t)1 << 25; // = 1*2^25 ( left shift)

typedef double NTYP; // re-type to hold a higher precision number  type ???

NTYP ESCAPEQ=1E+9; // squared escape radius for the critical orbit construction
const NTYP NULLEPSILON=1E-15; // arbitrary choice how close the (Squared) distance of two "equal" numbers must get



// ------------------------ structures -------------------------------------------------------

typedef struct Complex {
	NTYP re,im;

	Complex(const NTYP r_=0.0, const NTYP i_=0.0);

	char* str(char*);
	friend bool operator==(const Complex&,const Complex&);
	bool operator!=(const Complex);
	void output(FILE*);
	void output225(FILE*);

	friend Complex operator+(const Complex&,const Complex&);
	friend Complex operator*(const double&,const Complex&);
	Complex operator*(const Complex);
	Complex operator*(const double);
	Complex& operator=(const Complex&);
	friend Complex operator/(const Complex&,const Complex&);

	NTYP norm(void);
	NTYP normQ(void);
} ComplexT;

typedef struct Polynom {
	int32_t degree;
	Complex coeff[MAXDEGREE];
	char coeffnull[MAXDEGREE];

	Polynom();
	void setCoeff(const int32_t,const Complex);
	void setCoeff(const int32_t,const NTYP,const NTYP);
	void setCoeff(const int32_t,const NTYP);
	void eval_arg_f(const Complex,Complex&); // über Horner
	void clearCoeff(void);
	void output(FILE*,const char* comment);
} PolynomT;

typedef struct OrbitFix {
	int32_t anz;
	Complex* pts;

	OrbitFix();
	virtual ~OrbitFix();
	void setP(const int,const Complex);
} OrbitFixT;


// globals

ComplexT nullstellen[MAXDEGREE]; // array for critical points
int32_t nullstellenanz;
OrbitFixT *zyklus=NULL; 


// functions


/*
gives argument of complex number in turns 
https://en.wikipedia.org/wiki/Turn_(geometry)

*/

double GiveTurn( double x, double y){
	double t;

 	t =  atan2(y,x); //carg(Complex(x,y));
  	t /= 2*M_PI; // now in turns
  	if (t<0.0) t += 1.0; // map from (-1/2,1/2] to [0, 1) 
  	return (t);
}



char* chomp(char* s) {
	if (!s) return 0;
	for(int i=strlen(s);i>=0;i--) if (s[i]<32) s[i]=0; else break;
	return s;
}

// Complex

Complex operator/(const Complex& a,const Complex& b) {
    const NTYP n2=b.re*b.re+b.im*b.im;
    Complex ret(
		(a.re*b.re+a.im*b.im)/n2,
		(a.im*b.re-a.re*b.im) / n2
    );

    return ret;
}


void Complex::output(FILE* f) {
	fprintf(f,"%.20lg%+.20lgi",re,im);
}


// 
void Complex::output225(FILE* f) {
	//fprintf(f,"%ld,%ld",(int64_t)floor(DENOMINATOR*re),(int64_t)floor(DENOMINATOR*im)); // l // not works
	fprintf(f,"%.20f,%.20f",re,im); // l
}



// cabs(z) 
NTYP Complex::norm(void) {
	return sqrt(re*re+im*im);
}

NTYP Complex::normQ(void) {
	return (re*re+im*im);
}

char* Complex::str(char* function) {
	sprintf(function,"%lg%+lgi",re,im);
	return function;
}

Complex& Complex::operator=(const Complex& c) {
	if (this != &c) {
		re=c.re;
		im=c.im;
	}

	return *this;
}

bool operator!=(const Complex a, const Complex b) {
    return ( (a.re != b.re) || (a.im != b.im) );
}

Complex operator+(const Complex& a, const Complex& b) {
    Complex ret(a.re+b.re, a.im+b.im);
    return ret;
}

Complex operator*(const double& a, const Complex& b) {
    Complex ret=Complex(a*b.re,a*b.im);
    
    return ret;
}

Complex operator-(const Complex a, const Complex b) {
    Complex ret(a.re-b.re, a.im-b.im);
    return ret;
}

Complex Complex::operator*(const Complex a) {
    Complex ret(a.re*re-a.im*im, a.im*re+a.re*im);

    return ret;
}

Complex Complex::operator*(const double a) {
    Complex ret(a*re, a*im);

    return ret;
}

Complex::Complex(const NTYP r, const NTYP i) {
	re=r;
	im=i;
}

// Polynom

void Polynom::eval_arg_f(const Complex az,Complex& function) {
	function=coeff[degree];

	for(int32_t i=degree;i>0;i--) {
		function = function*az;
		function = function + coeff[i-1];
	}
}

Polynom::Polynom() {
	degree=0;
	clearCoeff();
}



// ------------------ setCoeff ----------------------------------------------
void Polynom::setCoeff(const int32_t aidx,const NTYP ar,const NTYP ai) {
	setCoeff(aidx,Complex(ar,ai));
}

void Polynom::setCoeff(const int32_t aidx,const NTYP ar) {
	setCoeff(aidx,Complex(ar,0.0));
}

void Polynom::setCoeff(const int32_t aidx,const Complex ac) {
	coeff[aidx]=ac;
	if (coeff[aidx].normQ() < 1E-20) {
		coeffnull[aidx]=1;
	} else {
		coeffnull[aidx]=0;
		if (aidx > degree) degree=aidx;
	}
}

// -----------------------------------------



void Polynom::clearCoeff(void) {
	for(int32_t i=0;i<MAXDEGREE;i++) {
		coeffnull[i]=1;
		coeff[i]=Complex(0.0,0.0);
	}
	degree=0;
}





//

void Polynom::output(FILE* f,const  char *comment) {
	
	if (degree>32) {
		printf("Program not works for degree %i > 32 \n",degree);
		return;
	}
	
	
	int32_t erster=1;
	
	fprintf(f,"%s", comment);
	for(int32_t i=degree;i>=0;i--) {
		if (coeffnull[i]==0) {
			if (erster) erster=0; else fprintf(f,"+");
			fprintf(f,"(");
			coeff[i].output(f);
			if (i != 0) fprintf(f,")*z^%i",i);
			else fprintf(f,")");
		}
	}
	fprintf(f,"\n\n");
	
	// output255
	
	//fprintf(f,"output255 : p(z)=");
	//if (degree>32) {
	//	printf("Program not works for degree %i > 32 \n",degree);
	//	return;
	//}
	/*
	erster=1;
	for(int32_t i=degree;i>=0;i--) {
		if (coeffnull[i]==0) {
			if (erster) erster=0; else fprintf(f,"+");
			fprintf(f,"(");
			coeff[i].output225(f);
			if (i != 0) fprintf(f,")*z^%i",i);
			else fprintf(f,")");
		}
	}
	fprintf(f,"\n");
	*/
	
}

/* compute derivative of a polynomial infunction
input inkft

*/
void ableitenFA(Polynom  & infunction,Polynom& functionabl) {
	functionabl.clearCoeff();

	for(int32_t i=1;i<=infunction.degree;i++) {
		if (infunction.coeffnull[i]==0) {
			functionabl.setCoeff(
				i-1,
				Complex(i,0.0)*infunction.coeff[i]
			);
		}
	}
}

// starting from a specific point: perform Newton iteration
int32_t newton(
	Polynom& afunction,Polynom& aabl,
	const Complex astart,Complex& function)
{
	Complex z=astart,f,fabl,zletzt,d;
	int32_t ret=0;

	for(int32_t i=1;i<maxit;i++) {
		zletzt=z;
		// z(n+1)=zn-f(zn) / f'(zn);
		afunction.eval_arg_f(z,f);
		aabl.eval_arg_f(z,fabl);
		z = z - f / fabl;
		d = z-zletzt;
		if (d.normQ() < NULLEPSILON) {
			function=z;
			ret=i;
			return i;
		}
	}

	return ret;
}

// search for a found root in the current list of
// already found roots => already known?
int32_t getNullstellenIdx(
	const Complex aw,
	const int32_t ait
) {
	for(int32_t i=0;i<nullstellenanz;i++) {
		Complex d=nullstellen[i] - aw;
		if (d.normQ() < NULLEPSILON) {
			return i;
		}
	}
	// new root found
	if (nullstellenanz > (MAXDEGREE-8)) {
		printf("Error. Too many roots\n");
		exit(99);
	}

	nullstellen[nullstellenanz]=aw;
	nullstellenanz++;
	return (nullstellenanz-1);
}



void OrbitFix::setP(const int idx,const Complex A) {
	pts[idx]=A;
}

OrbitFix::~OrbitFix() {
	delete[] pts;
}

OrbitFix::OrbitFix() {
	anz=0;
	pts=new Complex[MAXFIXORBITLEN];
}

void findeNullstellen(Polynom& function) {
	// region where all roots lie
	// Lagrange estimate
	// https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
	NTYP NN=function.coeff[function.degree].norm();
	NTYP ESCAPER=0.0;
	for(int32_t i=0;i<=function.degree;i++) {
		ESCAPER += (function.coeff[i].norm()/NN);
	}
	if (ESCAPER < 1.0) ESCAPER=1.0;
	ESCAPER=ceil(ESCAPER);

	Polynom derivative,abl2;
	ableitenFA(function,derivative);
	ableitenFA(derivative,abl2);

	// use a modified version of the universal set of
	// starting points for Newton iteration
	// by Schleicher and Sutherland
	// https://fractalforums.org/fractal-mathematics-and-new-theories/28/finding-all-roots-of-a-polynomial-with-mathematical-guarantee/2959
	// 3times the LAgrangian
	int32_t LEN=256; // spacing
	// the article provides a definite number of necessary
	// starting points on a circle
	NTYP SCRRE0=-3*ESCAPER;
	//NTYP SCRRE1= 3*ESCAPER;
	NTYP SCRIM0=-3*ESCAPER;
	NTYP SCRIM1= 3*ESCAPER;
	NTYP sk=SCRIM1-SCRIM0; sk /= LEN;

	Complex start;
	nullstellenanz=0;
	
	// check the border of a rectangle around the
	// Lagrange estimate
	#define SUCHE(X0,Y0,X1,Y1)\
	{\
		int32_t yd=1,xd=1;\
		if ((Y0) > (Y1)) yd=-1;\
		if ((X0) > (X1)) xd=-1;\
		for(int32_t y=(Y0);y<=(Y1);y+=yd) {\
			start.im=y*sk + SCRIM0;\
			for(int32_t x=(X0);x<=(X1);x+=xd) {\
				start.re=x*sk + SCRRE0; \
				\
				Complex nulls;\
				int32_t it=newton(function,derivative,start,nulls);\
				if (it > 0) {\
					getNullstellenIdx(nulls,it);\
					if (nullstellenanz >= function.degree) break;\
				} \
			} \
			if (nullstellenanz >= function.degree) break;\
		} \
	}

	SUCHE(0,0,0,LEN-1)
	if (nullstellenanz < function.degree) {
		SUCHE(0,LEN-1,LEN-1,LEN-1)
	}
	if (nullstellenanz < function.degree) {
		SUCHE(LEN-1,LEN-1,LEN-1,0)
	}
	if (nullstellenanz < function.degree) {
		SUCHE(LEN-1,0,0,0)
	}
}




// analyze critical points whether they are periodic
void AnalyseFunction(Polynom& function) {
	int32_t anzzyklus=0;
	zyklus= new OrbitFixT[MAXZYKLUS];

	ComplexT CP[MAXCP]; // array of critical points 

	PolynomT derivative;// compute polynomial derivative
	ableitenFA(function,derivative);
	derivative.output(stdout, "derivative dp/dz = "); // polynomial output
	

	// finding critical points = the points of an equation where the derivative is zero
	findeNullstellen(derivative);

	// copy zeros to CP-array
	int32_t cpanz=0; // number of critical points
	for(int32_t i=0;i<nullstellenanz;i++) {
		CP[cpanz]=nullstellen[i];
		cpanz++;
	}

	printf("%i critical points found\n",cpanz);

	OrbitFix orbit;

	// checking a cycle if already found
	#define ZYKLUSDA(PUNKT,SCHON)\
	{\
		SCHON=-1;\
		for(int32_t z=0;z<anzzyklus;z++) {\
			for(int32_t i=0;i<zyklus[z].anz;i++) {\
				Complex d=PUNKT - zyklus[z].pts[i];\
				if (d.normQ() < NULLEPSILON) {\
					SCHON=z;\
					break;\
				}\
			}\
			if (SCHON >= 0) break;\
		} \
	}

	// checking an orbit if periodic
	#define ZYKLISCH(ERG) \
	{\
		ERG=-1;\
		int32_t e=orbit.anz-MAXZYKLUSLEN;\
		if (e < 0) e=0;\
		for(int32_t i=(orbit.anz-2);i>=e;i--) {\
			Complex d=orbit.pts[orbit.anz-1] - orbit.pts[i];\
			if (d.normQ() < NULLEPSILON) {\
				ERG=i;\
				break;\
			}\
		} \
	}

	// add cycle if new
	#define ZYKLUSANBAUEN(PER0) \
	{\
		if (anzzyklus>=(MAXZYKLUS-2)) {\
			printf("Error. Too many cycles\n");\
			exit(99);\
		}\
		else {\
			zyklus[anzzyklus].anz=0;\
			int32_t a=0;\
			Complex der=Complex(1.0,0.0);\
			for(int32_t i=(PER0+1);i<orbit.anz;i++) {\
				zyklus[anzzyklus].setP(a,orbit.pts[i]);\
				Complex ab;\
				derivative.eval_arg_f(orbit.pts[i],ab);\
				der=der*ab;\
				a++;\
			}\
			/* multiplier l */\
			NTYP l=der.norm();\
			NTYP t = GiveTurn(der.re , der.im);\
			int32_t aus=0;\
			if (\
				(l <= 1.0) /* attracting and parabolic */\
			) { \
				aus=1; \
			}\
			if (aus>0) {\
				printf("and enters cycle #%i length=%i and it's stability = |multiplier|=%.5lg =",\
				anzzyklus,a,(double)l);\
				if (l<0.99) printf("attractive \n");\
				else {\
					if (l>1.0001) printf("repelling \n");\
						else printf("attracting or parabolic \n");\
					}\
				printf("\tinternal angle = %.20lg\n", t);\
				printf("cycle = {\n");\
				for(int32_t p=0;p<a;p++) {\
					printf("%.20lg,%.20lg ; ",\
					(double)zyklus[anzzyklus].pts[p].re,(double)zyklus[anzzyklus].pts[p].im);\
				}\
				printf("}\n");\
			}\
			if (aus>0) {\
				zyklus[anzzyklus].anz=a;\
				anzzyklus++;\
			}\
		}\
	}\
	
	
	
	

	for(int32_t cpi=0;cpi<cpanz;cpi++) {
		printf("\n\tcp#%i: %.20lg,%.20lg . It's critical orbit ",cpi,CP[cpi].re,CP[cpi].im);
		Complex z0=CP[cpi];

		int32_t esc=0;
		Complex zn=z0;
		for(int32_t i=0;i<maxit;i++) {
			orbit.setP(i,zn);
			orbit.anz=i+1;
			if (zn.normQ() > ESCAPEQ) {
				esc=1;
				break;
			}
			Complex znfunction;
			function.eval_arg_f(zn,znfunction);
			zn=znfunction;
		}

		if (esc>0) {
			printf("is escaping\n");
			continue;
		}
		printf("is bounded ");

		int32_t per0=-1;
		ZYKLISCH(per0)
		if (per0 < 0) {
			printf("is not (yet) periodic\n");
			continue;
		}

		int32_t schon;
		ZYKLUSDA(orbit.pts[orbit.anz-1],schon);

		if (schon >= 0) {
			printf(" and enters cycle #%i\n",schon);
			continue;
		}

		ZYKLUSANBAUEN(per0)

	}// cpi
	
	// to do 
	//printf("min distance between critical points = 
} // void AnalyseFunction












// read coefficients from text file
// problem : bad char in input string is not recognized
int32_t loadPolynom(Polynom  &function,  const char* FileName) {

	FILE *f=fopen(FileName,"rt");
	if (!f) {
		fprintf(stderr, "input file %s not found\n", FileName);
		return 0; 
		}
	
	
	char tmp[1024];
	function.clearCoeff();
	NTYP DN=1.0; 
	DN /= DENOMINATOR; //  = 1.0/2^25
	
	fprintf(stdout, "coefficients read from input file %s\n", FileName);
		

	while (!feof(f)) {
		int32_t idx; // degree of monomial
		int64_t re,im;
		fgets(tmp,1000,f);
		if (tmp[0]=='.') break;
		if (tmp[0]=='-') continue;

		chomp(tmp);
		
		if ((tmp[0]=='r')||(tmp[0]=='R')) // coefficient as 2 floating point numbers
		
		{
			double a,b;
			if (sscanf(&tmp[1],"%i,%lf,%lf",&idx,&a,&b) == 3) {
			function.setCoeff(
				idx,
				Complex(
					a,
					b
				)
						
			);
			fprintf(stdout, "\tdegree %i coefficient = ( %+.16lf %+.16lf*i) \n", idx, a,b);
			}
		// coefficients as a 2 rational numbers 	
		} else {
			if (sscanf(tmp,"%i,%ld,%ld",&idx,&re,&im) == 3) // l
			{
			function.setCoeff(
				idx,
				Complex(
					DN * re,  // re/(2^25)
					DN * im // im/(2^25)
				)
			);
			fprintf(stdout, "\tdegree %i coefficient = ( %+ld %+ld*i) / 2^25\n", idx, re, im);
			
			}
		}
	}
	fprintf(stdout, "\n");
	
	function.output(stdout, "Input polynomial p(z)="); // polynomial output
	fclose(f);

	return 1; 
}






int32_t main(int argc, char *argv[]) {


	

	char* FileName;
	PolynomT function;
	
	
	// read argument : input file 
	if( argc != 2 ) 
	{ 
		fprintf(stderr, "bad number of arguments ( should be one string = File Name)\n");
		fprintf(stderr, "to print the result to stdout use: \n");
		fprintf(stderr, "./a.out \"fileName\"  \n");
		fprintf(stderr, "or to print the result to the file use: \n");
		fprintf(stderr, "./a.out \"a.txt\" > a_out.txt  \n"); 
		return 1;
	} // check arguments
	
	FileName = argv[1]; // file name is an argument ot the program
	
	if (loadPolynom(function,FileName))
	//labirynth.txt"))// reads the text file for the polynomial to analyze
	{
		
		AnalyseFunction(function);
	}
	

	return 0;
}




