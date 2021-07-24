<p align="center">
  <a href="https://github.com/adammaj1/Describe-iterated-map-"><img src="https://img.shields.io/tokei/lines/github/adammaj1/Describe-iterated-map-" alt="Repository Size"></a>
  <a href="https://github.com/adammaj1/Describe-iterated-map-/releases"><img src="https://img.shields.io/github/downloads/adammaj1/Describe-iterated-map-/total?style=flat" alt="Github Release"></a>
</p>

[![Build status](https://github.com/Anaconda-Platform/anaconda-project/workflows/Build%20and%20test/badge.svg)](https://github.com/Anaconda-Platform/anaconda-project/actions)


# numerical periodicity detection of a polynomial Julia set

determine numerically the cycles 

# algorithm

## find critical points
first find [the critical points](https://en.wikipedia.org/wiki/Critical_point_(mathematics)) via [classical Newton method](https://en.wikibooks.org/wiki/Fractals/Mathematics/Newton_method) using a set of [equidistant starting points](https://gitlab.com/adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method) (  at a square the size 3 times [the Lagrangian estimate for the location of the roots of](https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds) f'(z) (usually a total of 1024 points). 

## compute the critical orbits
Then I construct the orbit of those using double precision (or sometimes float128), but without rounding control, for usually 25 000 iterations (escaping is tested using an escape radius either manually set or computed via the Douady estimate). 


## find periods
At the end of the orbit computation, I check backwards whether the last iterate has already been visited before (squared Euclidean distance of two points below 10^-15) to get the cycle length.

# Rational maps

One can use the code for polynomial cycles to also detect attracting/parabolic ones numerically in the rational case by making some changes:
* Load 2 polynomials f and g to denote the rational f/g.(routien "loadPolynom")
* Find the zeros of g (poles) using the "findeNullstellen" .routine.
* Construct symbolically the derivative as two polynomials (f'*g-g'*f) and (g*g) using "ableitenFA" for the polynomial derivative (I have a polynomial multiplication and addition routine here in case you haven't implemented those yourself)
* Find the zeros of (f'*g-g'*f) again using "findeNullstellen" (critical points)
* Construct numerically the orbits of all those values (critical points and poles) to get (attracting) cycles - this needs some work as it's best to implement this anew (from the routien "analysiereJulia"), implementing a check for inf/NaN as a numerical value. If needed,
evaluating the symbolic derivative polynomials to get the multiplier (function Polynom::eval_arg_f) 
* Max iteration and epsilon might need to be adjusted as the complex division constitutes another layer of numerical instability.




# use 
Program reads parameters from the input txt file and prints the output data (to the stdout or to the output txt file)


```c
g++ main.cpp -lm -Wall
./a.out a.txt > a_out.txt
```

or if one wants use all input files in the working directory use simply [make](./src/Makefile) :  

```bash  
make
```   

Example result:



```bash
chmod +x m.sh
./m.sh
program compiled without errors  
run the program for input txt file:
	almost_failed_out.txt is output file = skipped 
	almost_failed.txt = input file 	 almost_failed_out.txt = output file 
	cubic_p2_Bagula_out.txt is output file = skipped 
	cubic_p2_Bagula.txt = input file 	 cubic_p2_Bagula_out.txt = output file 
	cubic_parab_out.txt is output file = skipped 
	cubic_parab.txt = input file 	 cubic_parab_out.txt = output file 
	DancingRabbit_out.txt is output file = skipped 
	DancingRabbit.txt = input file 	 DancingRabbit_out.txt = output file 
	dumbbell_out.txt is output file = skipped 
	dumbbell.txt = input file 	 dumbbell_out.txt = output file 
	fractile3_out.txt is output file = skipped 
	fractile3.txt = input file 	 fractile3_out.txt = output file 
	ijon_b012_out.txt is output file = skipped 
	ijon_b012.txt = input file 	 ijon_b012_out.txt = output file 
	kawahira_sc_c3_out.txt is output file = skipped 
	kawahira_sc_c3.txt = input file 	 kawahira_sc_c3_out.txt = output file 
	labirynth_out.txt is output file = skipped 
	labirynth.txt = input file 	 labirynth_out.txt = output file 
	milnor_fig11_out.txt is output file = skipped 
	milnor_fig11.txt = input file 	 milnor_fig11_out.txt = output file 
	notanewrecord_out.txt is output file = skipped 
	notanewrecord.txt = input file 	 notanewrecord_out.txt = output file 
	p3near_out.txt is output file = skipped 
	p3near.txt = input file 	 p3near_out.txt = output file 
	q_p2_1_out.txt is output file = skipped 
	q_p2_1.txt = input file 	 q_p2_1_out.txt = output file 
	q_p2_2_out.txt is output file = skipped 
	q_p2_2.txt = input file 	 q_p2_2_out.txt = output file 
	q_primitive3o7_out.txt is output file = skipped 
	q_primitive3o7.txt = input file 	 q_primitive3o7_out.txt = output file 
	quadratic_out.txt is output file = skipped 
	quadratic_.txt = input file 	 quadratic__out.txt = output file 
	rabbit_out.txt is output file = skipped 
	rabbit.txt = input file 	 rabbit_out.txt = output file 
	spiral1_out.txt is output file = skipped 
	spiral1.txt = input file 	 spiral1_out.txt = output file 
	spiral6_out.txt is output file = skipped 
	spiral6.txt = input file 	 spiral6_out.txt = output file 
	spiral_marc_10_out.txt is output file = skipped 
	spiral_marc_10.txt = input file 	 spiral_marc_10_out.txt = output file 
	spiral_marc_out.txt is output file = skipped 
	spiral_marc.txt = input file 	 spiral_marc_out.txt = output file 
	square_spiral_out.txt is output file = skipped 
	square_spiral.txt = input file 	 square_spiral_out.txt = output file 
	thomasson85center_out.txt is output file = skipped 
	thomasson85center.txt = input file 	 thomasson85center_out.txt = output file 
	thomasson85_out.txt is output file = skipped 
	thomasson85.txt = input file 	 thomasson85_out.txt = output file 
	validationFailed_out2.txt = input file 	 validationFailed_out2_out.txt = output file 
	validationFailed_out.txt is output file = skipped 
	validationFailed.txt = input file 	 validationFailed_out.txt = output file 
end
```

# examples  
* see files in [txt directory](./txt):   
  * input files: *.txt
  * output files : *_out.txt
* [cubic Julia set - location by  Roger Lee Bagula](https://commons.wikimedia.org/wiki/File:Cubic_Julia_set_C_%3D-0.040000000000000036-0.78*I_with_internal_level_curves.png)
* [Critical orbit f(z) = z*z+-0.749413589136570+0.015312826507689*i](https://commons.wikimedia.org/wiki/File:Critical_orbit_f(z)_%3D_z*z%2Bc_and_c%3D-0.749413589136570%2B0.015312826507689*i.png)


# isues

## multiplicity

In case of multiple roots the numerical algorithm can show 2 roots close to each other.  See [output](./txt/cubic_p2_Bagula_out.txt) of the [cubic_p2_Bagula.txt](./txt/cubic_p2_Bagula.txt):

```txt
coefficients read from input file cubic_p2_Bagula.txt
	degree 3 coefficient = ( +1.0000000000000000 +0.0000000000000000*i) 
	degree 2 coefficient = ( +0.0000000000000000 +0.0000000000000000*i) 
	degree 1 coefficient = ( +0.0000000000000000 +0.0000000000000000*i) 
	degree 0 coefficient = ( -0.0400000000000000 -0.7800000000000000*i) 

Input polynomial p(z)=(1+0i)*z^3+(-0.040000000000000035527-0.78000000000000002665i)

derivative dp/dz = (3+0i)*z^2

2 critical points found

cp#0: -2.2351741790771484375e-08,-2.2351741790771484375e-08 . It's critical orbit is bounded and enters cycle #0 length=2 
      and it's stability = |multiplier|=0.95704 =attractive 
      internal angle = 0.045777099465243623055
      cycle = { 0.12845610612214106161,-0.42699144182978676643 ; -0.108141353107358687,-0.7232875185669475071 ; }

cp#1: -2.2351741790771484375e-08,9.4296410679817199707e-09 . It's critical orbit is bounded  and enters cycle #0
```

one can check using Maxima CAS that there is only one critical point z=0
```
%i1) f:z^3+c;
                                     3
(%o1)                               z  + c
(%i2) diff(f,z,1);
                                        2
(%o2)                                3 z
(%i3) solve(3*z^2=0);
(%o3)                               [z = 0]
```
 
## maxit

Some cases need a longer iteration. I recommend you change the following two lines in the code:

```cpp
int32_t maxit=25000;
```

to some other value - 50000 is working here.

If you want to use a higher maxit than 2^16 you also need to change the memory allocation constant:
```cpp
const int32_t MAXFIXORBITLEN=(1 << 16); // max length of a total orbit to some higher value than maxit.
```



# see also:
* [Periodic points of complex quadratic polynomial using Newton method in c](https://gitlab.com/adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method)
* [program describe ](https://en.wikibooks.org/wiki/Fractals/mandelbrot-numerics#m-describe) from mandelbrot-numerics  library by Claude Heiland-Allen
* [prm](https://github.com/raboehm/prm) - Find polynomial roots with multiplicities using mpmath by Bob Boehm

# credits
* basic stuff is the (modified) universal set of Newton starting points due to [Sutherland]
The initial points selected for the Newton iterations are based on a modification/simplification for practical purposes of the universal set of starting points by [Sutherland](http://pi.math.cornell.edu/~hubbard/NewtonInventiones.pdf). Sutherland et al use a set of circles that guarantee at least one converging starting point for every basin of attraction, the square used here (out of practical considerations) has no such guarantee. 
* [marcm200](https://github.com/marcm200) ( most of the code) 





# License
MIT License

Copyright (c) 2020 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE




# git

```git
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/Describe-iterated-map-.git
git push -u origin main
```

Subdirectory: 

```git
mkdir txt
git add *.txt
git mv  *.txt ./txt
git commit -m "move"
git push -u origin main
```
to overwrite

```git
git mv -f 
```

to remove file from remote repo

```git
git rm --cached file1.txt
git commit -m "remove"
git push -u origin main
```


local repo : ~Dokumenty/periodic

# shields

Shields.io is a service for concise, consistent, and legible badges in SVG and raster format, which can easily be included in GitHub readmes or any other web page.  

[github: badges/shields](https://github.com/badges/shields)

```
<p align="center">
  <a href="https://github.com/rust-fractal/rust-fractal-core/blob/master/LICENSE"><img src="https://img.shields.io/github/license/rust-fractal/rust-fractal-core" alt="Repository License"></a>
  <a href="https://github.com/rust-fractal/rust-fractal-core/"><img src="https://img.shields.io/tokei/lines/github/rust-fractal/rust-fractal-core" alt="Repository Size"></a>
  <a href="https://github.com/rust-fractal/rust-fractal-core/releases"><img src="https://img.shields.io/github/downloads/rust-fractal/rust-fractal-core/total?style=flat" alt="Github Release"></a>
</p>
```


