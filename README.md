<p align="center">
  <a href="https://github.com/adammaj1/Describe-iterated-map-"><img src="https://img.shields.io/tokei/lines/github/Describe-iterated-map-.git" alt="Repository Size"></a>
  <a href="https://github.com/adammaj1/Describe-iterated-map-/releases"><img src="https://img.shields.io/github/downloads/Describe-iterated-map-.git/total?style=flat" alt="Github Release"></a>
</p>


# numerical periodicity detection

determine numerically the cycles of a polynomial Julia set

# algorithm

## find critical points
first find [the critical points](https://en.wikipedia.org/wiki/Critical_point_(mathematics)) via [classical Newton method](https://en.wikibooks.org/wiki/Fractals/Mathematics/Newton_method) using a set of [equidistant starting points](https://gitlab.com/adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method) (  at a square the size 3 times [the Lagrangian estimate for the location of the roots of](https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds) f'(z) (usually a total of 1024 points). 

## compute the critical orbits
Then I construct the orbit of those using double precision (or sometimes float128), but without rounding control, for usually 25 000 iterations (escaping is tested using an escape radius either manually set or computed via the Douady estimate). 


## find periods
At the end of the orbit computation, I check backwards whether the last iterate has already been visited before (squared Euclidean distance of two points below 10^-15) to get the cycle length.



# use 
Program reads parameters from the input txt file and prints the output data (to the stdout ot the output txt file)


```c
g++ main.cpp -lm -Wall
./a.out rabbit.txt > rabbit_out.txt
```
# examples
* see files in [txt directory](./txt): 
  * input files: *.txt
  * output files : *_out.txt
* [cubic Julia set - location by  Roger Lee Bagula](https://commons.wikimedia.org/wiki/File:Cubic_Julia_set_C_%3D-0.040000000000000036-0.78*I_with_internal_level_curves.png)
* [Critical orbit f(z) = z*z+-0.749413589136570+0.015312826507689*i](https://commons.wikimedia.org/wiki/File:Critical_orbit_f(z)_%3D_z*z%2Bc_and_c%3D-0.749413589136570%2B0.015312826507689*i.png)

# see also:
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

local repo : ~Dokumenty/periodic

## shields

Shields.io is a service for concise, consistent, and legible badges in SVG and raster format, which can easily be included in GitHub readmes or any other web page.  

[github: badges/shields](https://github.com/badges/shields)

```
<p align="center">
  <a href="https://github.com/rust-fractal/rust-fractal-core/blob/master/LICENSE"><img src="https://img.shields.io/github/license/rust-fractal/rust-fractal-core" alt="Repository License"></a>
  <a href="https://github.com/rust-fractal/rust-fractal-core/"><img src="https://img.shields.io/tokei/lines/github/rust-fractal/rust-fractal-core" alt="Repository Size"></a>
  <a href="https://github.com/rust-fractal/rust-fractal-core/releases"><img src="https://img.shields.io/github/downloads/rust-fractal/rust-fractal-core/total?style=flat" alt="Github Release"></a>
</p>
```


