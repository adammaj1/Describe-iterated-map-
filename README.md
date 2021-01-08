# numerical periodicity detection

determine numerically the cycles of a polynomial Julia set

# algorithm

## find critical points
first find the critical points via classical Newton using a set of equidistant starting points at a square the size 3 times the Lagrangian estimate for the location of the roots of f'(z) (usually a total of 1024 points). 

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

# see also:
* [prm](https://github.com/raboehm/prm) - Find polynomial roots with multiplicities using mpmath by Bob Boehm

# contributors
* basic stuff is the (modified) universal set of Newton starting points due to Sutherland
* [marcm200](https://github.com/marcm200) ( most of the code) 





#License
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
