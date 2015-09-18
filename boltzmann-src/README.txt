RNAdegree - Clotelab - Boston College
using Turner 2004/1999 energy model with parameters from Vienna RNA Package

/******************************************************************************
 *   Copyright (C) 2015  Amir Bayegan, Peter Clote                            *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/

Program computes the expected network degree (expected number of
neighbors) for the directed graph G = (V,E) , where V is the set of all
secondary structures of a given RNA sequence, and E is the set of directed
edges s → t, where structure t is obtained from s by an element of move set
MS2 consisting of base pair additions, removals, and shifts, The algorithm is
not a simple extension of the algorithm for MS1, but requires entirely new
concepts and is surprisingly complex.

USAGE: ./RNAdegree rnaSequence [-ms1|-ms2] [-t temperature] [-e turner99|turner04|andronescu07] [-v] [-h]
where 
	rnaSequence: any sequence of upper or lowercase A,C,G,U characters
	-ms1 : moving set allowing addition and removal of base pairs.
	-ms2 : moving set allowing addition, removal and shift of base pairs. 
	-t : temperature in celsius at which the given RNA is folded. Default is 37.
	-e : select the energy parameters applied to fold the RNA sequences from Turner 1999/2004 or Andronescu 2007 . Default is turner04.
	-v : also prints partition function and the total number of structures in the output.
	-h : displays this information.  

Output:
  two floating point numbers(if -v is not set), where first is expected degree and second
  is expected degree divided by sequence length.


