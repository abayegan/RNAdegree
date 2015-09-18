RNAdegree - Clotelab - Boston College
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
This program computes the expected network degree (expected number of neighbors) for the directed graph G = (V,E) , 
where V is the set of all secondary structures of a given RNA sequence, and E is the set of directed edges s → t, 
where structure t is obtained either from s by an element of move set MS1 consisting of base pair additions and removals 
or move set MS2 consisting of base pair additions, removals, and shifts. The MS2 algorithm is not a simple extension of 
the algorithm for MS1, but requires entirely new concepts and is surprisingly complex.
The set of programs are described in the following articles:

    Move set MS1: 

    1. Clote P. Expected degree of RNA secondary structure networks. J Comput Chem. 2015 Jan 15;36(2):103-17. doi: 10.1002/jcc.23776. Epub 2014 Nov 7

    Move set MS2: 
 
    2. P. Clote, A. Bayegan. Network properties of the ensemble of RNA structures. submitted, 25 June 2015. 

****Boltzmann, Model C****
In the Boltzmann model E(s) indicates the Turner energy of s, which involves free energy parameters for stacked base pairs, 
hairpins, bulges, internal loops and multiloops. The partition function Z = Σsexp(-E(s)/RT) can be computed by the McCaskill
algorithm, and the probability of structure s is the usual Boltzmann probability P(s)=exp(-E(s)/RT)/Z. The expected number of 
neighbors in this models can be computed by Q1,n/Z where Qi,j = Σ s∈ss[i,j]BF(s).N(s) , BF(s) denotes Boltzmann factor of s, 
and N(s) denotes the number of neighbors of s with respect to MS1 or MS2.
The code for the Boltzmann model is available in 'Boltzmann-src' folder.

****Uniform, Model B****
In model B, the energy of each structure is zero, so the partition function Z is the total number of structures of the given RNA, 
and the probability P(s) of each structure s is 1/Z. The expected number of neighbors in this models can be computed by Q1,n/Z 
where Qi,j = Σ s∈ss[i,j]N(s) , and N(s) denotes the number of neighbors of s with respect to MS1 or MS2. 
The code for the uniform model is available in 'Uniform-src' folder.

****Python source****
The folder 'python-src' includes a prototype for the C-implementation, in the uniform probability case, drivers for testing, 
wrapper to enumerate all structures, etc. and subsequently determine the expected number of neighbors by exhaustion. 
This code requires wrappers to access Vienna RNA Package.

****contact order****
The source code for computing contact order is included in 'contactOrder-src'.
