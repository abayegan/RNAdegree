contactOrder - Clotelab - Boston College

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

The 3D absolute contact order for an RNA structure is defined by sum_{i<j}(j − i)/N , where the sum is 
over all N pairs of residues i, j that are in contact, taken here to mean that residues i, j each contain 
a heavy atom (non-hydrogen) within 6̊A, and that i, j are not consecutive (j!=i + 1). 3D contact order is 
calculated after parsing a pdb file and selecting all the atoms in contact(<6 angstrom).

The pseudoknot (pknot) absolute contact order is defined as sum_{i<j} (j − i)/N , where the sum is over all 
N base pairs (i, j) determined by RNAview, a program that determines hydrogen-bonded atoms of distinct nucleotides 
in a PDB file of RNA and additionally classifies the base pair with respect to the Leontis-Westhof classification. 

The 2D absolute contact order is defined as sum_{i<j} (j − i)/N , where the sum is over all N RNA contact order for 
the secondary and tertiary structures are computed. The secondary structure contact order is computed for a secondary 
structure obtainded from maximizing the number of base pairs predicted by RNAVIEW from pdb files. 

Input: a file containing paths to the pdbfiles  of interest. Each line of the file should include exactly one path.

Output: 1) A folder named RNAVIEW_out containing the output files of rnaview.
		2) Two folders named 3DcontactList and compatible3DcontactList. In these folders for every pdb, a file with extension '_3DcontactList'.
		or '_viewCompat3DContactList' is created indicating which atoms of the molecule are in contact with respect to the defenitions above.
		3) Three files with extensions '_2DcontactOrder' '_3DcontactOrder' and 'viewCompat3DcontactOrder' containing the final contact order values.


