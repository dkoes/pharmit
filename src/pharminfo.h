/*
Pharmit
Copyright (C) 2015  David Ryan Koes and the University of Pittsburgh

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

/*
 * pharminfo.h
 *
 *  Created on: August 11, 2015
 *      Author: dkoes
 *
 *  This file contains routines for reading and writing slightly condensed pharmacophore
 *  features for a whole molecule.  It is used by shape search to filter by pharmacophores.
 */
#include "pharmarec.h"
#include <cstdio>
#include <vector>

//write points to file and return the starting offset
unsigned long writePharmacophoreInfo(FILE *f, const vector<PharmaPoint>& points, const Pharmas& pharmas);

//returns true if each query point matches at least one member of the pharmacophore
//pointed to by pharmacophore - avoid copying into memory
bool pharmacophoreMatchesQuery(const char *pharmacophore, const vector<PharmaPoint>& querypoints, const Pharmas& pharmas);
