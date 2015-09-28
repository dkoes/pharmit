/*
Pharmit
Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

*/
#include "BoundingBox.h"

ostream& operator<<(ostream& out, const BoundingBox& box)
{
	out << "(" << box.minx << "," << box.maxx << ") ";
	out << "(" << box.miny << "," << box.maxy << ") ";
	out << "(" << box.minz << "," << box.maxz << ")";
	return out;
}
