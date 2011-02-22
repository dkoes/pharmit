/*
Pharmer: Efficient and Exact 3D Pharmacophore Search
Copyright (C) 2011  David Ryan Koes and the University of Pittsburgh

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
 * cgi.cpp
 *
 *  Created on: Mar 11, 2010
 *      Author: dkoes
 */

#include <iostream>
#include <stdexcept>
#include <cstdlib>

#include "cgi.h"

namespace cgicc {
//dkoes - return true if there is data with name
bool cgiTagExists(Cgicc& cgi, const string& name)
{
	form_iterator itr = cgi[name];
	return itr != cgi.getElements().end();
}

//get an integer value, 0 is returned if there's a problem
long cgiGetInt(Cgicc& cgi, const string& name)
{
	if (cgiTagExists(cgi, name))
		return cgi[name]->getIntegerValue();
	else
		return 0;
}

//empty string is returned if problems
string cgiGetString(Cgicc& cgi, const string& name)
{
	if (cgiTagExists(cgi, name))
		return cgi[name]->getValue();
	else
		return "";
}

//get a double value, 0 is returned if there's a problem
double cgiGetDouble(Cgicc& cgi, const string& name)
{
	if (cgiTagExists(cgi, name))
		return cgi[name]->getDoubleValue();
	else
		return 0;
}
}
cgicc::FastCgiIO::FastCgiIO(FCGX_Request& request) :
	std::ostream(&fOutBuf), fRequest(request), fOutBuf(request.out), fErrBuf(
			request.err), fErr(&fErrBuf)
{
	rdbuf(&fOutBuf);
	fErr.rdbuf(&fErrBuf);

	// Parse environment
	for (char **e = fRequest.envp; *e != NULL; ++e)
	{
		std::string s(*e);
		std::string::size_type i = s.find('=');
		if (i == std::string::npos)
			throw std::runtime_error("Illegally formed environment");
		fEnv[s.substr(0, i)] = s.substr(i + 1);
	}
}

cgicc::FastCgiIO::FastCgiIO(const FastCgiIO& io) :
	CgiInput(io), std::ostream(&fOutBuf), fRequest(io.fRequest),
			fErr(&fErrBuf), fEnv(io.fEnv)
{
	rdbuf(&fOutBuf);
	fErr.rdbuf(&fErrBuf);
}
