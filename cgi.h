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
 * cgi.h
 *
 *  Created on: Mar 11, 2010
 *      Author: dkoes
 */

#ifndef CGI_H_
#define CGI_H_

/*! \file FCgiIO.h
 * \brief Class that implements input and output through a FastCGI request.
 *
 * This class provides access to the input byte-stream and environment
 * variable interfaces of a FastCGI request.  It is fully compatible with the
 * Cgicc input API.
 *
 * It also provides access to the request's output and error streams, using a
 * similar interface.
 */

#include <ostream>
#include <string>
#include <map>

#include "fcgio.h"

#include "cgicc/CgiInput.h"
#include "cgicc/Cgicc.h"

using namespace std;

namespace cgicc
{

//dkoes - return true if there is data with name
bool cgiTagExists(Cgicc& cgi, const string& name);
//get an integer value, 0 is returned if there's a problem
long cgiGetInt(Cgicc& cgi, const string& name);

//empty string is returned if problems
string cgiGetString(Cgicc& cgi, const string& name);

//get a double value, 0 is returned if there's a problem
double cgiGetDouble(Cgicc& cgi, const string& name);

/*
 * This class provides access to the input byte-stream and environment
 * variable interfaces of a FastCGI request.  It is fully compatible with the
 * Cgicc input API.
 *
 * It also provides access to the request's output and error streams, using a
 * similar interface.
 */
class FastCgiIO : public cgicc::CgiInput, public std::ostream
{
public:

	FastCgiIO(FCGX_Request& request);
	FastCgiIO(const FastCgiIO& io);
	virtual inline ~FastCgiIO()
	{}

	virtual inline size_t read(char *data, size_t length)
	{
		return FCGX_GetStr(data, length, fRequest.in);
	}

	virtual inline std::string getenv(const char *varName)
	{
		return fEnv[varName];
	}

	inline std::ostream& err(void)
	{
		return fErr;
	}
	//@}

protected:
	FCGX_Request& fRequest;
	fcgi_streambuf fOutBuf;
	fcgi_streambuf fErrBuf;
	std::ostream fErr;
	std::map<std::string, std::string> fEnv;
};

} // namespace cgicc

#endif /* CGI_H_ */
