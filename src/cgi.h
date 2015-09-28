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

/*
 * cgi.h
 *
 *  Created on: Mar 11, 2010
 *      Author: dkoes
 */

#ifndef PHARMITSERVER_CGI_H_
#define PHARMITSERVER_CGI_H_

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

//dump the cgi request for debugging
string cgiDump(Cgicc& cgi);

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

#endif /* PHARMITSERVER_CGI_H_ */
