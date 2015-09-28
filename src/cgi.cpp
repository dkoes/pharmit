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
 * cgi.cpp
 *
 *  Created on: Mar 11, 2010
 *      Author: dkoes
 */

#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <sstream>
#include "cgi.h"

using namespace std;

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

string cgiDump(Cgicc& cgi)
{
	stringstream ret;
	ret << "cgiquery " << cgi.getEnvironment().getQueryString() << "\n"
	<< "cgipost " << cgi.getEnvironment().getPostData() << "\n";
	return ret.str();
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
