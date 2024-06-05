Pharmit
====

Copyright (c) David Ryan Koes, University of Pittsburgh and contributors.
All rights reserved.

Pharmit is licensed under both the BSD 3-clause license and the GNU
Public License version 2. Any use of the code that retains its reliance
on the GPL-licensed OpenBabel library is subject to the terms of the GPL2.

Use of the Pharmit code independently of OpenBabel (or any other
GPL2 licensed software) may choose between the BSD or GPL licenses.

See the LICENSE file provided with the distribution for more information.

Pharmit is available for searching and creating libraries at [http://pharmit.csb.pitt.edu](http://pharmit.csb.pitt.edu).  It generally isn't necessary to build it from source.

BUILDING
====

CMake is required to build Pharmit.  Starting from the src directory:
```
mkdir build
cd build
#most likely you will have to specify the location of smina nad lemon
cmake .. -DSMINA_DIR=$HOME/git/smina -DLEMON_DIR=/usr/lib/cmake/
make -j12
```
## UBUNTU 22.04

```bash
sudo apt install  git autoconf automake libtool ghostscript liblemon-dev libeigen3-dev libann-dev bmagic libcgicc-dev libgoogle-perftools-dev libglpk-dev coinor-* libjsoncpp-dev cmake libboost-dev swig libxml2-dev libcairo2-dev libboost-all-dev libcurl4-openssl-dev

git clone https://github.com/FastCGI-Archives/fcgi2.git
cd fcgi2
autoreconf -i
./configure
sudo make
sudo make install
cd ..

git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build
cd build
sudo cmake .. -DPYTHON_BINDINGS=1 -DRUN_SWIG=1  -DWITH_MAEPARSER=0 -DWITH_COORDGEN=0
sudo make -j12
sudo make install
cd ../..

git clone http://git.code.sf.net/p/smina/code smina
cd smina
mkdir build
cd build
sudo cmake ..
sudo make -j8
sudo cp libsmina.a /usr/local/lib/libsmina.a
cd ../..

git clone http://git.code.sf.net/p/pharmit/code pharmit
cd pharmit/src
mkdir build
cd build
#change `$path/smina` the smina directory to reflect the install location from above
sudo cmake .. -DSMINALIB=/usr/local/lib/libsmina.a -DSMINA_DIR=$path/smina/ -DLEMON_DIR=/usr/lib/x86_64-linux-gnu/cmake/lemon/
sudo make -j12
```

## UBUNTU 18.04

```
apt install  git autoconf automake libtool ghostscript liblemon-dev libeigen3-dev libann-dev bmagic libcgicc-dev libgoogle-perftools-dev libglpk-dev coinor-* libjsoncpp-dev cmake libboost-dev swig python-dev libxml2-dev  libcairo2-dev libboost-all-dev libcurl4-openssl-dev

#NOTE: Currently the liblemon-dev package has an incorrect cmake config file that does
#not set LEMON_LIBRARY to the correct location.  Either edit /usr/lib/cmake/LEMONConfig.cmake
#so the correct location (/usr/lib/x86_64-linux-gnu/) is set or set LEMON_LIBRARY directly
#on the cmake commandline

#fastcgi
git clone https://github.com/FastCGI-Archives/fcgi2.git
cd fcgi2
autoreconf -i
./configure
make
make install
cd ..

#openbabel 
#you can _probably_ use the libopenbabel-dev package, but it doesn't install cmake files by default, so you'd have to provide a custom FindOpenBabel module
#In general, I recommend using openbabel 3.0 or later as it fixed a number of bugs.
git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build
cd build
#presumably at some point maeparser and coordgen won't be broken...
cmake .. -DPYTHON_BINDINGS=1 -DRUN_SWIG=1  -DWITH_MAEPARSER=0 -DWITH_COORDGEN=0
make -j12
make install
cd ../..

#smina
git clone http://git.code.sf.net/p/smina/code smina
cd smina
mkdir build
cd build
cmake ..
make -j8

cd smina/build/linux/release/
make -j8
cp libsmina.a /usr/local/lib/
cd ../../../..

#pharmit - ensure SMINA_DIR is the install location above
git clone http://git.code.sf.net/p/pharmit/code pharmit
cd pharmit/src
mkdir build
cd build
#/usr/lib/cmake/LEMONConfig.cmake is editted to point LEMON_LIBRARY to /usr/lib/x86_64-linux-gnu/
cmake .. -DSMINALIB=/usr/local/lib/libsmina.a -DSMINA_DIR=$HOME/git/smina -DLEMON_DIR=/usr/lib/cmake/ 
make -j12
```

## CENTOS (these instructions may be out of date)
```
yum groupinstall 'Development Tools'
yum install git boost-devel autoconf automake libtool make cmake gperftools-libs gperftools-devel ghostscript libcurl-devel

yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm
yum install glpk-devel eigen3-devel fcgi  coin-or-*

#fastcgi
git clone https://github.com/FastCGI-Archives/fcgi2.git
cd fcgi2
autoreconf -i
./configure
make
make install
cd ..

#jsoncpp
git clone https://github.com/open-source-parsers/jsoncpp.git
cd jsoncpp
#using the older version mostly to avoid having to upgrade cmake wget
git checkout 1.7.7
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS=1 ..
make
make install
cd ../..

#cgicc
wget http://ftp.gnu.org/gnu/cgicc/cgicc-3.2.19.tar.gz
tar xvfz cgicc-3.2.19.tar.gz
cd cgicc-3.2.19
./configure
make 
make install
cd ..

#openbabel
#I'm getting the development version since it fixes some serious bugs with aromaticity detection,
#but the version in the epel repository is probably fine
git clone https://github.com/openbabel/openbabel.git
cd openbabel
mkdir build
cd build
cmake ..
make
make install
cd ../..

#eigen
wget http://bitbucket.org/eigen/eigen/get/3.3.5.tar.gz
tar xvfz 3.3.5.tar.g
cd eigen-eigen-b3f3d4950030/
mkdir build
cd build
cmake ..
make install
cd ../..

#lemon
wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.1.tar.gz
tar xvfz lemon-1.3.1.tar.gz 
cd lemon-1.3.1
mkdir build
#remove incompatible cmake policy
sed -i '/POLICY/d' CMakeLists.txt 
cd build
cmake ..
make
make install
cd ../..

#libann
wget https://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz
tar xvfz ann_1.1.2.tar.gz
cd ann_1.1.2
make linux-g++
cp -r include/ANN /usr/include/
cp lib/libANN.a /usr/lib/libann.a
cd ..

#bitmagic
git clone https://github.com/tlk00/BitMagic.git 
cd BitMagic
#ubuntu version is actualy 3.7
git checkout v3.8.0
#include as part of pharmit - can also install if you wish, but that's not what cmake is expecting
cp -r src/ ~/pharmit/src/bm

#smina
git clone http://git.code.sf.net/p/smina/code smina
cd smina/build/linux/release/
make
cp libsmina.a /usr/local/lib/
 
#pharmit
git clone http://git.code.sf.net/p/pharmit/code pharmit
cd pharmit/src
mkdir build
cmake .. -DCMAKE_CXX_FLAGS=-std=c++0x -DSMINA_DIR=$HOME/smina
#lots of warnings about missing directory which can be ignored for now
```

USING
====

The --help option will provide the following:
```
USAGE: pharmitserver [options] --cmd command [pharma, dbcreate, dbcreateserverdir, dbsearch, server]

OPTIONS:
  -cmd=<string>        - command [pharma, dbcreate, dbcreateserverdir, dbsearch, server]
  -dbdir=<string>      - database directory(s)
  -dbinfo=<string>     - [dbcreateserverdir] JSON file describing database subset
  -extra-info          - Output additional molecular properties.  Slower.
  -file-partition      - Partion database slices based on files
  -help                - Display available options (--help-hidden for more)
  -in=<string>         - input file(s)
  -ligs=<string>       - [dbcreateserverdir] Text file listing locations of molecules
  -logdir=<string>     - log directory for server
  -max-hits=<n>        - return at most n results
  -max-nrot=<uint>     - maximum allowed rotatable bonds
  -max-orient=<n>      - return at most n orientations of each conformation
  -max-rmsd=<number>   - maximum allowed RMSD; default max allowed by query
  -max-weight=<uint>   - maximum allowed molecular weight
  -min-nrot=<uint>     - minimum allowed rotatable bonds
  -min-port=<uint>     - port for minimization server
  -min-server=<string> - minimization server address
  -min-weight=<uint>   - minimum allowed molecular weight
  -noindex             - [dbcreateserverdir] Do not create indices
  -nthreads=<n>        - utilize n threads; default 1
  -out=<string>        - output file(s)
  -pharmaspec=<string> - pharmacophore specification
  -port=<uint>         - port for server to listen on
  -prefixes=<string>   - [dbcreateserverdir,server] File of directory prefixes to use for striping.
  -print               - print results
  -q                   - quiet; suppress informational messages
  -receptor=<string>   - Receptor file for interaction pharmacophroes
  -reduceconfs=<n>     - return at most n conformations for each molecule
  -show-query          - print query points
  -singledir=<string>  - Specify a single directory to recreate on a dbcreateserverdir command
  -sort-rmsd           - Sort results by RMSD.
  -timestamp=<string>  - Specify timestamp to use for server dirs (for use with single)
  -unweighted-rmsd     - Compute minimal RMSD without radius weights
```

To IDENTIFY the pharmacophore features of a molecule ex.pdb:
```
 pharmitserver pharma -in ex.pdb
 pharmitserver pharma -in ex.pdb -out out.sdf
 pharmitserver pharma -in ex.pdb -out out.json
```
For interactive, graphical editting of pharmacophore features, try
[http://pharmit.csb.pitt.edu](http://pharmit.csb.pitt.edu)
The saved session file is a json file that can be used as input for Pharmer.

To CREATE a database DB from library.sdf:
```
 pharmitserver dbcreate -dbdir DB -in library.sdf
```
If you have multiple disk drives, you can improve performance by striping the
database across the drives:
```
 pharmitserver dbcreate -dbdir /drive1/DB -dbdir /drive2/DB -in library.sdf
```
If you pre-split the input file things will go faster:
```
 pharmitserver dbcreate -dbdir /drive1/DB -dbdir /drive2/DB -file-partition -in library1.sdf -in library2.sdf
```

To SEARCH a database:
```
 pharmitserver dbsearch -dbdir DB -in query.json
 pharmitserver dbsearch -dbdir DB -in query.sdf
 pharmitserver dbsearch -dbdir DB -in query.ph4
```
When setting up a SERVER you create separate directories for each database to search:
```
 pharmitserver dbcreateserverdir  -ligs molport.ligs -prefixes ../dbprefixes -dbinfo dbinfo.json 
```
The ligs files provides the files to add with unique ids and names (see extractsmisubset.py).  e.g.:
```
/data22/conformers/541/5416265.sdf.gz 5416265 MolPort-020-216-564 MCULE-1493214818 PubChem-56905170 ZINC000072151660
/data17/conformers/53/534584.sdf.gz 534584 MolPort-000-887-730 MolPort-035-708-781 MCULE-2607511133 MCULE-2857980910 PubChem-19618431 PubChem-90484256 ZINC
000002535004 ZINC02535004
/data08/conformers/18890/188901458.sdf.gz 188901458 MolPort-007-641-322 MCULE-8587443370
```

The prefixes file specifies the directories (on different drives) to partition the database across. e.g.:
```
/data00/databases
/data01/databases
/data02/databases
```

The dbinfo.json file describes the library. e.g.:
```
{
"name" : "MolPort",
"html" : "MolPort<span><a target=\"_blank\" href=\"https://www.molport.com\" class=\"ui-icon-info ui-icon\"></a></span>",
"subdir" : "molport",
"prefix" : "MolPort"
}
```

To run a SERVER:
```
pharmitserver server -port 16000 -prefixes /home/dkoes/vendors/dbprefixes  -logdir /home/dkoes/log -min-server localhost -min-port 18000
```

To run a MINIMIZATION server:
```
~/git/smina/build/linux/release/server -port 18000 -logfile /home/dkoes/log/minlog
```

To HOST a server, add a fastcgi rule connecting the webserver hosting /web to the storage server running the backend.  e.g. (apache):
```
ln -s ~/pharmit/web /var/www/html
mkdir /var/www/html/fcgi-bin

apt install apache2 libapache2-mod-php php-mysql apache2-dev

git clone https://github.com/FastCGI-Archives/mod_fastcgi.git
cp Makefile.AP2 Makefile
#change /usr/local/apache2 to /usr/share/apache2 in Makefile
make
make install

#create /etc/apache2/mods-available/fastcgi.load
LoadModule fastcgi_module /usr/lib/apache2/modules/mod_fastcgi.so
#fastcgi.conf


<IfModule mod_fastcgi.c>
  AddHandler fastcgi-script .fcgi
  #FastCgiWrapper /usr/lib/apache2/suexec
  FastCgiIpcDir /var/lib/apache2/fastcgi

Alias /fcgi-bin/ /var/www/html/fcgi-bin/

<Directory /var/www/html/fcgi-bin/>
 SetHandler fastcgi-script
 Options +ExecCGI
</Directory>

FastCgiExternalServer /var/www/html/fcgi-bin/pharmitserv.fcgi -host 127.0.0.1:16000 -idle-timeout 300
FastCgiExternalServer /var/www/html/fcgi-bin/createlib.fcgi -host 127.0.0.1:11111 -idle-timeout 300

FastCgiConfig -autoUpdate -maxClassProcesses 1

</IfModule>

# enable
a2enmod fastcgi
mkdir /var/www/html/fcgi-bin  #this directory has to exist
```

LIBRARY CREATION
====

The library creation scripts require rdkit to be installed.
```
apt install python3-numpy  mysql-server mysql-client python3-mysqldb python3-pip
pip3 install flup psutil
mysql < ~/pharmit/scripts/pharmit.sql
mysql < ~/pharmit/scripts/conformers.sql
#also, run the mysql commands:
 CREATE USER 'pharmit'@'localhost';
 GRANT ALL PRIVILEGES ON pharmit.* TO 'pharmit'@'localhost';
 GRANT ALL PRIVILEGES ON conformers.* TO 'pharmit'@'localhost';

 INSERT INTO pharmit.users (email,name,maxprivatedbs,maxprivateconfs,maxdbs,maxconfs) VALUES ("guest","Anonymous",0,0,0,10000);

git clone https://github.com/rdkit/rdkit.git
cd rdkit
mkdir build
cd build
cmake .. -DPYTHON_EXECUTABLE=/usr/bin/python3
make -j12
```

BUGS (Not Really)
====
When working with large datasets spread across multiple hard drives,
you will likely run into a number of limits that will result in 
segmentation faults (type 6 or 11).  The following changes need to be in place
before building the database, or it will silently create incomplete 
databases that crash when you search them.

Make sure nofile is set high enough (1000000) in `/etc/security/limits.conf`

Make sure vm.max_map_count is set high enough - add the following to `/etc/sysctl.conf`
```
vm.max_map_count=1000000
```

Make sure you build against libcurl4-openssl-dev, not gnutls, or you
will get fortify_fail errors.

Continuing with curl bugs, DNS timeouts can cause a longjmp error:
http://stackoverflow.com/questions/9191668/error-longjmp-causes-uninitialized-stack-frame
to resolve these rebuild curl with --enable-ares

You can ignore all the curl issues entirely if you build with -DSKIP_REGISTERZINC
which is probably what you want unless you are using the full ZINC database.

Some versions of eigen3 will trigger an assert, OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG, at compile time.
This can be ignored in the latest versions of eigen3 with -DEIGEN_DISABLE_STACK_SIZE_ASSERT
If your version hasn't had this patch applied (which is currently true with Ubuntu 14.04),
then you will have to modify Eigen/src/Core/DenseStorage.h appropriately.

Currently, on Ubuntu 14.04, g++ 4.8 experiences an internal compiler error.
Use g++-4.6 to get around this.

BUGS (Really)
====

Send any bug reports to dkoes@pitt.edu with a complete test case.

CITING
====
Please use the citation from http://pubs.acs.org/doi/abs/10.1021/ci200097m

Default Pharmacophore Definitions
---
````
Aromatic 18 0 1.1 1 0.1 
a1aaaaa1
a1aaaa1

HydrogenDonor 1 1 0.5 1 0.1 
[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]
[#8!H0&!$([OH][C,S,P]=O)]
[#16!H0]

HydrogenAcceptor 89 2 0.5 1 0.1 
[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]
[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]

PositiveIon 7 3 0.75 0 0.1 
[+,+2,+3,+4]
[$(CC)](=N)N
[$(C(N)(N)=N)]
[$(n1cc[nH]c1)]

NegativeIon 8 4 0.75 0 0.1 
[-,-2,-3,-4]
C(=O)[O-,OH,OX1]
[$([S,P](=O)[O-,OH,OX1])]
c1[nH1]nnn1
c1nn[nH1]n1
C(=O)N[OH1,O-,OX1]
C(=O)N[OH1,O-]
CO(=N[OH1,O-])
[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]

Hydrophobic 6 5 1 0 2 
a1aaaaa1
a1aaaa1
[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]
[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]
*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]
[C&r3]1~[C&r3]~[C&r3]1
[C&r4]1~[C&r4]~[C&r4]~[C&r4]1
[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1
[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1
[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1
[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1
[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]
[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]
[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]
[$([S]~[#6])&!$(S~[!#6])]
````
