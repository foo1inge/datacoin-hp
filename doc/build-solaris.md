Copyright (c) 2009-2013 Bitcoin Developers

Distributed under the MIT/X11 software license, see the accompanying
file COPYING or http://www.opensource.org/licenses/mit-license.php.
This product includes software developed by the OpenSSL Project for use in the [OpenSSL Toolkit](http://www.openssl.org/). This product includes
cryptographic software written by Eric Young ([eay@cryptsoft.com](mailto:eay@cryptsoft.com)), and UPnP software written by Thomas Bernard.

UNIX BUILD NOTES
====================

To Build
---------------------

	cd src/
	gmake -f makefile.solaris		# Headless bitcoin

See readme-qt.rst for instructions on building Bitcoin-Qt, the graphical user interface.

Dependencies
---------------------

 Library     Purpose           Description
 -------     -------           -----------
 libssl      SSL Support       Secure communications
 libdb       Berkeley DB       Blockchain & wallet storage
 libboost    Boost             C++ Library
 miniupnpc   UPnP Support      Optional firewall-jumping support

[miniupnpc](http://miniupnp.free.fr/) may be used for UPnP port mapping.  It can be downloaded from [here](
http://miniupnp.tuxfamily.org/files/).  UPnP support is compiled in and
turned off by default.  Set USE_UPNP to a different value to control this:

	USE_UPNP=     No UPnP support miniupnp not required
	USE_UPNP=0    (the default) UPnP support turned off by default at runtime
	USE_UPNP=1    UPnP support turned on by default at runtime

IPv6 support may be disabled by setting:

	USE_IPV6=0    Disable IPv6 support

Licenses of statically linked libraries:
 Berkeley DB   New BSD license with additional requirement that linked
               software must be free open source
 Boost         MIT-like license
 miniupnpc     New (3-clause) BSD license

- Versions used in this release:
-  GCC           4.5.2
-  OpenSSL       1.0.1c
-  Berkeley DB   6.0.20
-  Boost         1.54
-  miniupnpc     1.6

Building dependencies
---------------------

first set environment (64-bit)

export CXXFLAGS="-m64 -march=native -mtune=native -I/usr/local/include"
export CPPFLAGS="-m64 -march=native -mtune=native -I/usr/local/include"
export CFLAGS="-m64 -march=native -mtune=native -I/usr/local/include"
export LDFLAGS="-m64 -L/usr/local/lib -R/usr/local/lib -L/usr/gnu/lib/amd64 -R/usr/gnu/lib/amd64"

GMP
---

autoreconf

./configure --enable-cxx
make
make check
sudo make install

OpenSSL
-------

./Configure solaris64-x86_64-gcc --prefix=/usr/local --openssldir=/usr/local/openssl $CXXFLAGS $LDFLAGS enable-ec_nistp_64_gcc_128 enable-gmp enable-md2 enable-rc5 enable-rfc3779 zlib shared
gmake depend
gmake
gmake test
sudo gmake install

BerkeleyDB
----------

dist/s_config
cd build_unix
../dist/configure --prefix=/usr/local --enable-cxx --enable-pthread_api --enable-o_direct --enable-dtrace
make
sudo make install

Boost
-----

You need to edit boost/cstdint.hpp and add "|| defined(__sun__)" to defines at lines before "#include <inittypes.h>"

./bootstrap.sh --prefix=/usr --libdir=/usr/lib/amd64
sudo ./b2 variant=release link=shared threading=multi address-model=64 cxxflags="$CXXFLAGS" linkflags="$LDFLAGS" architecture=x86 instruction-set=native install

Security
--------
To help make your bitcoin installation more secure by making certain attacks impossible to
exploit even if a vulnerability is found, you can take the following measures:

* Position Independent Executable
    Build position independent code to take advantage of Address Space Layout Randomization
    offered by some kernels. An attacker who is able to cause execution of code at an arbitrary
    memory location is thwarted if he doesn't know where anything useful is located.
    The stack and heap are randomly located by default but this allows the code section to be
    randomly located as well.

    On an Amd64 processor where a library was not compiled with -fPIC, this will cause an error
    such as: "relocation R_X86_64_32 against `......' can not be used when making a shared object;"

    To build with PIE, use:

    	make -f makefile.unix ... -e PIE=1

    To test that you have built PIE executable, install scanelf, part of paxutils, and use:

    	scanelf -e ./bitcoin

    The output should contain:
     TYPE
    ET_DYN

* Non-executable Stack
    If the stack is executable then trivial stack based buffer overflow exploits are possible if
    vulnerable buffers are found. By default, bitcoin should be built with a non-executable stack
    but if one of the libraries it uses asks for an executable stack or someone makes a mistake
    and uses a compiler extension which requires an executable stack, it will silently build an
    executable without the non-executable stack protection.

    To verify that the stack is non-executable after compiling use:
    `scanelf -e ./bitcoin`

    the output should contain:
	STK/REL/PTL
	RW- R-- RW-

    The STK RW- means that the stack is readable and writeable but not executable.
