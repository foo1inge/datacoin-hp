Primecoin High Performance Version 12
=====================================

This is the high performance version of Sunny King's Primecoin tree.

Features:
 * Extended sieve algorithm
 * Cache-friendly segmented sieve
 * Adjustable sieve parameters
 * Mining threads use GMP for faster bignum arithmetic

Donations are welcome if you want to support my work.

BTC: 1EaHwHBWeoJtSM2jEdx9Su1NcKvdXbsqxX
LTC: LPD1zDChmqcqKGHFHuLX2JWMMEC5jD5J4j
XPM: AJHjbkVzHhHugd5bpKDtddVDfhtEB8jQZ4

Sieve parameters
----------------

 * -sievesize determines how many numbers go through the sieve in one go. The
default value is 917504 (default value of -l1cachesize times 32). This value
has very little impact on performance.

 * -sievefilterprimes determines how many primes factors will be filtered out
by the sieve. The default value is 7849. This value has a minor impact on
performance.

 * -sieveextensions determines the number of sieve "extensions". Each extension
will effectively process an additional set of numbers through the sieve
relatively cheaply. The default value is 9. This value has a minor impact on
performance.

 * -l1cachesize determines the sieve segment size (in bytes). This should be
the size of the L1 or L2 cache (or some other number close to them). The
default value is 28672. This value has a minor impact on mining performance.

 * -primorial determines the primorial that is used as a multiplier in all
candidate numbers. This should a small prime number such as 47 or 53. Use the
parameters -debug -printprimorial to see which primorials are used by default.
This value has a minor impact on mining performance.

Some parameters can be changed on the fly using the following RPC commands:
 * setsievefilterprimes
 * setsievesize
 * setsieveextensions

Old and removed parameters:
 * -sieveroundpercentage
 * -sievepercentage

Primecoin 0.1.2 BETA
====================

Copyright (c) 2013 Primecoin Developers

Distributed under conditional MIT/X11 software license, see the accompanying
file COPYING.
This product includes software developed by the OpenSSL Project for use in the [OpenSSL Toolkit](http://www.openssl.org/). This product includes
cryptographic software written by Eric Young ([eay@cryptsoft.com](mailto:eay@cryptsoft.com)), and UPnP software written by Thomas Bernard.
This product includes the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

Intro
---------------------
Primecoin is a free open source cryptocurrency that implements the first
scientific computing proof-of-work for cryptocurrencies. The unique
proof-of-work design searches for rare prime formations, providing
experimental value for mathematicians to further understand the nature and
distribution related to prime number, a simple yet mysterious construct of
arithmetic that continues to baffle the top minds of mankind.

Upgrade
--------------------
First backup wallet. Then follow setup instructions. Double check balance
after completing setup and starting up client.

Setup
--------------------
You need the Qt4 run-time libraries to run Primecoin-Qt. On Debian or Ubuntu:
        `sudo apt-get install libqtgui4`

Unpack the files into a directory and run:

- bin/32/primecoin-qt (GUI, 32-bit)
- bin/32/primecoind (headless, 32-bit)
- bin/64/primecoin-qt (GUI, 64-bit)
- bin/64/primecoind (headless, 64-bit)

Website: http://primecoin.org
Forum: http://ppcointalk.org
Github (source code + sig + wiki): https://github.com/primecoin/primecoin
Sourceforge (release builds): https://sourceforge.net/projects/primecoin


Bitcoin 0.8.6 BETA
====================

Copyright (c) 2009-2013 Bitcoin Developers

Distributed under the MIT/X11 software license, see the accompanying
file COPYING or http://www.opensource.org/licenses/mit-license.php.
This product includes software developed by the OpenSSL Project for use in the [OpenSSL Toolkit](http://www.openssl.org/). This product includes
cryptographic software written by Eric Young ([eay@cryptsoft.com](mailto:eay@cryptsoft.com)), and UPnP software written by Thomas Bernard.


Intro
---------------------
Bitcoin is a free open source peer-to-peer electronic cash system that is
completely decentralized, without the need for a central server or trusted
parties.  Users hold the crypto keys to their own money and transact directly
with each other, with the help of a P2P network to check for double-spending.


Setup
---------------------
You need the Qt4 run-time libraries to run Bitcoin-Qt. On Debian or Ubuntu:
	`sudo apt-get install libqtgui4`

Unpack the files into a directory and run:

- bin/32/bitcoin-qt (GUI, 32-bit)
- bin/32/bitcoind (headless, 32-bit)
- bin/64/bitcoin-qt (GUI, 64-bit)
- bin/64/bitcoind (headless, 64-bit)

See the documentation at the [Bitcoin Wiki](https://en.bitcoin.it/wiki/Main_Page)
for help and more information.


Other Pages
---------------------
- [Unix Build Notes](build-unix.md)
- [OSX Build Notes](build-osx.md)
- [Windows Build Notes](build-msw.md)
- [Coding Guidelines](coding.md)
- [Release Process](release-process.md)
- [Release Notes](release-notes.md)
- [Multiwallet Qt Development](multiwallet-qt.md)
- [Unit Tests](unit-tests.md)
- [Translation Process](translation_process.md)
