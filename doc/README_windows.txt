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
This product includes software developed by the OpenSSL Project for use in
the OpenSSL Toolkit (http://www.openssl.org/).  This product includes
cryptographic software written by Eric Young (eay@cryptsoft.com).
This product includes the GNU MP Library.^M

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
First backup wallet. Uninstall previous version and then follow setup
instructions. Double check balance after completing setup and starting up
client.

Setup
-----
Complete windows setup procedure and run Primecoin (Qt).

Website: http://primecoin.org
Forum: http://ppcointalk.org
Github (source code + sig + wiki): https://github.com/primecoin/primecoin
Sourceforge (release builds): https://sourceforge.net/projects/primecoin



Bitcoin 0.8.6 BETA
==================

Copyright (c) 2009-2013 Bitcoin Developers

Distributed under the MIT/X11 software license, see the accompanying
file COPYING or http://www.opensource.org/licenses/mit-license.php.
This product includes software developed by the OpenSSL Project for use in
the OpenSSL Toolkit (http://www.openssl.org/).  This product includes
cryptographic software written by Eric Young (eay@cryptsoft.com).


Intro
-----
Bitcoin is a free open source peer-to-peer electronic cash system that is
completely decentralized, without the need for a central server or trusted
parties.  Users hold the crypto keys to their own money and transact directly
with each other, with the help of a P2P network to check for double-spending.


Setup
-----
Unpack the files into a directory and run bitcoin-qt.exe.

Bitcoin-Qt is the original Bitcoin client and it builds the backbone of the network.
However, it downloads and stores the entire history of Bitcoin transactions;
depending on the speed of your computer and network connection, the synchronization
process can take anywhere from a few hours to a day or more.

See the bitcoin wiki at:
  https://en.bitcoin.it/wiki/Main_Page
for more help and information.
