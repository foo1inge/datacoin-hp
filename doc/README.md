Primecoin High Performance Version
==================================

This is the high performance version of Sunny King's Primecoin tree.

Features:
 * Use GMP for bignum calculations in the mining threads
 * Replaced some bignum calculations with 64-bit arithmetic inside the sieve
 * Reduced the amount of memory allocations
 * L1 and L2 cache optimizations
 * Process only 10% of base primes when weaving the sieve
 * Configurable sieve size

Donations are welcome if you want to support my work.

BTC: 1EaHwHBWeoJtSM2jEdx9Su1NcKvdXbsqxX
LTC: LPD1zDChmqcqKGHFHuLX2JWMMEC5jD5J4j
XPM: AJHjbkVzHhHugd5bpKDtddVDfhtEB8jQZ4

Primecoin 0.1.1 BETA
====================

Copyright (c) 2013 Primecoin Developers

Distributed under conditional MIT/X11 software license, see the accompanying
file COPYING.
This product includes software developed by the OpenSSL Project for use in the [OpenSSL Toolkit](http://www.openssl.org/). This product includes
cryptographic software written by Eric Young ([eay@cryptsoft.com](mailto:eay@cryptsoft.com)), and UPnP software written by Thomas Bernard.

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



Bitcoin 0.8.3 BETA
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
