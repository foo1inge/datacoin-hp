// Copyright (c) 2013 Primecoin developers
// Distributed under conditional MIT/X11 software license,
// see the accompanying file COPYING

#include "prime.h"
#include <climits>

/**********************/
/* PRIMECOIN PROTOCOL */
/**********************/

// Prime Table
std::vector<unsigned int> vPrimes;
unsigned int nSieveSize = nDefaultSieveSize;
unsigned int nSieveFilterPrimes = nDefaultSieveFilterPrimes;
unsigned int nSieveExtensions = nDefaultSieveExtensions;
unsigned int nL1CacheSize = nDefaultL1CacheSize;

static unsigned int int_invert(unsigned int a, unsigned int nPrime);

void GeneratePrimeTable()
{
    const unsigned int nDefaultSieveExt = (fTestNet) ? nDefaultSieveExtensionsTestnet : nDefaultSieveExtensions;
    nSieveExtensions = (unsigned int)GetArg("-sieveextensions", nDefaultSieveExt);
    nSieveExtensions = std::max(std::min(nSieveExtensions, nMaxSieveExtensions), nMinSieveExtensions);
    nSieveSize = (unsigned int)GetArg("-sievesize", nDefaultSieveSize);
    nSieveSize = std::max(std::min(nSieveSize, nMaxSieveSize), nMinSieveSize);
    nSieveFilterPrimes = (unsigned int)GetArg("-sievefilterprimes", nDefaultSieveFilterPrimes);
    nSieveFilterPrimes = std::max(std::min(nSieveFilterPrimes, nMaxSieveFilterPrimes), nMinSieveFilterPrimes);
    nL1CacheSize = (unsigned int)GetArg("-l1cachesize", nDefaultL1CacheSize);
    nL1CacheSize = std::max(std::min(nL1CacheSize, nMaxL1CacheSize), nMinL1CacheSize);
    nL1CacheSize = nL1CacheSize / 8 * 8; // make it a multiple of 8
    printf("GeneratePrimeTable() : setting nSieveExtensions = %u, nSieveSize = %u, nSieveFilterPrimes = %u, nL1CacheSize = %u\n", nSieveExtensions, nSieveSize, nSieveFilterPrimes, nL1CacheSize);

    const unsigned nPrimeTableLimit = 1000000u;
    vPrimes.clear();
    // Generate prime table using sieve of Eratosthenes
    std::vector<bool> vfComposite (nPrimeTableLimit, false);
    for (unsigned int nFactor = 2; nFactor * nFactor < nPrimeTableLimit; nFactor++)
    {
        if (vfComposite[nFactor])
            continue;
        for (unsigned int nComposite = nFactor * nFactor; nComposite < nPrimeTableLimit; nComposite += nFactor)
            vfComposite[nComposite] = true;
    }
    for (unsigned int n = 2; n < nPrimeTableLimit; n++)
        if (!vfComposite[n])
            vPrimes.push_back(n);
    printf("GeneratePrimeTable() : prime table [1, %u] generated with %u primes\n", nPrimeTableLimit, (unsigned int) vPrimes.size());
}

// Mining statistics
uint64 nTotalTests;
unsigned int nTotalBlocksFound;
std::vector<uint64> vTotalChainsFound;
boost::timer::cpu_timer minerTimer;
int nSieveTargetLength = -1;

// Primecoin HP: Optional automatic donations with every block found
CBitcoinAddress donationAddress;
double dDonationPercentage;

void ResetMinerStatistics()
{
    nTotalTests = 0;
    nTotalBlocksFound = 0;
    vTotalChainsFound = std::vector<uint64> (nMaxChainLength, 0);
}

void InitPrimeMiner()
{
    ResetMinerStatistics();
    nSieveTargetLength = std::min((int)GetArg("-sievetargetlength", nDefaultSieveTargetLength), (int)nMaxChainLength);
    if (nSieveTargetLength > 0)
        printf("InitPrimeMiner() : Setting sieve target length to %d\n", nSieveTargetLength);

    // Primecoin HP: Optional automatic donations with every block found
    std::string strDonationPercentage = GetArg("-donationpercentage", "0.0");
    std::string strDonationAddress = GetArg("-donationaddress", !fTestNet ? strDefaultDonationAddress : strDefaultDonationAddressTestnet);
    dDonationPercentage = atof(strDonationPercentage.c_str());
    if (dDonationPercentage < dMinDonationPercentage)
        dDonationPercentage = 0.0;
    dDonationPercentage = std::min(dDonationPercentage, dMaxDonationPercentage);
    donationAddress = CBitcoinAddress(strDonationAddress);
    if (!donationAddress.IsValid())
    {
        dDonationPercentage = 0.0;
        printf("InitPrimeMiner(): Donation address is invalid, disabling donations\n");
    }
    if (dDonationPercentage > 0.001)
        printf("InitPrimeMiner(): Donating %2.2f%% of every block found to %s (thank you!)\n", dDonationPercentage, strDonationAddress.c_str());
    else
        printf("InitPrimeMiner(): Donations disabled\n");
}

void PrintMinerStatistics()
{
    printf("========================================================================\n");
    printf("Miner statistics\n");
    printf("========================================================================\n");

    boost::timer::cpu_times const elapsed_times(minerTimer.elapsed());
    int64 nRunningTime = elapsed_times.wall;
    double dRunningHours = (double)nRunningTime / 3600000000000.0;
    int64 nCPUTime = elapsed_times.system + elapsed_times.user;
    double dCPUHours = (double)nCPUTime / 3600000000000.0;
    printf("Running time: %.4f hours\n", dRunningHours);
    printf("CPU time: %.4f hours\n", dCPUHours);

    printf("Tests: %"PRI64u"\n", nTotalTests);
    printf("Blocks found: %u\n", nTotalBlocksFound);

    // Find the last non-zero chain count
    unsigned int nMaxPrintLength = nMaxChainLength;
    for (int i = nMaxChainLength - 1; i >= 0; i--)
    {
        if (vTotalChainsFound[i] > 0)
        {
            nMaxPrintLength = i + 1;
            break;
        }
    }

    printf("\n");
    printf("Chain statistics\n");
    for (unsigned int i = 0; i < nMaxPrintLength; i++)
        printf("%u-chains: %"PRI64u"\n", i + 1, vTotalChainsFound[i]);

    printf("========================================================================\n");

    // Reset statistics
    nHPSTimerStart = 0;
    ResetMinerStatistics();
}

void PrintCompactStatistics(volatile unsigned int vFoundChainCounter[nMaxChainLength])
{
    std::string strOutput;
    if (fLogTimestamps)
        strOutput = "chainstats ";
    else
        strOutput = strprintf("%s chainstats ", DateTimeStrFormat("%Y-%m-%d %H:%M:%S", GetTime()).c_str());
    for (unsigned int i = 0; i < nMaxChainLength; i++)
    {
        if (vFoundChainCounter[i])
            strOutput += strprintf(" %uch: %u", i + 1, vFoundChainCounter[i]);
    }
    printf("%s\n", strOutput.c_str());

    // Reset the statistics
    for (unsigned int i = 0; i < nMaxChainLength; i++)
        vFoundChainCounter[i] = 0;
}

// Get next prime number of p
bool PrimeTableGetNextPrime(unsigned int& p)
{
    BOOST_FOREACH(unsigned int nPrime, vPrimes)
    {
        if (nPrime > p)
        {
            p = nPrime;
            return true;
        }
    }
    return false;
}

// Get previous prime number of p
bool PrimeTableGetPreviousPrime(unsigned int& p)
{
    unsigned int nPrevPrime = 0;
    BOOST_FOREACH(unsigned int nPrime, vPrimes)
    {
        if (nPrime >= p)
            break;
        nPrevPrime = nPrime;
    }
    if (nPrevPrime)
    {
        p = nPrevPrime;
        return true;
    }
    return false;
}

// Compute Primorial number p#
void Primorial(unsigned int p, mpz_class& mpzPrimorial)
{
    unsigned long nPrimorial = 1;
    unsigned int i;
    if (sizeof(unsigned long) >= 8)
    {
        // Fast 64-bit loop if possible
        for (i = 0; i < 15; i++)
        {
            unsigned int nPrime = vPrimes[i];
            if (nPrime > p) break;
            nPrimorial *= nPrime;
        }
    }
    else
    {
        // Fast 32-bit loop first
        for (i = 0; i < 9; i++)
        {
            unsigned int nPrime = vPrimes[i];
            if (nPrime > p) break;
            nPrimorial *= nPrime;
        }
    }

    mpzPrimorial = nPrimorial;
    for (; i < vPrimes.size(); i++)
    {
        unsigned int nPrime = vPrimes[i];
        if (nPrime > p) break;
        mpzPrimorial *= nPrime;
    }
}

// Compute Primorial number p#
// Fast 32-bit version assuming that p <= 23
unsigned int PrimorialFast(unsigned int p)
{
    unsigned int nPrimorial = 1;
    BOOST_FOREACH(unsigned int nPrime, vPrimes)
    {
        if (nPrime > p) break;
        nPrimorial *= nPrime;
    }
    return nPrimorial;
}

// Compute first primorial number greater than or equal to pn
void PrimorialAt(mpz_class& bn, mpz_class& mpzPrimorial)
{
    mpzPrimorial = 1;
    BOOST_FOREACH(unsigned int nPrime, vPrimes)
    {
        mpzPrimorial *= nPrime;
        if (mpzPrimorial >= bn)
            return;
    }
}

// Proof-of-work Target (prime chain target):
//   format - 32 bit, 8 length bits, 24 fractional length bits

unsigned int nTargetInitialLength = 7; // initial chain length target
unsigned int nTargetMinLength = 6;     // minimum chain length target

unsigned int TargetGetLimit()
{
    return (nTargetMinLength << nFractionalBits);
}

unsigned int TargetGetInitial()
{
    return (nTargetInitialLength << nFractionalBits);
}

unsigned int TargetGetLength(unsigned int nBits)
{
    return ((nBits & TARGET_LENGTH_MASK) >> nFractionalBits);
}

bool TargetSetLength(unsigned int nLength, unsigned int& nBits)
{
    if (nLength >= 0xff)
        return error("TargetSetLength() : invalid length=%u", nLength);
    nBits &= TARGET_FRACTIONAL_MASK;
    nBits |= (nLength << nFractionalBits);
    return true;
}

static void TargetIncrementLength(unsigned int& nBits)
{
    nBits += (1 << nFractionalBits);
}

static void TargetDecrementLength(unsigned int& nBits)
{
    if (TargetGetLength(nBits) > nTargetMinLength)
        nBits -= (1 << nFractionalBits);
}

unsigned int TargetGetFractional(unsigned int nBits)
{
    return (nBits & TARGET_FRACTIONAL_MASK);
}

uint64 TargetGetFractionalDifficulty(unsigned int nBits)
{
    return (nFractionalDifficultyMax / (uint64) ((1llu<<nFractionalBits) - TargetGetFractional(nBits)));
}

bool TargetSetFractionalDifficulty(uint64 nFractionalDifficulty, unsigned int& nBits)
{
    if (nFractionalDifficulty < nFractionalDifficultyMin)
        return error("TargetSetFractionalDifficulty() : difficulty below min");
    uint64 nFractional = nFractionalDifficultyMax / nFractionalDifficulty;
    if (nFractional > (1u<<nFractionalBits))
        return error("TargetSetFractionalDifficulty() : fractional overflow: nFractionalDifficulty=%016"PRI64x, nFractionalDifficulty);
    nFractional = (1u<<nFractionalBits) - nFractional;
    nBits &= TARGET_LENGTH_MASK;
    nBits |= (unsigned int)nFractional;
    return true;
}

std::string TargetToString(unsigned int nBits)
{
    return strprintf("%02x.%06x", TargetGetLength(nBits), TargetGetFractional(nBits));
}

unsigned int TargetFromInt(unsigned int nLength)
{
    return (nLength << nFractionalBits);
}

// Get mint value from target
// Primecoin mint rate is determined by target
//   mint = 999 / (target length ** 2)
// Inflation is controlled via Moore's Law
bool TargetGetMint(unsigned int nBits, uint64& nMint)
{
    nMint = 0;
    static uint64 nMintLimit = 999llu * COIN;
    CBigNum bnMint = nMintLimit;
    if (TargetGetLength(nBits) < nTargetMinLength)
        return error("TargetGetMint() : length below minimum required, nBits=%08x", nBits);
    bnMint = (bnMint << nFractionalBits) / nBits;
    bnMint = (bnMint << nFractionalBits) / nBits;
    bnMint = (bnMint / CENT) * CENT;  // mint value rounded to cent
    nMint = bnMint.getuint256().Get64();
    if (nMint > nMintLimit)
    {
        nMint = 0;
        return error("TargetGetMint() : mint value over limit, nBits=%08x", nBits);
    }
    return true;
}

// Get next target value
bool TargetGetNext(unsigned int nBits, int64 nInterval, int64 nTargetSpacing, int64 nActualSpacing, unsigned int& nBitsNext)
{
    nBitsNext = nBits;
    // Convert length into fractional difficulty
    uint64 nFractionalDifficulty = TargetGetFractionalDifficulty(nBits);
    // Compute new difficulty via exponential moving toward target spacing
    CBigNum bnFractionalDifficulty = nFractionalDifficulty;
    bnFractionalDifficulty *= ((nInterval + 1) * nTargetSpacing);
    bnFractionalDifficulty /= ((nInterval - 1) * nTargetSpacing + nActualSpacing + nActualSpacing);
    if (bnFractionalDifficulty > nFractionalDifficultyMax)
        bnFractionalDifficulty = nFractionalDifficultyMax;
    if (bnFractionalDifficulty < nFractionalDifficultyMin)
        bnFractionalDifficulty = nFractionalDifficultyMin;
    uint64 nFractionalDifficultyNew = bnFractionalDifficulty.getuint256().Get64();
    if (fDebug && GetBoolArg("-printtarget"))
        printf("TargetGetNext() : nActualSpacing=%d nFractionDiff=%016"PRI64x" nFractionDiffNew=%016"PRI64x"\n", (int)nActualSpacing, nFractionalDifficulty, nFractionalDifficultyNew);
    // Step up length if fractional past threshold
    if (nFractionalDifficultyNew > nFractionalDifficultyThreshold)
    {
        nFractionalDifficultyNew = nFractionalDifficultyMin;
        TargetIncrementLength(nBitsNext);
    }
    // Step down length if fractional at minimum
    else if (nFractionalDifficultyNew == nFractionalDifficultyMin && TargetGetLength(nBitsNext) > nTargetMinLength)
    {
        nFractionalDifficultyNew = nFractionalDifficultyThreshold;
        TargetDecrementLength(nBitsNext);
    }
    // Convert fractional difficulty back to length
    if (!TargetSetFractionalDifficulty(nFractionalDifficultyNew, nBitsNext))
        return error("TargetGetNext() : unable to set fractional difficulty prev=0x%016"PRI64x" new=0x%016"PRI64x, nFractionalDifficulty, nFractionalDifficultyNew);
    return true;
}

// Check Fermat probable primality test (2-PRP): 2 ** (n-1) = 1 (mod n)
// true: n is probable prime
// false: n is composite; set fractional length in the nLength output
static bool FermatProbablePrimalityTest(const CBigNum& n, unsigned int& nLength)
{
    CAutoBN_CTX pctx;
    CBigNum a = 2; // base; Fermat witness
    CBigNum e = n - 1;
    CBigNum r;
    BN_mod_exp(&r, &a, &e, &n, pctx);
    if (r == 1)
        return true;
    // Failed Fermat test, calculate fractional length
    unsigned int nFractionalLength = (((n-r) << nFractionalBits) / n).getuint();
    if (nFractionalLength >= (1 << nFractionalBits))
        return error("FermatProbablePrimalityTest() : fractional assert");
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// Test probable primality of n = 2p +/- 1 based on Euler, Lagrange and Lifchitz
// fSophieGermain:
//   true:  n = 2p+1, p prime, aka Cunningham Chain of first kind
//   false: n = 2p-1, p prime, aka Cunningham Chain of second kind
// Return values
//   true: n is probable prime
//   false: n is composite; set fractional length in the nLength output
static bool EulerLagrangeLifchitzPrimalityTest(const CBigNum& n, bool fSophieGermain, unsigned int& nLength)
{
    CAutoBN_CTX pctx;
    CBigNum a = 2;
    CBigNum e = (n - 1) >> 1;
    CBigNum r;
    BN_mod_exp(&r, &a, &e, &n, pctx);
    CBigNum nMod8 = n % 8;
    bool fPassedTest = false;
    if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
        fPassedTest = (r == 1);
    else if (fSophieGermain && (nMod8 == 3)) // Lifchitz
        fPassedTest = ((r+1) == n);
    else if ((!fSophieGermain) && (nMod8 == 5)) // Lifchitz
        fPassedTest = ((r+1) == n);
    else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
        fPassedTest = (r == 1);
    else
        return error("EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = %d, %s", nMod8.getint(), (fSophieGermain? "first kind" : "second kind"));

    if (fPassedTest)
        return true;
    // Failed test, calculate fractional length
    r = (r * r) % n; // derive Fermat test remainder
    unsigned int nFractionalLength = (((n-r) << nFractionalBits) / n).getuint();
    if (nFractionalLength >= (1 << nFractionalBits))
        return error("EulerLagrangeLifchitzPrimalityTest() : fractional assert");
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength)
{
    return strprintf("%s%s", (nChainType==PRIME_CHAIN_CUNNINGHAM1)? "1CC" : ((nChainType==PRIME_CHAIN_CUNNINGHAM2)? "2CC" : "TWN"), TargetToString(nChainLength).c_str());
}

// primorial form of prime chain origin
std::string GetPrimeOriginPrimorialForm(CBigNum& bnPrimeChainOrigin)
{
    CBigNum bnNonPrimorialFactor = bnPrimeChainOrigin;
    unsigned int nPrimeSeq = 0;
    while (nPrimeSeq < vPrimes.size() && bnNonPrimorialFactor % vPrimes[nPrimeSeq] == 0)
    {
        bnNonPrimorialFactor /= vPrimes[nPrimeSeq];
        nPrimeSeq++;
    }
    return strprintf("%s*%u#", bnNonPrimorialFactor.ToString().c_str(), (nPrimeSeq > 0)? vPrimes[nPrimeSeq-1] : 0);
}

// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
// Return value:
//   true - Probable Cunningham Chain found (length at least 2)
//   false - Not Cunningham Chain
static bool ProbableCunninghamChainTest(const CBigNum& n, bool fSophieGermain, bool fFermatTest, unsigned int& nProbableChainLength)
{
    nProbableChainLength = 0;
    CBigNum N = n;

    // Fermat test for n first
    if (!FermatProbablePrimalityTest(N, nProbableChainLength))
        return false;

    // Euler-Lagrange-Lifchitz test for the following numbers in chain
    while (true)
    {
        TargetIncrementLength(nProbableChainLength);
        N = N + N + (fSophieGermain? 1 : (-1));
        if (fFermatTest)
        {
            if (!FermatProbablePrimalityTest(N, nProbableChainLength))
                break;
        }
        else
        {
            if (!EulerLagrangeLifchitzPrimalityTest(N, fSophieGermain, nProbableChainLength))
                break;
        }
    }

    return (TargetGetLength(nProbableChainLength) >= 2);
}

// Test probable prime chain for: nOrigin
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
bool ProbablePrimeChainTest(const CBigNum& bnPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int& nChainLengthCunningham1, unsigned int& nChainLengthCunningham2, unsigned int& nChainLengthBiTwin)
{
    nChainLengthCunningham1 = 0;
    nChainLengthCunningham2 = 0;
    nChainLengthBiTwin = 0;

    // Test for Cunningham Chain of first kind
    ProbableCunninghamChainTest(bnPrimeChainOrigin-1, true, fFermatTest, nChainLengthCunningham1);
    // Test for Cunningham Chain of second kind
    ProbableCunninghamChainTest(bnPrimeChainOrigin+1, false, fFermatTest, nChainLengthCunningham2);
    // Figure out BiTwin Chain length
    // BiTwin Chain allows a single prime at the end for odd length chain
    nChainLengthBiTwin =
        (TargetGetLength(nChainLengthCunningham1) > TargetGetLength(nChainLengthCunningham2))?
            (nChainLengthCunningham2 + TargetFromInt(TargetGetLength(nChainLengthCunningham2)+1)) :
            (nChainLengthCunningham1 + TargetFromInt(TargetGetLength(nChainLengthCunningham1)));

    return (nChainLengthCunningham1 >= nBits || nChainLengthCunningham2 >= nBits || nChainLengthBiTwin >= nBits);
}

// Check prime proof-of-work
bool CheckPrimeProofOfWork(uint256 hashBlockHeader, unsigned int nBits, const CBigNum& bnPrimeChainMultiplier, unsigned int& nChainType, unsigned int& nChainLength)
{
    // Check target
    if (TargetGetLength(nBits) < nTargetMinLength || TargetGetLength(nBits) > 99)
        return error("CheckPrimeProofOfWork() : invalid chain length target %s", TargetToString(nBits).c_str());

    // Check header hash limit
    if (hashBlockHeader < hashBlockHeaderLimit)
        return error("CheckPrimeProofOfWork() : block header hash under limit");
    // Check target for prime proof-of-work
    CBigNum bnPrimeChainOrigin = CBigNum(hashBlockHeader) * bnPrimeChainMultiplier;
    if (bnPrimeChainOrigin < bnPrimeMin)
        return error("CheckPrimeProofOfWork() : prime too small");
    // First prime in chain must not exceed cap
    if (bnPrimeChainOrigin > bnPrimeMax)
        return error("CheckPrimeProofOfWork() : prime too big");

    // Check prime chain
    unsigned int nChainLengthCunningham1 = 0;
    unsigned int nChainLengthCunningham2 = 0;
    unsigned int nChainLengthBiTwin = 0;
    if (!ProbablePrimeChainTest(bnPrimeChainOrigin, nBits, false, nChainLengthCunningham1, nChainLengthCunningham2, nChainLengthBiTwin))
        return error("CheckPrimeProofOfWork() : failed prime chain test target=%s length=(%s %s %s)", TargetToString(nBits).c_str(),
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str());
    if (nChainLengthCunningham1 < nBits && nChainLengthCunningham2 < nBits && nChainLengthBiTwin < nBits)
        return error("CheckPrimeProofOfWork() : prime chain length assert target=%s length=(%s %s %s)", TargetToString(nBits).c_str(),
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str());

    // Double check prime chain with Fermat tests only
    unsigned int nChainLengthCunningham1FermatTest = 0;
    unsigned int nChainLengthCunningham2FermatTest = 0;
    unsigned int nChainLengthBiTwinFermatTest = 0;
    if (!ProbablePrimeChainTest(bnPrimeChainOrigin, nBits, true, nChainLengthCunningham1FermatTest, nChainLengthCunningham2FermatTest, nChainLengthBiTwinFermatTest))
        return error("CheckPrimeProofOfWork() : failed Fermat test target=%s length=(%s %s %s) lengthFermat=(%s %s %s)", TargetToString(nBits).c_str(),
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str(),
            TargetToString(nChainLengthCunningham1FermatTest).c_str(), TargetToString(nChainLengthCunningham2FermatTest).c_str(), TargetToString(nChainLengthBiTwinFermatTest).c_str());
    if (nChainLengthCunningham1 != nChainLengthCunningham1FermatTest ||
        nChainLengthCunningham2 != nChainLengthCunningham2FermatTest ||
        nChainLengthBiTwin != nChainLengthBiTwinFermatTest)
        return error("CheckPrimeProofOfWork() : failed Fermat-only double check target=%s length=(%s %s %s) lengthFermat=(%s %s %s)", TargetToString(nBits).c_str(), 
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str(),
            TargetToString(nChainLengthCunningham1FermatTest).c_str(), TargetToString(nChainLengthCunningham2FermatTest).c_str(), TargetToString(nChainLengthBiTwinFermatTest).c_str());

    // Select the longest primechain from the three chain types
    nChainLength = nChainLengthCunningham1;
    nChainType = PRIME_CHAIN_CUNNINGHAM1;
    if (nChainLengthCunningham2 > nChainLength)
    {
        nChainLength = nChainLengthCunningham2;
        nChainType = PRIME_CHAIN_CUNNINGHAM2;
    }
    if (nChainLengthBiTwin > nChainLength)
    {
        nChainLength = nChainLengthBiTwin;
        nChainType = PRIME_CHAIN_BI_TWIN;
    }

    // Check that the certificate (bnPrimeChainMultiplier) is normalized
    if (bnPrimeChainMultiplier % 2 == 0 && bnPrimeChainOrigin % 4 == 0)
    {
        unsigned int nChainLengthCunningham1Extended = 0;
        unsigned int nChainLengthCunningham2Extended = 0;
        unsigned int nChainLengthBiTwinExtended = 0;
        if (ProbablePrimeChainTest(bnPrimeChainOrigin / 2, nBits, false, nChainLengthCunningham1Extended, nChainLengthCunningham2Extended, nChainLengthBiTwinExtended))
        { // try extending down the primechain with a halved multiplier
            if (nChainLengthCunningham1Extended > nChainLength || nChainLengthCunningham2Extended > nChainLength || nChainLengthBiTwinExtended > nChainLength)
                return error("CheckPrimeProofOfWork() : prime certificate not normalzied target=%s length=(%s %s %s) extend=(%s %s %s)",
                    TargetToString(nBits).c_str(),
                    TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str(),
                    TargetToString(nChainLengthCunningham1Extended).c_str(), TargetToString(nChainLengthCunningham2Extended).c_str(), TargetToString(nChainLengthBiTwinExtended).c_str());
        }
    }

    return true;
}

// prime target difficulty value for visualization
double GetPrimeDifficulty(unsigned int nBits)
{
    return ((double) nBits / (double) (1 << nFractionalBits));
}

// Estimate work transition target to longer prime chain
unsigned int EstimateWorkTransition(unsigned int nPrevWorkTransition, unsigned int nBits, unsigned int nChainLength)
{
    int64 nInterval = 500;
    int64 nWorkTransition = nPrevWorkTransition;
    unsigned int nBitsCeiling = 0;
    TargetSetLength(TargetGetLength(nBits)+1, nBitsCeiling);
    unsigned int nBitsFloor = 0;
    TargetSetLength(TargetGetLength(nBits), nBitsFloor);
    uint64 nFractionalDifficulty = TargetGetFractionalDifficulty(nBits);
    bool fLonger = (TargetGetLength(nChainLength) > TargetGetLength(nBits));
    if (fLonger)
        nWorkTransition = (nWorkTransition * (((nInterval - 1) * nFractionalDifficulty) >> 32) + 2 * ((uint64) nBitsFloor)) / ((((nInterval - 1) * nFractionalDifficulty) >> 32) + 2);
    else
        nWorkTransition = ((nInterval - 1) * nWorkTransition + 2 * ((uint64) nBitsCeiling)) / (nInterval + 1);
    return nWorkTransition;
}

/********************/
/* PRIMECOIN MINING */
/********************/

// Check Fermat probable primality test (2-PRP): 2 ** (n-1) = 1 (mod n)
// true: n is probable prime
// false: n is composite; set fractional length in the nLength output
static bool FermatProbablePrimalityTestFast(const mpz_class& n, unsigned int& nLength, CPrimalityTestParams& testParams, bool fFastFail = false)
{
    // Faster GMP version
    mpz_t& mpzE = testParams.mpzE;
    mpz_t& mpzR = testParams.mpzR;

    mpz_sub_ui(mpzE, n.get_mpz_t(), 1);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, n.get_mpz_t());
    if (mpz_cmp_ui(mpzR, 1) == 0)
        return true;
    if (fFastFail)
        return false;
    // Failed Fermat test, calculate fractional length
    mpz_sub(mpzE, n.get_mpz_t(), mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, n.get_mpz_t());
    unsigned int nFractionalLength = mpz_get_ui(mpzE);
    if (nFractionalLength >= (1 << nFractionalBits))
        return error("FermatProbablePrimalityTest() : fractional assert");
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// Test probable primality of n = 2p +/- 1 based on Euler, Lagrange and Lifchitz
// fSophieGermain:
//   true:  n = 2p+1, p prime, aka Cunningham Chain of first kind
//   false: n = 2p-1, p prime, aka Cunningham Chain of second kind
// Return values
//   true: n is probable prime
//   false: n is composite; set fractional length in the nLength output
static bool EulerLagrangeLifchitzPrimalityTestFast(const mpz_class& n, bool fSophieGermain, unsigned int& nLength, CPrimalityTestParams& testParams, bool fFastFail = false)
{
    // Faster GMP version
    mpz_t& mpzE = testParams.mpzE;
    mpz_t& mpzR = testParams.mpzR;
    mpz_t& mpzRplusOne = testParams.mpzRplusOne;

    mpz_sub_ui(mpzE, n.get_mpz_t(), 1);
    mpz_tdiv_q_2exp(mpzE, mpzE, 1);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, n.get_mpz_t());
    unsigned int nMod8 = mpz_get_ui(n.get_mpz_t()) % 8;
    bool fPassedTest = false;
    if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    else if (fSophieGermain && (nMod8 == 3)) // Lifchitz
    {
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, n.get_mpz_t());
    }
    else if ((!fSophieGermain) && (nMod8 == 5)) // Lifchitz
    {
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, n.get_mpz_t());
    }
    else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    else
        return error("EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = %d, %s", nMod8, (fSophieGermain? "first kind" : "second kind"));

    if (fPassedTest)
        return true;
    if (fFastFail)
        return false;

    // Failed test, calculate fractional length
    mpz_mul(mpzE, mpzR, mpzR);
    mpz_tdiv_r(mpzR, mpzE, n.get_mpz_t()); // derive Fermat test remainder

    mpz_sub(mpzE, n.get_mpz_t(), mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, n.get_mpz_t());
    unsigned int nFractionalLength = mpz_get_ui(mpzE);
    
    if (nFractionalLength >= (1 << nFractionalBits))
        return error("EulerLagrangeLifchitzPrimalityTest() : fractional assert");
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
}

// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
static void ProbableCunninghamChainTestFast(const mpz_class& n, bool fSophieGermain, unsigned int& nProbableChainLength, CPrimalityTestParams& testParams)
{
    nProbableChainLength = 0;

    // Fermat test for n first
    if (!FermatProbablePrimalityTestFast(n, nProbableChainLength, testParams, true))
        return;

    // Euler-Lagrange-Lifchitz test for the following numbers in chain
    mpz_class &N = testParams.mpzN;
    N = n;
    for (unsigned int nChainSeq = 1; true; nChainSeq++)
    {
        TargetIncrementLength(nProbableChainLength);
        N <<= 1;
        N += (fSophieGermain? 1 : (-1));
        bool fFastFail = nChainSeq < 4;
        if (!EulerLagrangeLifchitzPrimalityTestFast(N, fSophieGermain, nProbableChainLength, testParams, fFastFail))
            break;
    }
}

// Test Probable BiTwin Chain for: mpzOrigin
// Test the numbers in the optimal order for any given chain length
// Gives the correct length of a BiTwin chain even for short chains
static void ProbableBiTwinChainTestFast(const mpz_class& mpzOrigin, unsigned int& nProbableChainLength, CPrimalityTestParams& testParams)
{
    mpz_class& mpzOriginMinusOne = testParams.mpzOriginMinusOne;
    mpz_class& mpzOriginPlusOne = testParams.mpzOriginPlusOne;
    nProbableChainLength = 0;

    // Fermat test for origin-1 first
    mpzOriginMinusOne = mpzOrigin - 1;
    if (!FermatProbablePrimalityTestFast(mpzOriginMinusOne, nProbableChainLength, testParams, true))
        return;
    TargetIncrementLength(nProbableChainLength);

    // Fermat test for origin+1
    mpzOriginPlusOne = mpzOrigin + 1;
    if (!FermatProbablePrimalityTestFast(mpzOriginPlusOne, nProbableChainLength, testParams, true))
        return;
    TargetIncrementLength(nProbableChainLength);

    // Euler-Lagrange-Lifchitz test for the following numbers in chain
    for (unsigned int nChainSeq = 2; true; nChainSeq += 2)
    {
        mpzOriginMinusOne <<= 1;
        mpzOriginMinusOne++;
        bool fFastFail = nChainSeq < 4;
        if (!EulerLagrangeLifchitzPrimalityTestFast(mpzOriginMinusOne, true, nProbableChainLength, testParams, fFastFail))
            break;
        TargetIncrementLength(nProbableChainLength);

        mpzOriginPlusOne <<= 1;
        mpzOriginPlusOne--;
        if (!EulerLagrangeLifchitzPrimalityTestFast(mpzOriginPlusOne, false, nProbableChainLength, testParams, fFastFail))
            break;
        TargetIncrementLength(nProbableChainLength);
    }
}

// Test probable prime chain for: nOrigin
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
static bool ProbablePrimeChainTestFast(const mpz_class& mpzPrimeChainOrigin, CPrimalityTestParams& testParams)
{
    const unsigned int nBits = testParams.nBits;
    const unsigned int nCandidateType = testParams.nCandidateType;
    mpz_class& mpzOriginMinusOne = testParams.mpzOriginMinusOne;
    mpz_class& mpzOriginPlusOne = testParams.mpzOriginPlusOne;
    unsigned int& nChainLength = testParams.nChainLength;
    nChainLength = 0;

    // Test for Cunningham Chain of first kind
    if (nCandidateType == PRIME_CHAIN_CUNNINGHAM1)
    {
        mpzOriginMinusOne = mpzPrimeChainOrigin - 1;
        ProbableCunninghamChainTestFast(mpzOriginMinusOne, true, nChainLength, testParams);
    }
    else if (nCandidateType == PRIME_CHAIN_CUNNINGHAM2)
    {
        // Test for Cunningham Chain of second kind
        mpzOriginPlusOne = mpzPrimeChainOrigin + 1;
        ProbableCunninghamChainTestFast(mpzOriginPlusOne, false, nChainLength, testParams);
    }
    else if (nCandidateType == PRIME_CHAIN_BI_TWIN)
    {
        ProbableBiTwinChainTestFast(mpzPrimeChainOrigin, nChainLength, testParams);
    }

    return (nChainLength >= nBits);
}

// Perform Fermat test with trial division
// Return values:
//   true  - passes trial division test and Fermat test; probable prime
//   false - failed either trial division or Fermat test; composite
bool ProbablePrimalityTestWithTrialDivision(const mpz_class& mpzCandidate, unsigned int nTrialDivisionLimit, CPrimalityTestParams& testParams)
{
    unsigned int nDivisor = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23;
    unsigned int nDivisorPrimes = 9;

    // Fast trial division for the first few primes
    unsigned long nModulo = mpz_tdiv_ui(mpzCandidate.get_mpz_t(), nDivisor);
    for (unsigned int i = 0; i < nDivisorPrimes; i++)
    {
        if (nModulo % vPrimes[i] == 0)
            return false;
    }

    // Trial division
    for (unsigned int i = nDivisorPrimes; i < nTrialDivisionLimit; i++)
    {
        if (mpz_divisible_ui_p(mpzCandidate.get_mpz_t(), vPrimes[i]))
            return false;
    }
    unsigned int nLength = 0;
    return (FermatProbablePrimalityTestFast(mpzCandidate, nLength, testParams, true));
}

static void SieveDebugChecks(unsigned int nBits, unsigned int nTriedMultiplier, unsigned int nCandidateType, mpz_class& mpzHash, mpz_class& mpzFixedMultiplier, mpz_class& mpzChainOrigin)
{
    // Debugging code to verify the sieve output
    const unsigned int nTargetLength = TargetGetLength(nBits);
    mpz_class mpzChainN;
    mpz_class mpzChainNMod;
    if (nCandidateType == PRIME_CHAIN_CUNNINGHAM1 || nCandidateType == PRIME_CHAIN_BI_TWIN)
    {
        unsigned int nCC1Length = nTargetLength;
        if (nCandidateType == PRIME_CHAIN_BI_TWIN)
            nCC1Length = (nTargetLength + 1) / 2;
        for (unsigned int nChainPosition = 0; nChainPosition < nCC1Length; nChainPosition++)
        {
            mpzChainN = mpzChainOrigin << nChainPosition;
            mpzChainN--;
            for (unsigned int nPrimeSeq = 0; nPrimeSeq < nSieveFilterPrimes; nPrimeSeq++)
            {
                if (mpz_divisible_ui_p(mpzChainN.get_mpz_t(), vPrimes[nPrimeSeq]) > 0)
                {
                    std::string strHash = mpzHash.get_str();
                    std::string strFixedMultiplier = mpzFixedMultiplier.get_str();
                    printf("SIEVE BUG: %s * %s * %u * 2^%u - 1 is divisible by %u!\n", strHash.c_str(), strFixedMultiplier.c_str(), nTriedMultiplier, nChainPosition, vPrimes[nPrimeSeq]);
                }
            }
        }
    }
    if (nCandidateType == PRIME_CHAIN_CUNNINGHAM2 || nCandidateType == PRIME_CHAIN_BI_TWIN)
    {
        unsigned int nCC2Length = nTargetLength;
        if (nCandidateType == PRIME_CHAIN_BI_TWIN)
            nCC2Length = nTargetLength / 2;
        for (unsigned int nChainPosition = 0; nChainPosition < nCC2Length; nChainPosition++)
        {
            mpzChainN = mpzChainOrigin << nChainPosition;
            mpzChainN++;
            for (unsigned int nPrimeSeq = 0; nPrimeSeq < nSieveFilterPrimes; nPrimeSeq++)
            {
                if (mpz_divisible_ui_p(mpzChainN.get_mpz_t(), vPrimes[nPrimeSeq]) > 0)
                {
                    std::string strHash = mpzHash.get_str();
                    std::string strFixedMultiplier = mpzFixedMultiplier.get_str();
                    printf("SIEVE BUG: %s * %s * %u * 2^%u + 1 is divisible by %u!\n", strHash.c_str(), strFixedMultiplier.c_str(), nTriedMultiplier, nChainPosition, vPrimes[nPrimeSeq]);
                }
            }
        }
    }
    if (nCandidateType == 0)
    {
        std::string strHash = mpzHash.get_str();
        std::string strFixedMultiplier = mpzFixedMultiplier.get_str();
        printf("SIEVE BUG: %s * %s * %u has unknown type!\n", strHash.c_str(), strFixedMultiplier.c_str(), nTriedMultiplier);
    }
}

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(CBlock& block, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTests, unsigned int& nPrimesHit, mpz_class& mpzHash, CBlockIndex* pindexPrev, unsigned int vChainsFound[nMaxChainLength], CSieveOfEratosthenes& sieve, CPrimalityTestParams& testParams)
{
    nTests = 0;
    nPrimesHit = 0;

    // References to test parameters
    unsigned int& nBits = testParams.nBits;
    unsigned int& nChainLength = testParams.nChainLength;
    unsigned int& nCandidateType = testParams.nCandidateType;
    mpz_class& mpzHashFixedMult = testParams.mpzHashFixedMult;
    mpz_class& mpzChainOrigin = testParams.mpzChainOrigin;
    nBits = block.nBits;

    if (fNewBlock)
    {
        // Must rebuild the sieve
        sieve.Deplete();
    }
    fNewBlock = false;

    int64 nStart = 0; // microsecond timer
    if (!sieve.IsReady() || sieve.IsDepleted())
    {
        // Build sieve
        if (fDebug && GetBoolArg("-printmining"))
            nStart = GetTimeMicros();
        sieve.Reset(nSieveSize, nSieveFilterPrimes, nSieveExtensions, nL1CacheSize, nBits, mpzHash, mpzFixedMultiplier, pindexPrev);
        sieve.Weave();
        if (fDebug && GetBoolArg("-printmining"))
            printf("MineProbablePrimeChain() : new sieve (%u/%u@%u%%) ready in %uus\n", sieve.GetCandidateCount(), nSieveSize, sieve.GetProgressPercentage(), (unsigned int) (GetTimeMicros() - nStart));
        return false; // sieve generation takes time so return now
    }

    if (fDebug && GetBoolArg("-printmining2"))
        nStart = GetTimeMicros();

    // Number of candidates to be tested during a single call to this function
    const unsigned int nTestsAtOnce = 500;
    mpzHashFixedMult = mpzHash * mpzFixedMultiplier;

    // Process a part of the candidates
    while (nTests < nTestsAtOnce && pindexPrev == pindexBest)
    {
        unsigned int nTriedMultiplier = 0;
        if (!sieve.GetNextCandidateMultiplier(nTriedMultiplier, nCandidateType))
        {
            // power tests completed for the sieve
            if (fDebug && GetBoolArg("-printmining2"))
                printf("MineProbablePrimeChain() : %u tests (%u primes) in %uus\n", nTests, nPrimesHit, (unsigned int) (GetTimeMicros() - nStart));
            fNewBlock = true; // notify caller to change nonce
            return false;
        }
        nTests++;
        mpzChainOrigin = mpzHashFixedMult * nTriedMultiplier;
        bool fChainFound = ProbablePrimeChainTestFast(mpzChainOrigin, testParams);
        unsigned int nChainPrimeLength = TargetGetLength(nChainLength);

        if (fDebug && GetBoolArg("-debugsieve"))
            SieveDebugChecks(nBits, nTriedMultiplier, nCandidateType, mpzHash, mpzFixedMultiplier, mpzChainOrigin);

        // Collect mining statistics
        if(nChainPrimeLength >= 1)
        {
            nPrimesHit++;
            vChainsFound[nChainPrimeLength - 1]++;
        }

        // Check if a chain was found
        if (fChainFound)
        {
            mpz_class mpzPrimeChainMultiplier = mpzFixedMultiplier * nTriedMultiplier;
            CBigNum bnPrimeChainMultiplier;
            bnPrimeChainMultiplier.SetHex(mpzPrimeChainMultiplier.get_str(16));
            block.bnPrimeChainMultiplier = bnPrimeChainMultiplier;
            printf("nTriedMultiplier = %u\n", nTriedMultiplier); // Debugging
            printf("Probable prime chain found for block=%s!!\n  Target: %s\n  Chain: %s\n", block.GetHash().GetHex().c_str(),
                TargetToString(block.nBits).c_str(), GetPrimeChainName(nCandidateType, nChainLength).c_str());
            return true;
        }
    }
    
    if (fDebug && GetBoolArg("-printmining2"))
        printf("MineProbablePrimeChain() : %u tests (%u primes) in %uus\n", nTests, nPrimesHit, (unsigned int) (GetTimeMicros() - nStart));
    
    return false; // stop as new block arrived
}

// Get progress percentage of the sieve
unsigned int CSieveOfEratosthenes::GetProgressPercentage()
{
    return std::min(100u, (((nPrimeSeq >= vPrimes.size())? nSieveSize : vPrimes[nPrimeSeq]) * 100 / nSieveSize));
}

static unsigned int int_invert(unsigned int a, unsigned int nPrime)
{
    // Extended Euclidean algorithm to calculate the inverse of a in finite field defined by nPrime
    int rem0 = nPrime, rem1 = a % nPrime, rem2;
    int aux0 = 0, aux1 = 1, aux2;
    int quotient, inverse;

    while (1)
    {
        if (rem1 <= 1)
        {
            inverse = aux1;
            break;
        }

        rem2 = rem0 % rem1;
        quotient = rem0 / rem1;
        aux2 = -quotient * aux1 + aux0;

        if (rem2 <= 1)
        {
            inverse = aux2;
            break;
        }

        rem0 = rem1 % rem2;
        quotient = rem1 / rem2;
        aux0 = -quotient * aux2 + aux1;

        if (rem0 <= 1)
        {
            inverse = aux0;
            break;
        }

        rem1 = rem2 % rem0;
        quotient = rem2 / rem0;
        aux1 = -quotient * aux0 + aux2;
    }

    return (inverse + nPrime) % nPrime;
}

void CSieveOfEratosthenes::ProcessMultiplier(sieve_word_t *vfComposites, const unsigned int nMinMultiplier, const unsigned int nMaxMultiplier, const std::vector<unsigned int>& vPrimes, unsigned int *vMultipliers, unsigned int nLayerSeq)
{
    // Wipe the part of the array first
    if (nMinMultiplier < nMaxMultiplier)
    {
        memset(vfComposites + GetWordNum(nMinMultiplier), 0, (nMaxMultiplier - nMinMultiplier + nWordBits - 1) / nWordBits * sizeof(sieve_word_t));

        for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
        {
            const unsigned int nPrime = vPrimes[nPrimeSeq];
            unsigned int nVariableMultiplier = vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq];
            if (nVariableMultiplier < nMinMultiplier)
                nVariableMultiplier += (nMinMultiplier - nVariableMultiplier + nPrime - 1) / nPrime * nPrime;
#ifdef USE_ROTATE
            const unsigned int nRotateBits = nPrime % nWordBits;
            sieve_word_t lBitMask = GetBitMask(nVariableMultiplier);
            for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
            {
                vfComposites[GetWordNum(nVariableMultiplier)] |= lBitMask;
                lBitMask = rotate_left(lBitMask, nRotateBits);
            }
            vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq] = nVariableMultiplier;
#else
            for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
            {
                vfComposites[GetWordNum(nVariableMultiplier)] |= GetBitMask(nVariableMultiplier);
            }
            vMultipliers[nPrimeSeq * nSieveLayers + nLayerSeq] = nVariableMultiplier;
#endif
        }
    }
}

// Weave sieve for the next prime in table
// Return values:
//   True  - weaved another prime; nComposite - number of composites removed
//   False - sieve already completed
bool CSieveOfEratosthenes::Weave()
{
    // Check whether fixed multiplier fits in an unsigned long
    bool fUseLongForFixedMultiplier = mpzFixedMultiplier < ULONG_MAX;
    unsigned long nFixedMultiplier = mpzFixedMultiplier.get_ui();

    unsigned int nCombinedEndSeq = 1;
    unsigned int nFixedFactorCombinedMod = 0;

    for (unsigned int nPrimeSeqLocal = 1; nPrimeSeqLocal < nPrimes; nPrimeSeqLocal++)
    {
        if (pindexPrev != pindexBest)
            break;  // new block
        unsigned int nPrime = vPrimes[nPrimeSeqLocal];
        if (nPrimeSeqLocal >= nCombinedEndSeq)
        {
            // Combine multiple primes to produce a big divisor
            unsigned int nPrimeCombined = 1;
            while (nPrimeCombined < UINT_MAX / vPrimes[nCombinedEndSeq])
            {
                nPrimeCombined *= vPrimes[nCombinedEndSeq];
                nCombinedEndSeq++;
            }

            if (fUseLongForFixedMultiplier)
            {
                nFixedFactorCombinedMod = mpz_tdiv_ui(mpzHash.get_mpz_t(), nPrimeCombined);
                nFixedFactorCombinedMod = (uint64)nFixedFactorCombinedMod * (nFixedMultiplier % nPrimeCombined) % nPrimeCombined;
            }
            else
                nFixedFactorCombinedMod = mpz_tdiv_ui(mpzHashFixedMult.get_mpz_t(), nPrimeCombined);
        }

        unsigned int nFixedFactorMod = nFixedFactorCombinedMod % nPrime;
        if (nFixedFactorMod == 0)
        {
            // Nothing in the sieve is divisible by this prime
            continue;
        }
        // Find the modulo inverse of fixed factor
        unsigned int nFixedInverse = int_invert(nFixedFactorMod, nPrime);
        if (!nFixedInverse)
            return error("CSieveOfEratosthenes::Weave(): int_invert of fixed factor failed for prime #%u=%u", nPrimeSeqLocal, vPrimes[nPrimeSeqLocal]);
        unsigned int nTwoInverse = (nPrime + 1) / 2;

        // Check whether 32-bit arithmetic can be used for nFixedInverse
        const bool fUse32BArithmetic = (UINT_MAX / nTwoInverse) >= nPrime;

        if (fUse32BArithmetic)
        {
            // Weave the sieve for the prime
            for (unsigned int nChainSeq = 0; nChainSeq < nSieveLayers; nChainSeq++)
            {
                // Find the first number that's divisible by this prime
                vCunningham1Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nFixedInverse;
                vCunningham2Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nPrime - nFixedInverse;

                // For next number in chain
                nFixedInverse = nFixedInverse * nTwoInverse % nPrime;
            }
        }
        else
        {
            // Weave the sieve for the prime
            for (unsigned int nChainSeq = 0; nChainSeq < nSieveLayers; nChainSeq++)
            {
                // Find the first number that's divisible by this prime
                vCunningham1Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nFixedInverse;
                vCunningham2Multipliers[nPrimeSeqLocal * nSieveLayers + nChainSeq] = nPrime - nFixedInverse;

                // For next number in chain
                nFixedInverse = (uint64)nFixedInverse * nTwoInverse % nPrime;
            }
        }
    }

    // Process the array in chunks that fit the L1 cache
    const unsigned int nArrayRounds = (nSieveSize + nL1CacheElements - 1) / nL1CacheElements;

    // Calculate the number of CC1 and CC2 layers needed for BiTwin candidates
    const unsigned int nBiTwinCC1Layers = (nChainLength + 1) / 2;
    const unsigned int nBiTwinCC2Layers = nChainLength / 2;

    // Only 50% of the array is used in extensions
    const unsigned int nExtensionsMinMultiplier = nSieveSize / 2;
    const unsigned int nExtensionsMinWord = nExtensionsMinMultiplier / nWordBits;

    // Loop over each array one at a time for optimal L1 cache performance
    for (unsigned int j = 0; j < nArrayRounds; j++)
    {
        const unsigned int nMinMultiplier = nL1CacheElements * j;
        const unsigned int nMaxMultiplier = std::min(nL1CacheElements * (j + 1), nSieveSize);
        const unsigned int nExtMinMultiplier = std::max(nMinMultiplier, nExtensionsMinMultiplier);
        const unsigned int nMinWord = nMinMultiplier / nWordBits;
        const unsigned int nMaxWord = (nMaxMultiplier + nWordBits - 1) / nWordBits;
        const unsigned int nExtMinWord = std::max(nMinWord, nExtensionsMinWord);
        if (pindexPrev != pindexBest)
            break;  // new block

        // Loop over the layers
        for (unsigned int nLayerSeq = 0; nLayerSeq < nSieveLayers; nLayerSeq++) {
            if (pindexPrev != pindexBest)
                break;  // new block
            if (nLayerSeq < nChainLength)
            {
                ProcessMultiplier(vfCompositeLayerCC1, nMinMultiplier, nMaxMultiplier, vPrimes, vCunningham1Multipliers, nLayerSeq);
                ProcessMultiplier(vfCompositeLayerCC2, nMinMultiplier, nMaxMultiplier, vPrimes, vCunningham2Multipliers, nLayerSeq);
            }
            else
            {
                // Optimize: First halves of the arrays are not needed in the extensions
                ProcessMultiplier(vfCompositeLayerCC1, nExtMinMultiplier, nMaxMultiplier, vPrimes, vCunningham1Multipliers, nLayerSeq);
                ProcessMultiplier(vfCompositeLayerCC2, nExtMinMultiplier, nMaxMultiplier, vPrimes, vCunningham2Multipliers, nLayerSeq);
            }

            // Apply the layer to the primary sieve arrays
            if (nLayerSeq < nChainLength)
            {
                if (nLayerSeq < nBiTwinCC2Layers)
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                        vfCompositeBiTwin[nWord] |= vfCompositeLayerCC1[nWord] | vfCompositeLayerCC2[nWord];
                    }
                }
                else if (nLayerSeq < nBiTwinCC1Layers)
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                        vfCompositeBiTwin[nWord] |= vfCompositeLayerCC1[nWord];
                    }
                }
                else
                {
                    for (unsigned int nWord = nMinWord; nWord < nMaxWord; nWord++)
                    {
                        vfCompositeCunningham1[nWord] |= vfCompositeLayerCC1[nWord];
                        vfCompositeCunningham2[nWord] |= vfCompositeLayerCC2[nWord];
                    }
                }
            }

            // Apply the layer to extensions
            for (unsigned int nExtensionSeq = 0; nExtensionSeq < nSieveExtensions; nExtensionSeq++)
            {
                const unsigned int nLayerOffset = nExtensionSeq + 1;
                if (nLayerSeq >= nLayerOffset && nLayerSeq < nChainLength + nLayerOffset)
                {
                    const unsigned int nLayerExtendedSeq = nLayerSeq - nLayerOffset;
                    sieve_word_t *vfExtCC1 = vfExtendedCompositeCunningham1 + nExtensionSeq * nCandidatesWords;
                    sieve_word_t *vfExtCC2 = vfExtendedCompositeCunningham2 + nExtensionSeq * nCandidatesWords;
                    sieve_word_t *vfExtTWN = vfExtendedCompositeBiTwin + nExtensionSeq * nCandidatesWords;
                    if (nLayerExtendedSeq < nBiTwinCC2Layers)
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                            vfExtTWN[nWord] |= vfCompositeLayerCC1[nWord] | vfCompositeLayerCC2[nWord];
                        }
                    }
                    else if (nLayerExtendedSeq < nBiTwinCC1Layers)
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                            vfExtTWN[nWord] |= vfCompositeLayerCC1[nWord];
                        }
                    }
                    else
                    {
                        for (unsigned int nWord = nExtMinWord; nWord < nMaxWord; nWord++)
                        {
                            vfExtCC1[nWord] |= vfCompositeLayerCC1[nWord];
                            vfExtCC2[nWord] |= vfCompositeLayerCC2[nWord];
                        }
                    }
                }
            }
        }

        // Combine the bitsets
        // vfCandidates = ~(vfCompositeCunningham1 & vfCompositeCunningham2 & vfCompositeBiTwin)
        for (unsigned int i = nMinWord; i < nMaxWord; i++)
            vfCandidates[i] = ~(vfCompositeCunningham1[i] & vfCompositeCunningham2[i] & vfCompositeBiTwin[i]);

        // Combine the extended bitsets
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = nExtMinWord; i < nMaxWord; i++)
                vfExtendedCandidates[j * nCandidatesWords + i] = ~(
                    vfExtendedCompositeCunningham1[j * nCandidatesWords + i] &
                    vfExtendedCompositeCunningham2[j * nCandidatesWords + i] &
                    vfExtendedCompositeBiTwin[j * nCandidatesWords + i]);
    }

    // The sieve has been partially weaved
    this->nPrimeSeq = nPrimes - 1;

    return false;
}

static const double dLogTwo = log(2.0);
static const double dLogOneAndHalf = log(1.5);

// Estimate the probability of primality for a number in a candidate chain
double EstimateCandidatePrimeProbability(unsigned int nPrimorialMultiplier, unsigned int nChainPrimeNum, unsigned int nMiningProtocol)
{
    // h * q# / r# * s is prime with probability 1/log(h * q# / r# * s),
    //   (prime number theorem)
    //   here s ~ max sieve size / 2,
    //   h ~ 2^255 * 1.5,
    //   r = 7 (primorial multiplier embedded in the hash)
    // Euler product to p ~ 1.781072 * log(p)   (Mertens theorem)
    // If sieve is weaved up to p, a number in a candidate chain is a prime
    // with probability
    //     (1/log(h * q# / r# * s)) / (1/(1.781072 * log(p)))
    //   = 1.781072 * log(p) / (255 * log(2) + log(1.5) + log(q# / r#) + log(s))
    //
    // This model assumes that the numbers on a chain being primes are
    // statistically independent after running the sieve, which might not be
    // true, but nontheless it's a reasonable model of the chances of finding
    // prime chains.
    const unsigned int nSieveWeaveOptimalPrime = vPrimes[nSieveFilterPrimes - 1];
    const unsigned int nAverageCandidateMultiplier = nSieveSize / 2;
    double dFixedMultiplier = 1.0;
    for (unsigned int i = 0; vPrimes[i] <= nPrimorialMultiplier; i++)
        dFixedMultiplier *= vPrimes[i];
    if (nMiningProtocol < 2)
    {
        for (unsigned int i = 0; vPrimes[i] <= nPrimorialHashFactor; i++)
            dFixedMultiplier /= vPrimes[i];
    }

    double dExtendedSieveWeightedSum = 0.5 * nSieveSize;
    double dExtendedSieveCandidates = nSieveSize;
    for (unsigned int i = 0; i < nSieveExtensions; i++)
    {
        dExtendedSieveWeightedSum += 0.75 * (nSieveSize * (2 << i));
        dExtendedSieveCandidates += nSieveSize / 2;
    }
    const double dExtendedSieveAverageMultiplier = dExtendedSieveWeightedSum / dExtendedSieveCandidates;

    return (1.781072 * log((double)std::max(1u, nSieveWeaveOptimalPrime)) / (255.0 * dLogTwo + dLogOneAndHalf + log(dFixedMultiplier) + log(nAverageCandidateMultiplier) + dLogTwo * nChainPrimeNum + log(dExtendedSieveAverageMultiplier)));
}

// Esimate the prime probablity of numbers that haven't been sieved
double EstimateNormalPrimeProbability(unsigned int nPrimorialMultiplier, unsigned int nChainPrimeNum, unsigned int nMiningProtocol)
{
    const unsigned int nAverageCandidateMultiplier = nSieveSize / 2;
    double dFixedMultiplier = 1.0;
    for (unsigned int i = 0; vPrimes[i] <= nPrimorialMultiplier; i++)
        dFixedMultiplier *= vPrimes[i];
    if (nMiningProtocol < 2)
    {
        for (unsigned int i = 0; vPrimes[i] <= nPrimorialHashFactor; i++)
            dFixedMultiplier /= vPrimes[i];
    }

    double dExtendedSieveWeightedSum = 0.5 * nSieveSize;
    double dExtendedSieveCandidates = nSieveSize;
    for (unsigned int i = 0; i < nSieveExtensions; i++)
    {
        dExtendedSieveWeightedSum += 0.75 * (nSieveSize * (2 << i));
        dExtendedSieveCandidates += nSieveSize / 2;
    }
    const double dExtendedSieveAverageMultiplier = dExtendedSieveWeightedSum / dExtendedSieveCandidates;

    // The primorial is implicitly filtering out the first few prime factors
    double dPrimorialBoost = 1.0;
    for (unsigned int i = 0; vPrimes[i] <= nPrimorialMultiplier; i++)
        dPrimorialBoost *= (double)vPrimes[i] / (vPrimes[i] - 1);

    return (dPrimorialBoost / (255.0 * dLogTwo + dLogOneAndHalf + log(dFixedMultiplier) + log(nAverageCandidateMultiplier) + dLogTwo * nChainPrimeNum + log(dExtendedSieveAverageMultiplier)));
}
