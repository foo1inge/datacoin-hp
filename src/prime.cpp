// Copyright (c) 2013 Primecoin developers
// Distributed under conditional MIT/X11 software license,
// see the accompanying file COPYING

#include "prime.h"

// Prime Table
std::vector<unsigned int> vPrimes;
std::vector<unsigned int> vTwoInverses;
static unsigned int nSieveSize = nDefaultSieveSize;

void GeneratePrimeTable()
{
    nSieveSize = (int)GetArg("-sievesize", nDefaultSieveSize);
    nSieveSize = std::max(std::min(nSieveSize, nMaxSieveSize), nMinSieveSize);
    const unsigned nPrimeTableLimit = nSieveSize;
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
    printf("GeneratePrimeTable() : prime table [1, %d] generated with %lu primes", nPrimeTableLimit, vPrimes.size());
    //BOOST_FOREACH(unsigned int nPrime, vPrimes)
    //    printf(" %u", nPrime);
    printf("\n");
    
    mpz_t p;
    mpz_t mpzTwoInverse;
    mpz_init(p);
    mpz_init(mpzTwoInverse);
    const unsigned int nPrimes = vPrimes.size();
    vTwoInverses = std::vector<unsigned int> (nPrimes, 0);
    for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
    {
        mpz_set_ui(p, vPrimes[nPrimeSeq]);
        if (!mpz_invert(mpzTwoInverse, mpzTwo.get_mpz_t(), p))
            printf("GeneratePrimeTable(): mpz_invert of 2 failed for prime #%u=%u", nPrimeSeq, vPrimes[nPrimeSeq]);
        vTwoInverses[nPrimeSeq] = mpz_get_ui(mpzTwoInverse);
    }
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
    mpzPrimorial = 1;
    BOOST_FOREACH(unsigned int nPrime, vPrimes)
    {
        if (nPrime > p) break;
        mpzPrimorial *= nPrime;
    }
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

// Check Fermat probable primality test (2-PRP): 2 ** (n-1) = 1 (mod n)
// true: n is probable prime
// false: n is composite; set fractional length in the nLength output
static bool FermatProbablePrimalityTest(const mpz_class& n, unsigned int& nLength)
{
    // Faster GMP version
    
    mpz_t mpzN;
    mpz_t mpzE;
    mpz_t mpzR;
    
    mpz_init_set(mpzN, n.get_mpz_t());
    mpz_init(mpzE);
    mpz_sub_ui(mpzE, mpzN, 1);
    mpz_init(mpzR);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, mpzN);
    if (mpz_cmp_ui(mpzR, 1) == 0)
    {
        mpz_clear(mpzN);
        mpz_clear(mpzE);
        mpz_clear(mpzR);
        return true;
    }
    // Failed Fermat test, calculate fractional length
    mpz_sub(mpzE, mpzN, mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, mpzN);
    unsigned int nFractionalLength = mpz_get_ui(mpzE);
    mpz_clear(mpzN);
    mpz_clear(mpzE);
    mpz_clear(mpzR);
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
static bool EulerLagrangeLifchitzPrimalityTest(const mpz_class& n, bool fSophieGermain, unsigned int& nLength)
{
    // Faster GMP version
    mpz_t mpzN;
    mpz_t mpzE;
    mpz_t mpzR;
    
    mpz_init_set(mpzN, n.get_mpz_t());
    mpz_init(mpzE);
    mpz_sub_ui(mpzE, mpzN, 1);
    mpz_tdiv_q_2exp(mpzE, mpzE, 1);
    mpz_init(mpzR);
    mpz_powm(mpzR, mpzTwo.get_mpz_t(), mpzE, mpzN);
    unsigned int nMod8 = mpz_tdiv_ui(mpzN, 8);
    bool fPassedTest = false;
    if (fSophieGermain && (nMod8 == 7)) // Euler & Lagrange
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    else if (fSophieGermain && (nMod8 == 3)) // Lifchitz
    {
        mpz_t mpzRplusOne;
        mpz_init(mpzRplusOne);
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, mpzN);
        mpz_clear(mpzRplusOne);
    }
    else if ((!fSophieGermain) && (nMod8 == 5)) // Lifchitz
    {
        mpz_t mpzRplusOne;
        mpz_init(mpzRplusOne);
        mpz_add_ui(mpzRplusOne, mpzR, 1);
        fPassedTest = !mpz_cmp(mpzRplusOne, mpzN);
        mpz_clear(mpzRplusOne);
    }
    else if ((!fSophieGermain) && (nMod8 == 1)) // LifChitz
    {
        fPassedTest = !mpz_cmp_ui(mpzR, 1);
    }
    else
    {
        mpz_clear(mpzN);
        mpz_clear(mpzE);
        mpz_clear(mpzR);
        return error("EulerLagrangeLifchitzPrimalityTest() : invalid n %% 8 = %d, %s", nMod8, (fSophieGermain? "first kind" : "second kind"));
    }
    
    if (fPassedTest)
    {
        mpz_clear(mpzN);
        mpz_clear(mpzE);
        mpz_clear(mpzR);
        return true;
    }
    
    // Failed test, calculate fractional length
    mpz_mul(mpzE, mpzR, mpzR);
    mpz_tdiv_r(mpzR, mpzE, mpzN); // derive Fermat test remainder

    mpz_sub(mpzE, mpzN, mpzR);
    mpz_mul_2exp(mpzR, mpzE, nFractionalBits);
    mpz_tdiv_q(mpzE, mpzR, mpzN);
    unsigned int nFractionalLength = mpz_get_ui(mpzE);
    mpz_clear(mpzN);
    mpz_clear(mpzE);
    mpz_clear(mpzR);
    
    if (nFractionalLength >= (1 << nFractionalBits))
        return error("EulerLagrangeLifchitzPrimalityTest() : fractional assert");
    nLength = (nLength & TARGET_LENGTH_MASK) | nFractionalLength;
    return false;
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

void TargetIncrementLength(unsigned int& nBits)
{
    nBits += (1 << nFractionalBits);
}

void TargetDecrementLength(unsigned int& nBits)
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

// Test Probable Cunningham Chain for: n
// fSophieGermain:
//   true - Test for Cunningham Chain of first kind (n, 2n+1, 4n+3, ...)
//   false - Test for Cunningham Chain of second kind (n, 2n-1, 4n-3, ...)
// Return value:
//   true - Probable Cunningham Chain found (length at least 2)
//   false - Not Cunningham Chain
static bool ProbableCunninghamChainTest(const mpz_class& n, bool fSophieGermain, bool fFermatTest, unsigned int& nProbableChainLength)
{
    nProbableChainLength = 0;
    mpz_class N = n;

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
bool ProbablePrimeChainTest(const mpz_class& mpzPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int& nChainLengthCunningham1, unsigned int& nChainLengthCunningham2, unsigned int& nChainLengthBiTwin)
{
    nChainLengthCunningham1 = 0;
    nChainLengthCunningham2 = 0;
    nChainLengthBiTwin = 0;

    // Test for Cunningham Chain of first kind
    ProbableCunninghamChainTest(mpzPrimeChainOrigin-1, true, fFermatTest, nChainLengthCunningham1);
    // Test for Cunningham Chain of second kind
    ProbableCunninghamChainTest(mpzPrimeChainOrigin+1, false, fFermatTest, nChainLengthCunningham2);
    // Figure out BiTwin Chain length
    // BiTwin Chain allows a single prime at the end for odd length chain
    nChainLengthBiTwin =
        (TargetGetLength(nChainLengthCunningham1) > TargetGetLength(nChainLengthCunningham2))?
            (nChainLengthCunningham2 + TargetFromInt(TargetGetLength(nChainLengthCunningham2)+1)) :
            (nChainLengthCunningham1 + TargetFromInt(TargetGetLength(nChainLengthCunningham1)));

    return (nChainLengthCunningham1 >= nBits || nChainLengthCunningham2 >= nBits || nChainLengthBiTwin >= nBits);
}

// Sieve for mining
boost::thread_specific_ptr<CSieveOfEratosthenes> psieve;

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(CBlock& block, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTriedMultiplier, unsigned int& nProbableChainLength, unsigned int& nTests, unsigned int& nPrimesHit, unsigned int& nChainsHit, mpz_class& mpzHash)
{
    nProbableChainLength = 0;
    nTests = 0;
    nPrimesHit = 0;
    nChainsHit = 0;

    if (fNewBlock && psieve.get() != NULL)
    {
        // Must rebuild the sieve
        psieve.reset();
    }
    fNewBlock = false;

    int64 nStart, nCurrent; // microsecond timer
    CBlockIndex* pindexPrev = pindexBest;
    if (psieve.get() == NULL)
    {
        // Build sieve
        nStart = GetTimeMicros();
        CSieveOfEratosthenes *lpsieve = new CSieveOfEratosthenes(nSieveSize, block.nBits, mpzHash, mpzFixedMultiplier, pindexPrev);
        int64 nSieveRoundLimit = (int)GetArg("-gensieveroundlimitms", 1000);
        while (lpsieve->Weave() && pindexPrev == pindexBest && (GetTimeMicros() - nStart < 1000 * nSieveRoundLimit));
        if (fDebug && GetBoolArg("-printmining"))
            printf("MineProbablePrimeChain() : new sieve (%u/%u@%u%%) ready in %uus\n", lpsieve->GetCandidateCount(), nSieveSize, lpsieve->GetProgressPercentage(), (unsigned int) (GetTimeMicros() - nStart));
        psieve.reset(lpsieve);
    }

    mpz_class mpzChainOrigin;

    nStart = GetTimeMicros();
    nCurrent = nStart;

    while (nCurrent - nStart < 10000 && nCurrent >= nStart && pindexPrev == pindexBest)
    {
        nTests++;
        if (!psieve->GetNextCandidateMultiplier(nTriedMultiplier))
        {
            // power tests completed for the sieve
            psieve.reset();
            fNewBlock = true; // notify caller to change nonce
            return false;
        }
        mpzChainOrigin = mpzHash * mpzFixedMultiplier * nTriedMultiplier;
        unsigned int nChainLengthCunningham1 = 0;
        unsigned int nChainLengthCunningham2 = 0;
        unsigned int nChainLengthBiTwin = 0;
        if (ProbablePrimeChainTest(mpzChainOrigin, block.nBits, false, nChainLengthCunningham1, nChainLengthCunningham2, nChainLengthBiTwin))
        {
            mpz_class mpzPrimeChainMultiplier = mpzFixedMultiplier * nTriedMultiplier;
            CBigNum bnPrimeChainMultiplier;
            bnPrimeChainMultiplier.SetHex(mpzPrimeChainMultiplier.get_str(16));
            block.bnPrimeChainMultiplier = bnPrimeChainMultiplier;
            printf("Probable prime chain found for block=%s!!\n  Target: %s\n  Length: (%s %s %s)\n", block.GetHash().GetHex().c_str(),
            TargetToString(block.nBits).c_str(), TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str());
            nProbableChainLength = std::max(std::max(nChainLengthCunningham1, nChainLengthCunningham2), nChainLengthBiTwin);
            return true;
        }
        nProbableChainLength = std::max(std::max(nChainLengthCunningham1, nChainLengthCunningham2), nChainLengthBiTwin);
        if(TargetGetLength(nProbableChainLength) >= 1)
            nPrimesHit++;
        if(TargetGetLength(nProbableChainLength) >= nStatsChainLength)
            nChainsHit++;

        nCurrent = GetTimeMicros();
    }
    return false; // stop as timed out
}

// Check prime proof-of-work
bool CheckPrimeProofOfWork(uint256& hashBlockHeader, unsigned int nBits, const mpz_class& mpzPrimeChainMultiplier, unsigned int& nChainType, unsigned int& nChainLength)
{
    // Check target
    if (TargetGetLength(nBits) < nTargetMinLength || TargetGetLength(nBits) > 99)
        return error("CheckPrimeProofOfWork() : invalid chain length target %s", TargetToString(nBits).c_str());

    // Check header hash limit
    if (hashBlockHeader < hashBlockHeaderLimit)
        return error("CheckPrimeProofOfWork() : block header hash under limit");
    mpz_class mpzHashBlockHeader;
    mpz_set_uint256(mpzHashBlockHeader.get_mpz_t(), hashBlockHeader);
    // Check target for prime proof-of-work
    mpz_class mpzPrimeChainOrigin = mpzHashBlockHeader * mpzPrimeChainMultiplier;
    if (mpzPrimeChainOrigin < mpzPrimeMin)
        return error("CheckPrimeProofOfWork() : prime too small");
    // First prime in chain must not exceed cap
    if (mpzPrimeChainOrigin > mpzPrimeMax)
        return error("CheckPrimeProofOfWork() : prime too big");

    // Check prime chain
    unsigned int nChainLengthCunningham1 = 0;
    unsigned int nChainLengthCunningham2 = 0;
    unsigned int nChainLengthBiTwin = 0;
    if (!ProbablePrimeChainTest(mpzPrimeChainOrigin, nBits, false, nChainLengthCunningham1, nChainLengthCunningham2, nChainLengthBiTwin))
        return error("CheckPrimeProofOfWork() : failed prime chain test target=%s length=(%s %s %s)", TargetToString(nBits).c_str(),
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str());
    if (nChainLengthCunningham1 < nBits && nChainLengthCunningham2 < nBits && nChainLengthBiTwin < nBits)
        return error("CheckPrimeProofOfWork() : prime chain length assert target=%s length=(%s %s %s)", TargetToString(nBits).c_str(),
            TargetToString(nChainLengthCunningham1).c_str(), TargetToString(nChainLengthCunningham2).c_str(), TargetToString(nChainLengthBiTwin).c_str());

    // Double check prime chain with Fermat tests only
    unsigned int nChainLengthCunningham1FermatTest = 0;
    unsigned int nChainLengthCunningham2FermatTest = 0;
    unsigned int nChainLengthBiTwinFermatTest = 0;
    if (!ProbablePrimeChainTest(mpzPrimeChainOrigin, nBits, true, nChainLengthCunningham1FermatTest, nChainLengthCunningham2FermatTest, nChainLengthBiTwinFermatTest))
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

    // Check that the certificate (mpzPrimeChainMultiplier) is normalized
    if (mpzPrimeChainMultiplier % 2 == 0 && mpzPrimeChainOrigin % 4 == 0)
    {
        unsigned int nChainLengthCunningham1Extended = 0;
        unsigned int nChainLengthCunningham2Extended = 0;
        unsigned int nChainLengthBiTwinExtended = 0;
        if (ProbablePrimeChainTest(mpzPrimeChainOrigin / 2, nBits, false, nChainLengthCunningham1Extended, nChainLengthCunningham2Extended, nChainLengthBiTwinExtended))
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

// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength)
{
    return strprintf("%s%s", (nChainType==PRIME_CHAIN_CUNNINGHAM1)? "1CC" : ((nChainType==PRIME_CHAIN_CUNNINGHAM2)? "2CC" : "TWN"), TargetToString(nChainLength).c_str());
}


// Get progress percentage of the sieve
unsigned int CSieveOfEratosthenes::GetProgressPercentage()
{
    return std::min(100u, (((nPrimeSeq >= vPrimes.size())? nSieveSize : vPrimes[nPrimeSeq]) * 100 / nSieveSize));
}

// Weave sieve for the next prime in table
// Return values:
//   True  - weaved another prime; nComposite - number of composites removed
//   False - sieve already completed
bool CSieveOfEratosthenes::Weave()
{
    // Faster GMP version
    const unsigned int nChainLength = TargetGetLength(nBits);
    const unsigned int nHalfChainLength = (nChainLength + 1) / 2;
    const unsigned int nTotalPrimes = vPrimes.size();

    // Keep all variables local for max performance
    CBlockIndex* pindexPrev = this->pindexPrev;
    unsigned int nSieveSize = this->nSieveSize;

    // Process only 10% of the primes
    // Most composites are still found
    const unsigned int nPrimes = (uint64)nTotalPrimes * 10 / 100;

    mpz_t mpzFixedFactor; // fixed factor to derive the chain
    mpz_t mpzFixedFactorMod;
    mpz_t p;
    mpz_t mpzFixedInverse;

    unsigned int vCunningham1AMultipliers[nPrimes][nHalfChainLength];
    unsigned int vCunningham1BMultipliers[nPrimes][nHalfChainLength];
    unsigned int vCunningham2AMultipliers[nPrimes][nHalfChainLength];
    unsigned int vCunningham2BMultipliers[nPrimes][nHalfChainLength];
    
    mpz_init_set(mpzFixedFactor, this->mpzFixedFactor.get_mpz_t());
    mpz_init(mpzFixedFactorMod);
    mpz_init(p);
    mpz_init(mpzFixedInverse);
    
    memset(vCunningham1AMultipliers, 0, sizeof(vCunningham1AMultipliers));
    memset(vCunningham1BMultipliers, 0, sizeof(vCunningham1BMultipliers));
    memset(vCunningham2AMultipliers, 0, sizeof(vCunningham2AMultipliers));
    memset(vCunningham2BMultipliers, 0, sizeof(vCunningham2BMultipliers));

    for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
    {
        if (pindexPrev != pindexBest)
            break;  // new block
        unsigned int nPrime = vPrimes[nPrimeSeq];
        unsigned long nFixedFactorMod = mpz_tdiv_r_ui(mpzFixedFactorMod, mpzFixedFactor, nPrime);
        if (nFixedFactorMod == 0)
        {
            // Nothing in the sieve is divisible by this prime
            continue;
        }
        mpz_set_ui(p, nPrime);
        // Find the modulo inverse of fixed factor
        if (!mpz_invert(mpzFixedInverse, mpzFixedFactorMod, p))
            return error("CSieveOfEratosthenes::Weave(): mpz_invert of fixed factor failed for prime #%u=%u", nPrimeSeq, vPrimes[nPrimeSeq]);
        unsigned long nFixedInverse = mpz_get_ui(mpzFixedInverse);
        unsigned long nTwoInverse = vTwoInverses[nPrimeSeq];

        // Weave the sieve for the prime
        for (unsigned int nBiTwinSeq = 0; nBiTwinSeq < 2 * nChainLength; nBiTwinSeq++)
        {
            // Find the first number that's divisible by this prime
            int nDelta = ((nBiTwinSeq % 2 == 0)? (-1) : 1);
            unsigned int nSolvedMultiplier = (uint64)nFixedInverse * (nPrime - nDelta) % nPrime;
            if (nBiTwinSeq % 2 == 1)
                nFixedInverse = (uint64)nFixedInverse * nTwoInverse % nPrime;

            if (nBiTwinSeq < nChainLength)
            {
                if (((nBiTwinSeq & 1u) == 0))
                    vCunningham1AMultipliers[nPrimeSeq][nBiTwinSeq / 2] = nSolvedMultiplier;
                else
                    vCunningham2AMultipliers[nPrimeSeq][nBiTwinSeq / 2] = nSolvedMultiplier;
            } else {
                if (((nBiTwinSeq & 1u) == 0))
                    vCunningham1BMultipliers[nPrimeSeq][(nBiTwinSeq - nChainLength) / 2] = nSolvedMultiplier;
                else
                    vCunningham2BMultipliers[nPrimeSeq][(nBiTwinSeq - nChainLength) / 2] = nSolvedMultiplier;
            }
        }
    }
    
    // Number of elements that are likely to fit in L1 cache
    const unsigned int nL1CacheElements = 200000;
    const unsigned int nArrayRounds = (nSieveSize + nL1CacheElements - 1) / nL1CacheElements;

    // Loop over each array one at a time for optimal L1 cache performance
    for (unsigned int j = 0; j < nArrayRounds; j++)
    {
        const unsigned int nMinMultiplier = nL1CacheElements * j;
        const unsigned int nMaxMultiplier = std::min(nL1CacheElements * (j + 1), nSieveSize);
        
        if (pindexPrev != pindexBest)
            break;  // new block
        
        for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
        {
            unsigned int nPrime = vPrimes[nPrimeSeq];
            for (unsigned int i = 0; i < nHalfChainLength; i++)
            {
                unsigned int nVariableMultiplier = vCunningham1AMultipliers[nPrimeSeq][i];
                if (nVariableMultiplier == 0) break;
                for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
                    vfCompositeCunningham1A[nVariableMultiplier] = true;
                vCunningham1AMultipliers[nPrimeSeq][i] = nVariableMultiplier;
            }
        }
        
        for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
        {
            unsigned int nPrime = vPrimes[nPrimeSeq];
            for (unsigned int i = 0; i < nHalfChainLength; i++)
            {
                unsigned int nVariableMultiplier = vCunningham1BMultipliers[nPrimeSeq][i];
                if (nVariableMultiplier == 0) break;
                for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
                    vfCompositeCunningham1B[nVariableMultiplier] = true;
                vCunningham1BMultipliers[nPrimeSeq][i] = nVariableMultiplier;
            }
        }
        
        for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
        {
            unsigned int nPrime = vPrimes[nPrimeSeq];
            for (unsigned int i = 0; i < nHalfChainLength; i++)
            {
                unsigned int nVariableMultiplier = vCunningham2AMultipliers[nPrimeSeq][i];
                if (nVariableMultiplier == 0) break;
                for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
                    vfCompositeCunningham2A[nVariableMultiplier] = true;
                vCunningham2AMultipliers[nPrimeSeq][i] = nVariableMultiplier;
            }
        }
        
        for (unsigned int nPrimeSeq = 1; nPrimeSeq < nPrimes; nPrimeSeq++)
        {
            unsigned int nPrime = vPrimes[nPrimeSeq];
            for (unsigned int i = 0; i < nHalfChainLength; i++)
            {
                unsigned int nVariableMultiplier = vCunningham2BMultipliers[nPrimeSeq][i];
                if (nVariableMultiplier == 0) break;
                for (; nVariableMultiplier < nMaxMultiplier; nVariableMultiplier += nPrime)
                    vfCompositeCunningham2B[nVariableMultiplier] = true;
                vCunningham2BMultipliers[nPrimeSeq][i] = nVariableMultiplier;
            }
        }
        
        // Combine all the bitsets
        // vfCompositeCunningham1 = vfCompositeCunningham1A | vfCompositeCunningham1B
        // vfCompositeCunningham2 = vfCompositeCunningham2A | vfCompositeCunningham2B
        // vfCompositeBiTwin = vfCompositeCunningham1A | vfCompositeCunningham2A
        // vfCandidates = ~(vfCompositeCunningham1 & vfCompositeCunningham2 & vfCompositeBiTwin)
        {
            // Fast version
            const unsigned int nBytes = (nMaxMultiplier - nMinMultiplier) / 8;
            unsigned long *lCandidates = (unsigned long *)&vfCandidates + (nMinMultiplier / 8 / sizeof(unsigned long));
            unsigned long *lCompositeCunningham1A = (unsigned long *)&vfCompositeCunningham1A + (nMinMultiplier / 8 / sizeof(unsigned long));
            unsigned long *lCompositeCunningham1B = (unsigned long *)&vfCompositeCunningham1B + (nMinMultiplier / 8 / sizeof(unsigned long));
            unsigned long *lCompositeCunningham2A = (unsigned long *)&vfCompositeCunningham2A + (nMinMultiplier / 8 / sizeof(unsigned long));
            unsigned long *lCompositeCunningham2B = (unsigned long *)&vfCompositeCunningham2B + (nMinMultiplier / 8 / sizeof(unsigned long));
            const unsigned int nLongs = (nBytes + sizeof(unsigned long) + 1) / sizeof(unsigned long);
            for (unsigned int i = 0; i < nLongs; i++)
            {
                lCandidates[i] = ~((lCompositeCunningham1A[i] | lCompositeCunningham1B[i]) &
                                (lCompositeCunningham2A[i] | lCompositeCunningham2B[i]) &
                                (lCompositeCunningham1A[i] | lCompositeCunningham2A[i]));
            }
        }
    }
    
    this->nPrimeSeq = nPrimes - 1;
    
    mpz_clear(mpzFixedFactor);
    mpz_clear(mpzFixedFactorMod);
    mpz_clear(p);
    mpz_clear(mpzFixedInverse);
    
    return false;
}

