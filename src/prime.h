// Copyright (c) 2013 Primecoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef PRIMECOIN_PRIME_H
#define PRIMECOIN_PRIME_H

#include "main.h"
#include "base58.h"

#include <gmp.h>
#include <gmpxx.h>
#include <bitset>
#include <boost/timer/timer.hpp>

/**********************/
/* PRIMECOIN PROTOCOL */
/**********************/

extern std::vector<unsigned int> vPrimes;
static const unsigned int nMaxSieveExtensions = 20;
static const unsigned int nMinSieveExtensions = 0;
static const unsigned int nDefaultSieveExtensions = 9;
static const unsigned int nDefaultSieveExtensionsTestnet = 4;
extern unsigned int nSieveExtensions;
static const unsigned int nMaxSieveFilterPrimes = 78498u; // size of prime table
static const unsigned int nDefaultSieveFilterPrimes = 7849u;
static const unsigned int nMinSieveFilterPrimes = 1000u;
extern unsigned int nSieveFilterPrimes;
static const unsigned int nMaxSieveSize = 10000000u;
static const unsigned int nDefaultSieveSize = 917504u;
static const unsigned int nMinSieveSize = 100000u;
extern unsigned int nSieveSize;
static const unsigned int nMaxL1CacheSize = 128000u;
static const unsigned int nDefaultL1CacheSize = 28672u;
static const unsigned int nMinL1CacheSize = 12000u;
extern unsigned int nL1CacheSize;
static const uint256 hashBlockHeaderLimit = (uint256(1) << 255);
static const CBigNum bnOne = 1;
static const CBigNum bnPrimeMax = (bnOne << 2000) - 1;
static const CBigNum bnPrimeMin = (bnOne << 255);
static const mpz_class mpzOne = 1;
static const mpz_class mpzTwo = 2;
static const mpz_class mpzPrimeMax = (mpzOne << 2000) - 1;
static const mpz_class mpzPrimeMin = (mpzOne << 255);
static const unsigned int nPrimorialHashFactor = 7;
static const unsigned int nInitialPrimorialMultiplier = 47;
static const unsigned int nInitialPrimorialMultiplierTestnet = 17;

// Mining statistics
static const unsigned int nMaxChainLength = 24;
extern uint64 nTotalTests;
extern unsigned int nTotalBlocksFound;
extern std::vector<uint64> vTotalChainsFound;
extern boost::timer::cpu_timer minerTimer;
static const unsigned int nDefaultSieveTargetLength = -1;
extern int nSieveTargetLength;

// Primecoin HP: Optional automatic donations with every block found
static const std::string strDefaultDonationPercentage = "0.0";
static const double dMinDonationPercentage = 0.1;
static const double dMaxDonationPercentage = 99.9;
static const std::string strDefaultDonationAddress = "APRimehpExTQJ5fssWa5r84nYEnwGYffXK";
static const std::string strDefaultDonationAddressTestnet = "mpriMeHPGWdXYZNtM3ib4WqaViKnrZafkY";
extern CBitcoinAddress donationAddress;
extern double dDonationPercentage;

extern unsigned int nTargetInitialLength;
extern unsigned int nTargetMinLength;

// Generate small prime table
void GeneratePrimeTable();
// Reset the miner statistics
void ResetMinerStatistics();
// Initialize the miner
void InitPrimeMiner();
// Print miner statistics
void PrintMinerStatistics();
// Print compact statistics
void PrintCompactStatistics(volatile unsigned int vFoundChainCounter[nMaxChainLength]);
// Get next prime number of p
bool PrimeTableGetNextPrime(unsigned int& p);
// Get previous prime number of p
bool PrimeTableGetPreviousPrime(unsigned int& p);

// Compute primorial number p#
void Primorial(unsigned int p, mpz_class& mpzPrimorial);
// Compute Primorial number p#
// Fast 32-bit version assuming that p <= 23
unsigned int PrimorialFast(unsigned int p);
// Compute the first primorial number greater than or equal to bn
void PrimorialAt(mpz_class& bn, mpz_class& mpzPrimorial);

// Test probable prime chain for: bnPrimeChainOrigin
// fFermatTest
//   true - Use only Fermat tests
//   false - Use Fermat-Euler-Lagrange-Lifchitz tests
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
bool ProbablePrimeChainTest(const CBigNum& bnPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int& nChainLengthCunningham1, unsigned int& nChainLengthCunningham2, unsigned int& nChainLengthBiTwin);

static const unsigned int nFractionalBits = 24;
static const unsigned int TARGET_FRACTIONAL_MASK = (1u<<nFractionalBits) - 1;
static const unsigned int TARGET_LENGTH_MASK = ~TARGET_FRACTIONAL_MASK;
static const uint64 nFractionalDifficultyMax = (1llu << (nFractionalBits + 32));
static const uint64 nFractionalDifficultyMin = (1llu << 32);
static const uint64 nFractionalDifficultyThreshold = (1llu << (8 + 32));
static const unsigned int nWorkTransitionRatio = 32;
static const unsigned int nWorkTransitionRatioLog = 5; // log_2(32) = 5
unsigned int TargetGetLimit();
unsigned int TargetGetInitial();
unsigned int TargetGetLength(unsigned int nBits);
bool TargetSetLength(unsigned int nLength, unsigned int& nBits);
unsigned int TargetGetFractional(unsigned int nBits);
uint64 TargetGetFractionalDifficulty(unsigned int nBits);
bool TargetSetFractionalDifficulty(uint64 nFractionalDifficulty, unsigned int& nBits);
std::string TargetToString(unsigned int nBits);
unsigned int TargetFromInt(unsigned int nLength);
bool TargetGetMint(unsigned int nBits, uint64& nMint);
bool TargetGetNext(unsigned int nBits, int64 nInterval, int64 nTargetSpacing, int64 nActualSpacing, unsigned int& nBitsNext);

// Check prime proof-of-work
enum // prime chain type
{
    PRIME_CHAIN_CUNNINGHAM1 = 1u,
    PRIME_CHAIN_CUNNINGHAM2 = 2u,
    PRIME_CHAIN_BI_TWIN     = 3u,
};
bool CheckPrimeProofOfWork(uint256 hashBlockHeader, unsigned int nBits, const CBigNum& bnPrimeChainMultiplier, unsigned int& nChainType, unsigned int& nChainLength);

// prime target difficulty value for visualization
double GetPrimeDifficulty(unsigned int nBits);
// Estimate work transition target to longer prime chain
unsigned int EstimateWorkTransition(unsigned int nPrevWorkTransition, unsigned int nBits, unsigned int nChainLength);
// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength);
// primorial form of prime chain origin
std::string GetPrimeOriginPrimorialForm(CBigNum& bnPrimeChainOrigin);


/********************/
/* PRIMECOIN MINING */
/********************/

class CSieveOfEratosthenes;
class CPrimalityTestParams;

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(CBlock& block, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTests, unsigned int& nPrimesHit, mpz_class& mpzHash, CBlockIndex* pindexPrev, unsigned int vChainsFound[nMaxChainLength], CSieveOfEratosthenes& sieve, CPrimalityTestParams& testParams);

// Perform Fermat test with trial division
// Return values:
//   true  - passes trial division test and Fermat test; probable prime
//   false - failed either trial division or Fermat test; composite
bool ProbablePrimalityTestWithTrialDivision(const mpz_class& mpzCandidate, unsigned int nTrialDivisionLimit, CPrimalityTestParams& testParams);

// Estimate the probability of primality for a number in a candidate chain
double EstimateCandidatePrimeProbability(unsigned int nPrimorialMultiplier, unsigned int nChainPrimeNum, unsigned int nMiningProtocol);
// Esimate the prime probablity of numbers that haven't been sieved
double EstimateNormalPrimeProbability(unsigned int nPrimorialMultiplier, unsigned int nChainPrimeNum, unsigned int nMiningProtocol);

/*
 * Use GCC-style builtin functions such as
 * __builtin_popcountl and
 * __sync_add_and_fetch
 */
#if defined(__GNUC__) || defined(__clang__)
#    define USE_GCC_BUILTINS
#endif

#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__x86_64__) || defined(_M_X64)
#    define USE_ROTATE
#    if defined(__GNUC__) || defined(__clang__)
#        define USE_ASM
#    endif
#endif

// Check if the target platform is 64-bit
#if defined(__x86_64__) || defined(_M_X64)
typedef unsigned long long sieve_word_t;
#else
typedef unsigned int sieve_word_t;
#endif

#ifdef USE_ROTATE
inline sieve_word_t rotate_left(sieve_word_t bits, unsigned int count)
{
#ifdef USE_ASM
    asm("rol %1, %0"
        : "+r" (bits)
        : "c" ((uint8_t)count)
        : "cc");
    return bits;
#else
    // NOTE: At least LLVM doesn't always recognize this pattern
    return (bits << count) | (bits >> (nWordBits - count));
#endif
}
#endif

class CPrimalityTestParams
{
public:
    // GMP variables
    mpz_t mpzE;
    mpz_t mpzR;
    mpz_t mpzRplusOne;

    // GMP C++ variables
    mpz_class mpzHashFixedMult;
    mpz_class mpzChainOrigin;
    mpz_class mpzOriginMinusOne;
    mpz_class mpzOriginPlusOne;
    mpz_class mpzN;

    // Values specific to a round
    unsigned int nBits;
    unsigned int nCandidateType;

    // Results
    unsigned int nChainLength;

    CPrimalityTestParams()
    {
        nBits = 0;
        nCandidateType = 0;
        nChainLength = 0;
        mpz_init(mpzE);
        mpz_init(mpzR);
        mpz_init(mpzRplusOne);
    }

    ~CPrimalityTestParams()
    {
        mpz_clear(mpzE);
        mpz_clear(mpzR);
        mpz_clear(mpzRplusOne);
    }
};

// Sieve of Eratosthenes for proof-of-work mining
//
// Includes the sieve extension feature from jhPrimeminer by jh000
//
// A layer of the sieve determines whether the CC1 or CC2 chain members near the
// origin fixed_multiplier * candidate_multiplier * 2^k are known to be
// composites.
//
// The default sieve is composed of layers 1 .. nChainLength.
//
// An extension i is composed of layers i .. i + nChainLength. The candidates
// indexes from the extensions are multiplied by 2^i. The first half of the
// candidates are covered by the default sieve and previous extensions.
//
// The larger numbers in the extensions have a slightly smaller probability of
// being primes and take slightly longer to test but they can be calculated very
// efficiently because the layers overlap.
class CSieveOfEratosthenes
{
    unsigned int nSieveSize; // size of the sieve
    unsigned int nSieveFilterPrimes; // filter a certain number of primes
    unsigned int nSieveExtensions; // extend the sieve a given number of times
    unsigned int nBits; // target of the prime chain to search for
    mpz_class mpzHash; // hash of the block header
    mpz_class mpzFixedMultiplier; // fixed round multiplier
    mpz_class mpzHashFixedMult; // mpzHash * mpzFixedMultiplier

    // final set of candidates for probable primality checking
    sieve_word_t *vfCandidates;
    sieve_word_t *vfCompositeBiTwin;
    sieve_word_t *vfCompositeCunningham1;
    sieve_word_t *vfCompositeCunningham2;
    sieve_word_t *vfCompositeLayerCC1;
    sieve_word_t *vfCompositeLayerCC2;

    // extended sets
    sieve_word_t *vfExtendedCandidates;
    sieve_word_t *vfExtendedCompositeBiTwin;
    sieve_word_t *vfExtendedCompositeCunningham1;
    sieve_word_t *vfExtendedCompositeCunningham2;

    // divisible multipliers
    unsigned int *vCunningham1Multipliers;
    unsigned int *vCunningham2Multipliers;

    static const unsigned int nWordBits = 8 * sizeof(sieve_word_t);
    unsigned int nCandidatesWords;
    unsigned int nCandidatesBytes;

    unsigned int nPrimeSeq; // prime sequence number currently being processed
    unsigned int nCandidateCount; // cached total count of candidates
    unsigned int nCandidateMultiplier; // current candidate for power test
    unsigned int nCandidateIndex; // internal candidate index
    bool fCandidateIsExtended; // is the current candidate in the extended part
    unsigned int nCandidateActiveExtension; // which extension is active

    unsigned int nChainLength; // target chain length
    unsigned int nSieveLayers; // sieve layers
    unsigned int nPrimes; // number of times to weave the sieve
    unsigned int nL1CacheElements; // number of bits that can be stored in L1 cache

    CBlockIndex* pindexPrev;

    // previous parameters
    unsigned int nCandidatesBytesPrev;
    unsigned int nSieveExtensionsPrev;
    unsigned int nMultiplierBytesPrev;

    bool fIsReady;
    bool fIsDepleted;

    unsigned int GetWordNum(unsigned int nBitNum) {
        return nBitNum / nWordBits;
    }

    sieve_word_t GetBitMask(unsigned int nBitNum) {
        return (sieve_word_t)1 << (nBitNum % nWordBits);
    }

    void ProcessMultiplier(sieve_word_t *vfComposites, const unsigned int nMinMultiplier, const unsigned int nMaxMultiplier, const std::vector<unsigned int>& vPrimes, unsigned int *vMultipliers, unsigned int nLayerSeq);

    void freeArrays()
    {
        if (vfCandidates)
            free(vfCandidates);
        if (vfCompositeBiTwin)
            free(vfCompositeBiTwin);
        if (vfCompositeCunningham1)
            free(vfCompositeCunningham1);
        if (vfCompositeCunningham2)
            free(vfCompositeCunningham2);
        if (vfCompositeLayerCC1)
            free(vfCompositeLayerCC1);
        if (vfCompositeLayerCC2)
            free(vfCompositeLayerCC2);
        if (vfExtendedCandidates)
            free(vfExtendedCandidates);
        if (vfExtendedCompositeBiTwin)
            free(vfExtendedCompositeBiTwin);
        if (vfExtendedCompositeCunningham1)
            free(vfExtendedCompositeCunningham1);
        if (vfExtendedCompositeCunningham2)
            free(vfExtendedCompositeCunningham2);
        if (vCunningham1Multipliers)
            free(vCunningham1Multipliers);
        if (vCunningham2Multipliers)
            free(vCunningham2Multipliers);
        vfCandidates = NULL;
        vfCompositeBiTwin = NULL;
        vfCompositeCunningham1 = NULL;
        vfCompositeCunningham2 = NULL;
        vfCompositeLayerCC1 = NULL;
        vfCompositeLayerCC2 = NULL;
        vfExtendedCandidates = NULL;
        vfExtendedCompositeBiTwin = NULL;
        vfExtendedCompositeCunningham1 = NULL;
        vfExtendedCompositeCunningham2 = NULL;
    }

public:
    CSieveOfEratosthenes()
    {
        nSieveSize = 0;
        nSieveFilterPrimes = 0;
        nSieveExtensions = 0;
        nBits = 0;
        mpzHash = 0;
        mpzFixedMultiplier = 0;
        mpzHashFixedMult = 0;
        vfCandidates = NULL;
        vfCompositeBiTwin = NULL;
        vfCompositeCunningham1 = NULL;
        vfCompositeCunningham2 = NULL;
        vfCompositeLayerCC1 = NULL;
        vfCompositeLayerCC2 = NULL;
        vfExtendedCandidates = NULL;
        vfExtendedCompositeBiTwin = NULL;
        vfExtendedCompositeCunningham1 = NULL;
        vfExtendedCompositeCunningham2 = NULL;
        vCunningham1Multipliers = NULL;
        vCunningham2Multipliers = NULL;
        nCandidatesWords = 0;
        nCandidatesBytes = 0;
        nCandidatesBytesPrev = 0;
        nSieveExtensionsPrev = 0;
        nMultiplierBytesPrev = 0;
        nPrimeSeq = 0;
        nCandidateCount = 0;
        nCandidateMultiplier = 0;
        nCandidateIndex = 0;
        fCandidateIsExtended = false;
        nCandidateActiveExtension = 0;
        nChainLength = 0;
        nSieveLayers = 0;
        nPrimes = 0;
        nL1CacheElements = 0;
        pindexPrev = NULL;
        fIsReady = false;
        fIsDepleted = true;
    }

    ~CSieveOfEratosthenes()
    {
        freeArrays();
    }

    void Reset(unsigned int nSieveSize, unsigned int nSieveFilterPrimes, unsigned int nSieveExtensions, unsigned int nL1CacheSize, unsigned int nBits, mpz_class& mpzHash, mpz_class& mpzFixedMultiplier, CBlockIndex* pindexPrev)
    {
        this->nSieveSize = nSieveSize;
        this->nSieveFilterPrimes = nSieveFilterPrimes;
        this->nSieveExtensions = nSieveExtensions;
        nL1CacheElements = nL1CacheSize * 8;
        this->nBits = nBits;
        this->mpzHash = mpzHash;
        this->mpzFixedMultiplier = mpzFixedMultiplier;
        this->pindexPrev = pindexPrev;
        mpzHashFixedMult = mpzHash * mpzFixedMultiplier;
        nPrimeSeq = 0;
        nCandidateCount = 0;
        nCandidateMultiplier = 0;
        nCandidateIndex = 0;
        fCandidateIsExtended = false;
        nCandidateActiveExtension = 0;
        nCandidatesWords = (nSieveSize + nWordBits - 1) / nWordBits;
        nCandidatesBytes = nCandidatesWords * sizeof(sieve_word_t);
        nChainLength = TargetGetLength(nBits);

        // Override target length if requested
        if (nSieveTargetLength > 0)
            nChainLength = nSieveTargetLength;
        nSieveLayers = nChainLength + nSieveExtensions;

        // Filter only a certain number of prime factors
        // Most composites are still found
        nPrimes = nSieveFilterPrimes;
        const unsigned int nMultiplierBytes = nPrimes * nSieveLayers * sizeof(unsigned int);

        // Allocate arrays if parameters have changed
        if (nCandidatesBytes != nCandidatesBytesPrev || nSieveExtensions != nSieveExtensionsPrev || nMultiplierBytes != nMultiplierBytesPrev)
        {
            nCandidatesBytesPrev = nCandidatesBytes;
            nSieveExtensionsPrev = nSieveExtensions;
            nMultiplierBytesPrev = nMultiplierBytes;
            freeArrays();
            vfCandidates = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeBiTwin = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeCunningham1 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeCunningham2 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeLayerCC1 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfCompositeLayerCC2 = (sieve_word_t *)malloc(nCandidatesBytes);
            vfExtendedCandidates = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeBiTwin = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeCunningham1 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vfExtendedCompositeCunningham2 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes);
            vCunningham1Multipliers = (unsigned int *)malloc(nMultiplierBytes);
            vCunningham2Multipliers = (unsigned int *)malloc(nMultiplierBytes);
        }

        // Initialize arrays
        memset(vfCandidates, 0, nCandidatesBytes);
        memset(vfCompositeBiTwin, 0, nCandidatesBytes);
        memset(vfCompositeCunningham1, 0, nCandidatesBytes);
        memset(vfCompositeCunningham2, 0, nCandidatesBytes);
        memset(vfCompositeLayerCC1, 0, nCandidatesBytes);
        memset(vfCompositeLayerCC2, 0, nCandidatesBytes);
        memset(vfExtendedCandidates, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeBiTwin, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeCunningham1, 0, nSieveExtensions * nCandidatesBytes);
        memset(vfExtendedCompositeCunningham2, 0, nSieveExtensions * nCandidatesBytes);
        memset(vCunningham1Multipliers, 0xFF, nMultiplierBytes);
        memset(vCunningham2Multipliers, 0xFF, nMultiplierBytes);

        fIsReady = true;
        fIsDepleted = false;
    }

    // Get total number of candidates for power test
    unsigned int GetCandidateCount()
    {
        if (nCandidateCount)
            return nCandidateCount;

        unsigned int nCandidates = 0;
#ifdef USE_GCC_BUILTINS
        for (unsigned int i = 0; i < nCandidatesWords; i++)
            nCandidates += __builtin_popcountl(vfCandidates[i]);
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = nCandidatesWords / 2; i < nCandidatesWords; i++)
                nCandidates += __builtin_popcountl(vfExtendedCandidates[j * nCandidatesWords + i]);
#else
        for (unsigned int i = 0; i < nCandidatesWords; i++)
        {
            sieve_word_t lBits = vfCandidates[i];
            for (unsigned int j = 0; j < nWordBits; j++)
            {
                nCandidates += (lBits & 1);
                lBits >>= 1;
            }
        }
        for (unsigned int j = 0; j < nSieveExtensions; j++)
        {
            for (unsigned int i = nCandidatesWords / 2; i < nCandidatesWords; i++)
            {
                sieve_word_t lBits = vfExtendedCandidates[j * nCandidatesWords + i];
                for (unsigned int j = 0; j < nWordBits; j++)
                {
                    nCandidates += (lBits & 1);
                    lBits >>= 1;
                }
            }
        }
#endif
        nCandidateCount = nCandidates;
        return nCandidates;
    }

    // Scan for the next candidate multiplier (variable part)
    // Return values:
    //   True - found next candidate; nVariableMultiplier has the candidate
    //   False - scan complete, no more candidate and reset scan
    bool GetNextCandidateMultiplier(unsigned int& nVariableMultiplier, unsigned int& nCandidateType)
    {
        sieve_word_t *vfActiveCandidates;
        sieve_word_t *vfActiveCompositeTWN;
        sieve_word_t *vfActiveCompositeCC1;
        sieve_word_t *vfActiveCompositeCC2;

        if (fCandidateIsExtended)
        {
            vfActiveCandidates = vfExtendedCandidates + nCandidateActiveExtension * nCandidatesWords;
            vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nCandidateActiveExtension * nCandidatesWords;
            vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nCandidateActiveExtension * nCandidatesWords;
            vfActiveCompositeCC2 = vfExtendedCompositeCunningham2 + nCandidateActiveExtension * nCandidatesWords;
        }
        else
        {
            vfActiveCandidates = vfCandidates;
            vfActiveCompositeTWN = vfCompositeBiTwin;
            vfActiveCompositeCC1 = vfCompositeCunningham1;
            vfActiveCompositeCC2 = vfCompositeCunningham2;
        }

        // Acquire the current word from the bitmap
        sieve_word_t lBits = vfActiveCandidates[GetWordNum(nCandidateIndex)];

        loop
        {
            nCandidateIndex++;
            if (nCandidateIndex >= nSieveSize)
            {
                // Check if extensions are available
                if (!fCandidateIsExtended && nSieveExtensions > 0)
                {
                    fCandidateIsExtended = true;
                    nCandidateActiveExtension = 0;
                    nCandidateIndex = nSieveSize / 2;
                }
                else if (fCandidateIsExtended && nCandidateActiveExtension + 1 < nSieveExtensions)
                {
                    nCandidateActiveExtension++;
                    nCandidateIndex = nSieveSize / 2;
                }
                else
                {
                    // Out of candidates
                    fCandidateIsExtended = false;
                    nCandidateActiveExtension = 0;
                    nCandidateIndex = 0;
                    nCandidateMultiplier = 0;
                    fIsDepleted = true;
                    return false;
                }

                // Fix the pointers
                if (fCandidateIsExtended)
                {
                    vfActiveCandidates = vfExtendedCandidates + nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCompositeCC2 = vfExtendedCompositeCunningham2 + nCandidateActiveExtension * nCandidatesWords;
                }
                else
                {
                    vfActiveCandidates = vfCandidates;
                    vfActiveCompositeTWN = vfCompositeBiTwin;
                    vfActiveCompositeCC1 = vfCompositeCunningham1;
                    vfActiveCompositeCC2 = vfCompositeCunningham2;
                }

                // Update the current word
                lBits = vfActiveCandidates[GetWordNum(nCandidateIndex)];
            }

            if (nCandidateIndex % nWordBits == 0)
            {
                // Update the current word
                lBits = vfActiveCandidates[GetWordNum(nCandidateIndex)];

                // Check if any bits are set
                if (lBits == 0)
                {
                    // Skip an entire word
                    nCandidateIndex += nWordBits - 1;
                    continue;
                }
            }

            if (lBits & GetBitMask(nCandidateIndex))
            {
                if (fCandidateIsExtended)
                    nCandidateMultiplier = nCandidateIndex * (2 << nCandidateActiveExtension);
                else
                    nCandidateMultiplier = nCandidateIndex;
                nVariableMultiplier = nCandidateMultiplier;
                if (~vfActiveCompositeTWN[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex))
                    nCandidateType = PRIME_CHAIN_BI_TWIN;
                else if (~vfActiveCompositeCC1[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex))
                    nCandidateType = PRIME_CHAIN_CUNNINGHAM1;
                else if (~vfActiveCompositeCC2[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex))
                    nCandidateType = PRIME_CHAIN_CUNNINGHAM2;
                else
                    nCandidateType = 0; // unknown
                return true;
            }
        }
    }

    // Get progress percentage of the sieve
    unsigned int GetProgressPercentage();

    // Weave the sieve for the next prime in table
    // Return values:
    //   True  - weaved another prime; nComposite - number of composites removed
    //   False - sieve already completed
    bool Weave();

    bool IsReady() { return fIsReady; }
    bool IsDepleted() { return fIsDepleted; }
    void Deplete() { fIsDepleted = true; }
};

inline void mpz_set_uint256(mpz_t r, uint256& u)
{
    mpz_import(r, 32 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, &u);
}

#endif
