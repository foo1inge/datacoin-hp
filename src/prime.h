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
static const unsigned int nDefaultSieveExtensions = 10;
static const unsigned int nDefaultSieveExtensionsTestnet = 4;
extern unsigned int nSieveExtensions;
static const unsigned int nMaxSieveFilterPrimes = 78498u; // size of prime table
static const unsigned int nDefaultSieveFilterPrimes = 14000u;
static const unsigned int nMinSieveFilterPrimes = 1000u;
extern unsigned int nSieveFilterPrimes;
static const unsigned int nMaxSieveSize = 10000000u;
static const unsigned int nDefaultSieveSize = 1376256u;
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
static const std::string strDefaultDonationPercentage = "1.0";
static const double dMinDonationPercentage = 0.1;
static const double dMaxDonationPercentage = 99.9;
static const std::string strDefaultDonationAddress = "DSSHsB1R8mrZd1ujhxcPqQaqAu2cNZsCNn";
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

#ifdef USE_GCC_BUILTINS
#    define likely(x)       __builtin_expect(!!(x),1)
#    define unlikely(x)     __builtin_expect(!!(x),0)
#else
#    define likely(x)       (x)
#    define unlikely(x)     (x)
#endif

#if defined(__i386__) || defined(_M_IX86) || defined(_X86_) || defined(__x86_64__) || defined(_M_X64)
#    define USE_ROTATE
#    if defined(__GNUC__) || defined(__clang__)
#        define USE_ASM
#        define USE_INTRINSICS
#        define USE_BTS
#    endif
#endif

// Unrolled loops only on x86-64
// Probably too much register pressure on x86
#if defined(__x86_64__) || defined(_M_X64)
#    define USE_UNROLLED_LOOPS
#endif

#ifdef USE_INTRINSICS
#    include <immintrin.h>
#endif

// Check if the target platform is 64-bit
#if defined(__x86_64__) || defined(_M_X64)
#define USE_64BIT
typedef unsigned long long sieve_word_t;
typedef signed long long sieve_signed_t;
#else
typedef unsigned int sieve_word_t;
typedef signed int sieve_signed_t;
#endif

static const unsigned int nWordBits = 8 * sizeof(sieve_word_t);
static const unsigned int nWordBitsLog = (nWordBits == 64) ? 6 : 5;

#ifdef USE_BTS
// The bts instruction should be very fast starting from Intel Core 2
inline sieve_word_t bts(sieve_word_t nr, sieve_word_t bits)
{
    asm ("bts %1, %0"
                    : "+r" (bits)
                    : "Ir" (nr)
                    : "cc");
    return bits;
}
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

#if defined(USE_ASM) && defined(__BMI2__)
#   define USE_BMI2
#endif

#if defined(USE_INTRINSICS) && defined(__SSE2__)
#    define USE_SSE2
#endif

#if defined(USE_INTRINSICS) && defined(__AVX2__)
#    define USE_AVX2
#endif

#if defined(USE_AVX2)
const unsigned int nRequiredAlignment = 256;
#elif defined(USE_SSE2)
const unsigned int nRequiredAlignment = 128;
#elif defined(USE_64BIT)
const unsigned int nRequiredAlignment = 64;
#else
const unsigned int nRequiredAlignment = 32;
#endif

#ifdef USE_BMI2
inline sieve_word_t shrx(sieve_word_t bits, sieve_word_t count)
{
    sieve_word_t result;
    asm("shrx %2, %1, %0"
        :"=r" (result)
        :"r" (bits), "r" (count));
    return result;
}

inline sieve_word_t shlx(sieve_word_t bits, sieve_word_t count)
{
    sieve_word_t result;
    asm("shlx %2, %1, %0"
        :"=r" (result)
        :"r" (bits), "r" (count));
    return result;
}
#endif

class CPrimalityTestParams
{
public:
    // GMP C++ variables
    mpz_class mpzHashFixedMult;
    mpz_class mpzChainOrigin;
    mpz_class mpzOriginMinusOne;
    mpz_class mpzOriginPlusOne;
    mpz_class mpzN;
    mpz_class mpzNMinusOne;
    mpz_class mpzBase;
    mpz_class mpzR;
    mpz_class mpzR2;
    mpz_class mpzE;
    mpz_class mpzFrac;

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
    }
};

// Sieve of Eratosthenes for proof-of-work mining
//
// Includes the sieve extension feature from jhPrimeminer by jh000
//
// Layer k of the sieve determines whether numbers of the form
// hash * primorial * candidate * 2^k - 1 and
// hash * primorial * candidate * 2^k + 1
// are known to be composites.
//
// The default sieve is composed of layers 1 .. nChainLength. The multipliers
// in the default sieve are all odd.
//
// An extension i is composed of layers i .. i + nChainLength. The candidates
// indexes from the extensions are multiplied by 2^i. Therefore, the resulting
// multipliers are always even.
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

    // raw and possibly unaligned pointers
    sieve_word_t *vfRawCandidates;
    sieve_word_t *vfRawCompositeBiTwin;
    sieve_word_t *vfRawCompositeCunningham1;
    sieve_word_t *vfRawCompositeCunningham2;
    sieve_word_t *vfRawCompositeLayerCC1;
    sieve_word_t *vfRawCompositeLayerCC2;
    sieve_word_t *vfRawExtendedCandidates;
    sieve_word_t *vfRawExtendedCompositeBiTwin;
    sieve_word_t *vfRawExtendedCompositeCunningham1;
    sieve_word_t *vfRawExtendedCompositeCunningham2;

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
    unsigned int nMinPrimeSeq; // smallest prime which will be used for sieving

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

    void ProcessMultiplier(sieve_word_t *vfComposites, const unsigned int nMinMultiplier, const unsigned int nMaxMultiplier, unsigned int *vMultipliers, unsigned int nLayerSeq);

    void freeArrays()
    {
        if (vfRawCandidates)
            free(vfRawCandidates);
        if (vfRawCompositeBiTwin)
            free(vfRawCompositeBiTwin);
        if (vfRawCompositeCunningham1)
            free(vfRawCompositeCunningham1);
        if (vfRawCompositeCunningham2)
            free(vfRawCompositeCunningham2);
        if (vfRawCompositeLayerCC1)
            free(vfRawCompositeLayerCC1);
        if (vfRawCompositeLayerCC2)
            free(vfRawCompositeLayerCC2);
        if (vfRawExtendedCandidates)
            free(vfRawExtendedCandidates);
        if (vfRawExtendedCompositeBiTwin)
            free(vfRawExtendedCompositeBiTwin);
        if (vfRawExtendedCompositeCunningham1)
            free(vfRawExtendedCompositeCunningham1);
        if (vfRawExtendedCompositeCunningham2)
            free(vfRawExtendedCompositeCunningham2);
        if (vCunningham1Multipliers)
            free(vCunningham1Multipliers);
        if (vCunningham2Multipliers)
            free(vCunningham2Multipliers);
        vfRawCandidates = NULL;
        vfRawCompositeBiTwin = NULL;
        vfRawCompositeCunningham1 = NULL;
        vfRawCompositeCunningham2 = NULL;
        vfRawCompositeLayerCC1 = NULL;
        vfRawCompositeLayerCC2 = NULL;
        vfRawExtendedCandidates = NULL;
        vfRawExtendedCompositeBiTwin = NULL;
        vfRawExtendedCompositeCunningham1 = NULL;
        vfRawExtendedCompositeCunningham2 = NULL;
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
        vfRawCandidates = NULL;
        vfRawCompositeBiTwin = NULL;
        vfRawCompositeCunningham1 = NULL;
        vfRawCompositeCunningham2 = NULL;
        vfRawCompositeLayerCC1 = NULL;
        vfRawCompositeLayerCC2 = NULL;
        vfRawExtendedCandidates = NULL;
        vfRawExtendedCompositeBiTwin = NULL;
        vfRawExtendedCompositeCunningham1 = NULL;
        vfRawExtendedCompositeCunningham2 = NULL;
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
        nMinPrimeSeq = 0;
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
        nCandidatesWords = (nSieveSize + nRequiredAlignment - 1) / nRequiredAlignment * (nRequiredAlignment / nWordBits);
        nCandidatesBytes = nCandidatesWords * sizeof(sieve_word_t);
        nChainLength = TargetGetLength(nBits);
        nMinPrimeSeq = 0;

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
            vfRawCandidates = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawCompositeBiTwin = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawCompositeCunningham1 = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawCompositeCunningham2 = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawCompositeLayerCC1 = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawCompositeLayerCC2 = (sieve_word_t *)malloc(nCandidatesBytes + nRequiredAlignment);
            vfRawExtendedCandidates = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes + nRequiredAlignment);
            vfRawExtendedCompositeBiTwin = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes + nRequiredAlignment);
            vfRawExtendedCompositeCunningham1 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes + nRequiredAlignment);
            vfRawExtendedCompositeCunningham2 = (sieve_word_t *)malloc(nSieveExtensions * nCandidatesBytes + nRequiredAlignment);
            vCunningham1Multipliers = (unsigned int *)malloc(nMultiplierBytes);
            vCunningham2Multipliers = (unsigned int *)malloc(nMultiplierBytes);
#define ALIGN_PTR(x)    ((void *)(((uintptr_t)(x) + nRequiredAlignment - 1) / nRequiredAlignment * nRequiredAlignment))
            vfCandidates = (sieve_word_t *)ALIGN_PTR(vfRawCandidates);
            vfCompositeBiTwin = (sieve_word_t *)ALIGN_PTR(vfRawCompositeBiTwin);
            vfCompositeCunningham1 = (sieve_word_t *)ALIGN_PTR(vfRawCompositeCunningham1);
            vfCompositeCunningham2 = (sieve_word_t *)ALIGN_PTR(vfRawCompositeCunningham2);
            vfCompositeLayerCC1 = (sieve_word_t *)ALIGN_PTR(vfRawCompositeLayerCC1);
            vfCompositeLayerCC2 = (sieve_word_t *)ALIGN_PTR(vfRawCompositeLayerCC2);
            vfExtendedCandidates = (sieve_word_t *)ALIGN_PTR(vfRawExtendedCandidates);
            vfExtendedCompositeBiTwin = (sieve_word_t *)ALIGN_PTR(vfRawExtendedCompositeBiTwin);
            vfExtendedCompositeCunningham1 = (sieve_word_t *)ALIGN_PTR(vfRawExtendedCompositeCunningham1);
            vfExtendedCompositeCunningham2 = (sieve_word_t *)ALIGN_PTR(vfRawExtendedCompositeCunningham2);
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
#    ifdef USE_64BIT
        for (unsigned int i = 0; i < nCandidatesWords; i++)
            nCandidates += __builtin_popcountll(vfCandidates[i]);
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = 0; i < nCandidatesWords; i++)
                nCandidates += __builtin_popcountll(vfExtendedCandidates[j * nCandidatesWords + i]);
#    else
        for (unsigned int i = 0; i < nCandidatesWords; i++)
            nCandidates += __builtin_popcount(vfCandidates[i]);
        for (unsigned int j = 0; j < nSieveExtensions; j++)
            for (unsigned int i = 0; i < nCandidatesWords; i++)
                nCandidates += __builtin_popcount(vfExtendedCandidates[j * nCandidatesWords + i]);
#    endif
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
            for (unsigned int i = 0; i < nCandidatesWords; i++)
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
        const sieve_word_t nMaxWord = (nSieveSize + nWordBits - 1) / nWordBits;
        const unsigned int _nSieveSize = nSieveSize;

        if (fCandidateIsExtended)
        {
            const sieve_word_t nExtOffset = nCandidateActiveExtension * nCandidatesWords;
            vfActiveCandidates = vfExtendedCandidates + nExtOffset;
            vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nExtOffset;
            vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nExtOffset;
            vfActiveCompositeCC2 = vfExtendedCompositeCunningham2 + nExtOffset;
        }
        else
        {
            vfActiveCandidates = vfCandidates;
            vfActiveCompositeTWN = vfCompositeBiTwin;
            vfActiveCompositeCC1 = vfCompositeCunningham1;
            vfActiveCompositeCC2 = vfCompositeCunningham2;
        }

        loop
        {
            nCandidateIndex++;
            if (unlikely(nCandidateIndex >= _nSieveSize))
            {
                // Check if extensions are available
                if (!fCandidateIsExtended && nSieveExtensions > 0)
                {
                    fCandidateIsExtended = true;
                    nCandidateActiveExtension = 0;
                    nCandidateIndex = 0;
                }
                else if (fCandidateIsExtended && nCandidateActiveExtension + 1 < nSieveExtensions)
                {
                    nCandidateActiveExtension++;
                    nCandidateIndex = 0;
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
                    const sieve_word_t nExtOffset = nCandidateActiveExtension * nCandidatesWords;
                    vfActiveCandidates = vfExtendedCandidates + nExtOffset;
                    vfActiveCompositeTWN = vfExtendedCompositeBiTwin + nExtOffset;
                    vfActiveCompositeCC1 = vfExtendedCompositeCunningham1 + nExtOffset;
                    vfActiveCompositeCC2 = vfExtendedCompositeCunningham2 + nExtOffset;
                }
                else
                {
                    vfActiveCandidates = vfCandidates;
                    vfActiveCompositeTWN = vfCompositeBiTwin;
                    vfActiveCompositeCC1 = vfCompositeCunningham1;
                    vfActiveCompositeCC2 = vfCompositeCunningham2;
                }
            }

            if (unlikely(nCandidateIndex % nWordBits == 0))
            {
                sieve_word_t nWord = nCandidateIndex / nWordBits;
                // Fast loop to skip empty words
                for (; likely(nWord < nMaxWord); nWord++)
                {
                    if (vfActiveCandidates[nWord] != 0)
                        break;
                }
                nCandidateIndex = nWord * nWordBits;
                if (nCandidateIndex >= _nSieveSize)
                    continue;
            }

            if (unlikely((vfActiveCandidates[GetWordNum(nCandidateIndex)] & GetBitMask(nCandidateIndex)) != 0))
            {
                if (fCandidateIsExtended)
                    nCandidateMultiplier = (2 * nCandidateIndex + 1) * (2 << nCandidateActiveExtension);
                else
                    nCandidateMultiplier = 2 * nCandidateIndex + 1;
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
