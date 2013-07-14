// Copyright (c) 2013 Primecoin developers
// Distributed under the MIT/X11 software license, see the accompanying
// file COPYING or http://www.opensource.org/licenses/mit-license.php.

#ifndef PRIMECOIN_PRIME_H
#define PRIMECOIN_PRIME_H

#include "main.h"

#include <gmp.h>
#include <gmpxx.h>
#include <bitset>

static const unsigned int nMaxSieveSize = 10000000u;
static const unsigned int nDefaultSieveSize = 1000000u;
static const unsigned int nMinSieveSize = 100000u;
static const uint256 hashBlockHeaderLimit = (uint256(1) << 255);
static const CBigNum bnOne = 1;
static const CBigNum bnPrimeMax = (bnOne << 2000) - 1;
static const CBigNum bnPrimeMin = (bnOne << 255);
static const mpz_class mpzOne = 1;
static const mpz_class mpzTwo = 2;
static const mpz_class mpzPrimeMax = (mpzOne << 2000) - 1;
static const mpz_class mpzPrimeMin = (mpzOne << 255);

// Estimate how many 5-chains are found per hour
static const unsigned int nStatsChainLength = 5;

extern unsigned int nTargetInitialLength;
extern unsigned int nTargetMinLength;

// Generate small prime table
void GeneratePrimeTable();
// Get next prime number of p
bool PrimeTableGetNextPrime(unsigned int& p);
// Get previous prime number of p
bool PrimeTableGetPreviousPrime(unsigned int& p);

// Compute primorial number p#
void Primorial(unsigned int p, mpz_class& mpzPrimorial);
// Compute the first primorial number greater than or equal to bn
void PrimorialAt(mpz_class& bn, mpz_class& mpzPrimorial);

// Test probable prime chain for: bnPrimeChainOrigin
// fFermatTest
//   true - Use only Fermat tests
//   false - Use Fermat-Euler-Lagrange-Lifchitz tests
// Return value:
//   true - Probable prime chain found (one of nChainLength meeting target)
//   false - prime chain too short (none of nChainLength meeting target)
bool ProbablePrimeChainTest(const mpz_class& mpzPrimeChainOrigin, unsigned int nBits, bool fFermatTest, unsigned int& nChainLengthCunningham1, unsigned int& nChainLengthCunningham2, unsigned int& nChainLengthBiTwin);

static const unsigned int nFractionalBits = 24;
static const unsigned int TARGET_FRACTIONAL_MASK = (1u<<nFractionalBits) - 1;
static const unsigned int TARGET_LENGTH_MASK = ~TARGET_FRACTIONAL_MASK;
static const uint64 nFractionalDifficultyMax = (1llu << (nFractionalBits + 32));
static const uint64 nFractionalDifficultyMin = (1llu << 32);
static const uint64 nFractionalDifficultyThreshold = (1llu << (8 + 32));
static const unsigned int nWorkTransitionRatio = 32;
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

// Mine probable prime chain of form: n = h * p# +/- 1
bool MineProbablePrimeChain(CBlock& block, mpz_class& mpzFixedMultiplier, bool& fNewBlock, unsigned int& nTriedMultiplier, unsigned int& nProbableChainLength, unsigned int& nTests, unsigned int& nPrimesHit, unsigned int& nChainsHit, mpz_class& mpzHash);

// Check prime proof-of-work
enum // prime chain type
{
    PRIME_CHAIN_CUNNINGHAM1 = 1u,
    PRIME_CHAIN_CUNNINGHAM2 = 2u,
    PRIME_CHAIN_BI_TWIN     = 3u,
};
bool CheckPrimeProofOfWork(uint256& hashBlockHeader, unsigned int nBits, const mpz_class& mpzPrimeChainMultiplier, unsigned int& nChainType, unsigned int& nChainLength);

// prime target difficulty value for visualization
double GetPrimeDifficulty(unsigned int nBits);
// Estimate work transition target to longer prime chain
unsigned int EstimateWorkTransition(unsigned int nPrevWorkTransition, unsigned int nBits, unsigned int nChainLength);
// prime chain type and length value
std::string GetPrimeChainName(unsigned int nChainType, unsigned int nChainLength);

// Sieve of Eratosthenes for proof-of-work mining
class CSieveOfEratosthenes
{
    unsigned int nSieveSize; // size of the sieve
    unsigned int nBits; // target of the prime chain to search for
    mpz_class mpzFixedFactor; // fixed factor to derive the chain

    // bitmaps of the sieve, index represents the variable part of multiplier
    std::bitset<nMaxSieveSize> vfCompositeCunningham1;
    std::bitset<nMaxSieveSize> vfCompositeCunningham2;
    std::bitset<nMaxSieveSize> vfCompositeBiTwin;
    std::bitset<nMaxSieveSize> vfCandidates;

    unsigned int nPrimeSeq; // prime sequence number currently being processed
    unsigned int nCandidateMultiplier; // current candidate for power test
    
    CBlockIndex* pindexPrev;

public:
    CSieveOfEratosthenes(unsigned int nSieveSize, unsigned int nBits, mpz_class& mpzHash, mpz_class& mpzFixedMultiplier, CBlockIndex* pindexPrev)
    {
        this->nSieveSize = nSieveSize;
        this->nBits = nBits;
        this->mpzFixedFactor = mpzFixedMultiplier * mpzHash;
        this->pindexPrev = pindexPrev;
        nPrimeSeq = 0;
        nCandidateMultiplier = 0;
    }
    
    void CombineCandidates()
    {
        //vfCandidates = ~(vfCompositeCunningham1 & vfCompositeCunningham2 & vfCompositeBiTwin);

        // Faster version (not exactly clean but much faster)
        unsigned char *cCandidates = (unsigned char *)&vfCandidates;
        unsigned char *cCompositeCunningham1 = (unsigned char *)&vfCompositeCunningham1;
        unsigned char *cCompositeCunningham2 = (unsigned char *)&vfCompositeCunningham2;
        unsigned char *cCompositeBiTwin = (unsigned char *)&vfCompositeBiTwin;
        const unsigned int nBytes = nSieveSize / 8;

        unsigned long *lCandidates = (unsigned long *)cCandidates;
        unsigned long *lCompositeCunningham1 = (unsigned long *)cCompositeCunningham1;
        unsigned long *lCompositeCunningham2 = (unsigned long *)cCompositeCunningham2;
        unsigned long *lCompositeBiTwin = (unsigned long *)cCompositeBiTwin;
        const unsigned int nLongs = nBytes / sizeof(unsigned long);
        for (unsigned int i = 0; i < nLongs; i++) {
            lCandidates[i] = ~(lCompositeCunningham1[i] & lCompositeCunningham2[i] & lCompositeBiTwin[i]);
        }
        const unsigned int nBytesProcessed = sizeof(unsigned long) * nLongs;

        for (unsigned int i = nBytesProcessed; i < nBytes; i++) {
            cCandidates[i] = ~(cCompositeCunningham1[i] & cCompositeCunningham2[i] & cCompositeBiTwin[i]);
        }
    }

    // Get total number of candidates for power test
    unsigned int GetCandidateCount()
    {
        unsigned int nCandidates = 0;
        for (unsigned int nMultiplier = 0; nMultiplier < nSieveSize; nMultiplier++)
        {
            if (vfCandidates[nMultiplier])
                nCandidates++;
        }
        return nCandidates;
    }

    // Scan for the next candidate multiplier (variable part)
    // Return values:
    //   True - found next candidate; nVariableMultiplier has the candidate
    //   False - scan complete, no more candidate and reset scan
    bool GetNextCandidateMultiplier(unsigned int& nVariableMultiplier)
    {
        loop
        {
            nCandidateMultiplier++;
            if (nCandidateMultiplier >= nSieveSize)
            {
                nCandidateMultiplier = 0;
                return false;
            }
            if (vfCandidates[nCandidateMultiplier])
            {
                nVariableMultiplier = nCandidateMultiplier;
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
};

inline void mpz_set_uint256(mpz_t r, uint256& u)
{
    mpz_import(r, 32 / sizeof(unsigned long), -1, sizeof(unsigned long), -1, 0, &u);
}

#endif
