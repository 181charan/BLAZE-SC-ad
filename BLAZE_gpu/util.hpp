#ifndef UTIL_HPP // Header guard
#define UTIL_HPP

#include <unordered_map>
#include <utility>
#include <cstdio>

#include <boost/functional/hash.hpp>

// Declarations here
constexpr unsigned kAlphabetSize = 28;
constexpr unsigned kmerSize = 3;
constexpr unsigned table = 21952; // Alphabet^3 28^3

int matrix[28][28] = {{-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768},
{-32768,4,-2,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-1,-2,-1,0,-4,-1,-1},
{-32768,-2,4,-3,4,1,-3,-1,0,-3,0,-4,-3,4,-2,0,-1,0,-1,-3,-4,-1,-3,0,-3,-4,-1,-3},
{-32768,0,-3,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-1,-2,-3,9,-4,-1,-1},
{-32768,-2,4,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-1,-3,1,-3,-4,-1,-3},
{-32768,-1,1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-1,-2,4,-4,-4,-1,-3},
{-32768,-2,-3,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,-1,3,-3,-2,-4,-1,0},
{-32768,0,-1,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-1,-3,-2,-3,-4,-1,-4},
{-32768,-2,0,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,-1,2,0,-3,-4,-1,-3},
{-32768,-1,-3,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,-1,-3,-1,-4,-1,3},
{-32768,-1,0,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-1,-2,1,-3,-4,-1,-3},
{-32768,-1,-4,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,-1,-3,-1,-4,-1,3},
{-32768,-1,-3,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,-1,-1,-1,-4,-1,2},
{-32768,-2,4,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-1,-2,0,-3,-4,-1,-3},
{-32768,-1,-2,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-1,-3,-1,-3,-4,-1,-3},
{-32768,-1,0,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,-1,4,-3,-4,-1,-2},
{-32768,-1,-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-1,-2,0,-3,-4,-1,-2},
{-32768,1,0,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-1,-2,0,-1,-4,-1,-2},
{-32768,0,-1,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-1,-2,-1,-1,-4,-1,-1},
{-32768,0,-3,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,-1,-2,-1,-4,-1,2},
{-32768,-3,-4,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,-1,2,-2,-2,-4,-1,-2},
{-32768,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-1,-1},
{-32768,-2,-3,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,-1,7,-2,-2,-4,-1,-1},
{-32768,-1,0,-3,1,4,-3,-2,0,-3,1,-3,-1,0,-1,4,0,0,-1,-2,-2,-1,-2,4,-3,-4,-1,-3},
{-32768,0,-3,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-1,-2,-3,9,-4,-1,-1},
{-32768,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1,-4,-4},
{-32768,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-1,-1},
{-32768,-1,-3,-1,-3,-3,0,-4,-3,3,-3,3,2,-3,-3,-2,-2,-2,-1,2,-2,-1,-1,-3,-1,-4,-1,3}};

int indx_calc(int a, int b, int c){
    //return(((a << 5*2) | (b << 5*1) | (c << 5*0)));
    return ((a * pow(28,2)) + (b * pow(28,1)) + (c*(pow(28,0))));
}

std::vector<int> neighbouring_words(int x, int y, int z){

    std::vector<int> neighbors;
    int score = 0;

    for(int a=0;a<28;a++){
        for(int b=0;b<28;b++){
            for(int c=0;c<28;c++){
                score = int(matrix[x][a] + matrix[y][b] + matrix[z][c]);
                if ( score >= 11  ){
                    if (x == a and y == b and z == c){
                        continue;
                    }
                    neighbors.push_back(indx_calc(a,b,c));
                    //neighbors.push_back(((a << 5*2) | (b << 5*1) | (c << 5*0)));
                }
            }
        }
    }

    return neighbors;
}


struct cellType {
    unsigned a;
    unsigned b;
    unsigned c;
};

std::vector<cellType> lut_generator(){
    std::vector<cellType> lut(21952);

    for(int a=0;a<28;a++){
        for(int b=0;b<28;b++){
            for(int c=0;c<28;c++){
                uint idx = indx_calc(a,b,c);
                lut[idx].a = a;
                lut[idx].b = b;
                lut[idx].c = c;
            }
        }
    }

    return lut;

}

#define T_PER_PROB 64

// #define TQ_VALUES  \
//   X(32) X(64) X(96) X(128) X(160) X(192) X(224) X(256) \
//   X(288) X(320) X(352) X(384) X(416) X(448) X(480) X(512) \
//   X(544) X(576) X(608) X(640) X(672) X(704) X(736) X(768) \
//   X(800) X(832) X(864) X(896) X(928) X(960) X(992) X(1024) \
//   X(1056) X(1088) X(1120) X(1152) X(1184) X(1216) X(1248) X(1280) \
//   X(1312) X(1344) X(1376) X(1408) X(1440) X(1472) X(1504) X(1536) \
//   X(1568) X(1600) X(1632) X(1664) X(1696) X(1728) X(1760) X(1792) \
//   X(1824) X(1856) X(1888) X(1920) X(1952) X(1984) X(2016) X(2048)

  #define TQ_VALUES  \
  X(32) X(64) X(96) X(128) X(160) X(192) X(224) X(256) \
  X(288) X(320) X(352) X(384) X(416) X(448) X(480) X(512) \
  X(544) X(576) X(608) X(640) 

// #define TQ_VALUES  \
//   X(640)

// #define SUBJECT_LEN_BINS \
//   Y(64) Y(128) Y(192) Y(256) Y(320) \
//   Y(384) Y(448) Y(512) Y(576) Y(1024) Y(2048)

#define SUBJECT_LEN_BINS \
  Y(64) Y(128) Y(192) Y(256) Y(320) \
  Y(384) Y(448) Y(512) Y(576) Y(1024) Y(2048)

template<unsigned int T_QUERY_BITSET_SIZE, unsigned int T_SUBJECT_LEN, unsigned int T_T_PER_PROB>
__global__ void gpuBLAZE_full(uint8_t * gpuDbArray, uint32_t * gpuSequenceLength, uint8_t * gpuQuerydbArray, char * gpuQuery, size_t gpuQueryLength, uint32_t streamID, uint8_t * gpuNumSeeds, int ungapped_cutoff, int gapped_cutoff, uint32_t * gappedSeeds, uint32_t * gpuHSPLen, uint16_t * query_lookup, uint16_t *query_prefix, uint16_t *seed_lookup_table);
template<unsigned int T_QUERY_BITSET_SIZE, unsigned int T_SUBJECT_LEN, unsigned int T_T_PER_PROB>
__global__ void gpuBLAZE_everything(uint8_t * gpuDbArray, uint32_t * gpuSequenceLength, uint8_t * gpuQuerydbArray, char * gpuQuery, size_t gpuQueryLength, uint32_t streamID, uint8_t * gpuNumSeeds, int ungapped_cutoff, int gapped_cutoff, uint32_t * gappedSeeds, uint32_t * gpuHSPLen);


#define INSTANTIATE(Q, S) \
  template __global__ void gpuBLAZE_full<Q, S, T_PER_PROB>(uint8_t*, uint32_t* , uint8_t*, char*, size_t, uint32_t, uint8_t*, int, int, uint32_t*, uint32_t*, uint16_t* , uint16_t*, uint16_t*); \
  template __global__ void gpuBLAZE_everything<Q, S, T_PER_PROB>(uint8_t*, uint32_t* , uint8_t*, char*, size_t, uint32_t, uint8_t*, int, int, uint32_t*, uint32_t*);

#define X(Q) \
  INSTANTIATE(Q, 64) \
  INSTANTIATE(Q, 128) \
  INSTANTIATE(Q, 192) \
  INSTANTIATE(Q, 256) \
  INSTANTIATE(Q, 320) \
  INSTANTIATE(Q, 384) \
  INSTANTIATE(Q, 448) \
  INSTANTIATE(Q, 512) \
  INSTANTIATE(Q, 576) \
  INSTANTIATE(Q, 1024) \
  INSTANTIATE(Q, 2048)
TQ_VALUES
#undef X

typedef void (*KernelFn)(uint8_t*, uint32_t* , uint8_t*, char*, size_t, uint32_t, uint8_t*, int, int, uint32_t*, uint32_t*, uint16_t*, uint16_t*, uint16_t*);
struct KernelSpec { int qSize; int sLen; KernelFn fn; };

typedef void (*KernelFn1)(uint8_t*, uint32_t* , uint8_t*, char*, size_t, uint32_t, uint8_t*, int, int, uint32_t*, uint32_t*);
struct KernelSpec1 { int qSize; int sLen; KernelFn1 fn; };

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};


std::unordered_map<std::pair<int, int>, KernelFn, pair_hash> kernelMap = {
    #define X(Q) \
        {{Q, 64}, gpuBLAZE_full<Q, 64, T_PER_PROB>}, \
        {{Q, 128}, gpuBLAZE_full<Q, 128, T_PER_PROB>}, \
        {{Q, 192}, gpuBLAZE_full<Q, 192, T_PER_PROB>}, \
        {{Q, 256}, gpuBLAZE_full<Q, 256, T_PER_PROB>}, \
        {{Q, 320}, gpuBLAZE_full<Q, 320, T_PER_PROB>}, \
        {{Q, 384}, gpuBLAZE_full<Q, 384, T_PER_PROB>}, \
        {{Q, 448}, gpuBLAZE_full<Q, 448, T_PER_PROB>}, \
        {{Q, 512}, gpuBLAZE_full<Q, 512, T_PER_PROB>}, \
        {{Q, 576}, gpuBLAZE_full<Q, 576, T_PER_PROB>}, \
        {{Q, 1024}, gpuBLAZE_full<Q, 1024, T_PER_PROB>},\
        {{Q, 2048}, gpuBLAZE_full<Q, 2048, T_PER_PROB>},  
    TQ_VALUES
    #undef X
};

std::unordered_map<std::pair<int, int>, KernelFn1, pair_hash> kernelMapEverything = {
    #define X(Q) \
        {{Q, 64}, gpuBLAZE_everything<Q, 64, T_PER_PROB>}, \
        {{Q, 128}, gpuBLAZE_everything<Q, 128, T_PER_PROB>}, \
        {{Q, 192}, gpuBLAZE_everything<Q, 192, T_PER_PROB>}, \
        {{Q, 256}, gpuBLAZE_everything<Q, 256, T_PER_PROB>}, \
        {{Q, 320}, gpuBLAZE_everything<Q, 320, T_PER_PROB>}, \
        {{Q, 384}, gpuBLAZE_everything<Q, 384, T_PER_PROB>}, \
        {{Q, 448}, gpuBLAZE_everything<Q, 448, T_PER_PROB>}, \
        {{Q, 512}, gpuBLAZE_everything<Q, 512, T_PER_PROB>}, \
        {{Q, 576}, gpuBLAZE_everything<Q, 576, T_PER_PROB>}, \
        {{Q, 1024}, gpuBLAZE_everything<Q, 1024, T_PER_PROB>},\
        {{Q, 2048}, gpuBLAZE_everything<Q, 2048, T_PER_PROB>}, 
    TQ_VALUES
    #undef X
};

#endif // End Header guard