#include "Rtypes.h" // for Double32_t

using namespace std;

void SetBitOn(UInt_t & var, UInt_t shift) {var = (var | (1 << shift));}
void SetBitOff(UInt_t & var, UInt_t shift) {var = (var & (~(1 << shift)));}
UInt_t GetBit(UInt_t var, UInt_t shift) {return ((var >> shift) & 1);}
