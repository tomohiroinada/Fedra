#ifndef _SYSAL_TRACK_TYPE_2_
#define _SYSAL_TRACK_TYPE_2_

#include "tvectors.h"

#pragma pack(push)
#pragma pack(4)

typedef struct
{
	unsigned short Area;
	float X;
	float Y;
	float Z;
	} Grain;

typedef struct
{
	int Field;
	unsigned Grains;
	unsigned AreaSum;
	Grain *pGrains;
	BYTE *pCorrection;
	TVector Intercept;
	TVector Slope;
	float Sigma;
	float FirstZ;
	float LastZ;
	boolean Valid;
	} Track2;

#pragma pack(pop)

#endif