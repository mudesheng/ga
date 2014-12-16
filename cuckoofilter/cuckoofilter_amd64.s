
TEXT ·CompareAndSwapUint16(SB),4,$0-13
	MOVQ	addr+0(FP), BP
	MOVW	old+8(FP), AX
	MOVW	new+10(FP), CX
	LOCK
	CMPXCHGL	CX, 0(BP)
	SETEQ	swapped+12(FP)
	RET
