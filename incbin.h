#pragma once
#define INCBIN(name, filename) \
    extern const unsigned char _binary_##name##_start[]; \
    extern const unsigned char _binary_##name##_end[]; \
    __asm__( \
        ".section .rodata\n" \
        ".global _binary_" #name "_start\n" \
        "_binary_" #name "_start:\n" \
        ".incbin \"" filename "\"\n" \
        ".global _binary_" #name "_end\n" \
        "_binary_" #name "_end:\n" \
    );
