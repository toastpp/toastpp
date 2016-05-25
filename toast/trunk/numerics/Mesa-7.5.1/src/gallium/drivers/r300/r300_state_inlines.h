/*
 * Copyright 2009 Joakim Sindholt <opensource@zhasha.com>
 *                Corbin Simpson <MostAwesomeDude@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * on the rights to use, copy, modify, merge, publish, distribute, sub
 * license, and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHOR(S) AND/OR THEIR SUPPLIERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE. */

#ifndef R300_STATE_INLINES_H
#define R300_STATE_INLINES_H

#include "pipe/p_format.h"

#include "r300_reg.h"

/* Some maths. These should probably find their way to u_math, if needed. */

static INLINE int pack_float_16_6x(float f) {
    return ((int)(f * 6.0) & 0xffff);
}

/* Blend state. */

static INLINE uint32_t r300_translate_blend_function(int blend_func)
{
    switch (blend_func) {
        case PIPE_BLEND_ADD:
            return R300_COMB_FCN_ADD_CLAMP;
        case PIPE_BLEND_SUBTRACT:
            return R300_COMB_FCN_SUB_CLAMP;
        case PIPE_BLEND_REVERSE_SUBTRACT:
            return R300_COMB_FCN_RSUB_CLAMP;
        case PIPE_BLEND_MIN:
            return R300_COMB_FCN_MIN;
        case PIPE_BLEND_MAX:
            return R300_COMB_FCN_MAX;
        default:
            debug_printf("r300: Unknown blend function %d\n", blend_func);
            break;
    }
    return 0;
}

/* XXX we can also offer the D3D versions of some of these... */
static INLINE uint32_t r300_translate_blend_factor(int blend_fact)
{
    switch (blend_fact) {
        case PIPE_BLENDFACTOR_ONE:
            return R300_BLEND_GL_ONE;
        case PIPE_BLENDFACTOR_SRC_COLOR:
            return R300_BLEND_GL_SRC_COLOR;
        case PIPE_BLENDFACTOR_SRC_ALPHA:
            return R300_BLEND_GL_SRC_ALPHA;
        case PIPE_BLENDFACTOR_DST_ALPHA:
            return R300_BLEND_GL_DST_ALPHA;
        case PIPE_BLENDFACTOR_DST_COLOR:
            return R300_BLEND_GL_DST_COLOR;
        case PIPE_BLENDFACTOR_SRC_ALPHA_SATURATE:
            return R300_BLEND_GL_SRC_ALPHA_SATURATE;
        case PIPE_BLENDFACTOR_CONST_COLOR:
            return R300_BLEND_GL_CONST_COLOR;
        case PIPE_BLENDFACTOR_CONST_ALPHA:
            return R300_BLEND_GL_CONST_ALPHA;
        /* XXX WTF are these?
        case PIPE_BLENDFACTOR_SRC1_COLOR:
        case PIPE_BLENDFACTOR_SRC1_ALPHA: */
        case PIPE_BLENDFACTOR_ZERO:
            return R300_BLEND_GL_ZERO;
        case PIPE_BLENDFACTOR_INV_SRC_COLOR:
            return R300_BLEND_GL_ONE_MINUS_SRC_COLOR;
        case PIPE_BLENDFACTOR_INV_SRC_ALPHA:
            return R300_BLEND_GL_ONE_MINUS_SRC_ALPHA;
        case PIPE_BLENDFACTOR_INV_DST_ALPHA:
            return R300_BLEND_GL_ONE_MINUS_DST_ALPHA;
        case PIPE_BLENDFACTOR_INV_DST_COLOR:
            return R300_BLEND_GL_ONE_MINUS_DST_COLOR;
        case PIPE_BLENDFACTOR_INV_CONST_COLOR:
            return R300_BLEND_GL_ONE_MINUS_CONST_COLOR;
        case PIPE_BLENDFACTOR_INV_CONST_ALPHA:
            return R300_BLEND_GL_ONE_MINUS_CONST_ALPHA;
        /* XXX see above
        case PIPE_BLENDFACTOR_INV_SRC1_COLOR:
        case PIPE_BLENDFACTOR_INV_SRC1_ALPHA: */
        default:
            debug_printf("r300: Unknown blend factor %d\n", blend_fact);
            break;
    }
    return 0;
}

/* DSA state. */

static INLINE uint32_t r300_translate_depth_stencil_function(int zs_func)
{
    switch (zs_func) {
        case PIPE_FUNC_NEVER:
            return R300_ZS_NEVER;
        case PIPE_FUNC_LESS:
            return R300_ZS_LESS;
        case PIPE_FUNC_EQUAL:
            return R300_ZS_EQUAL;
        case PIPE_FUNC_LEQUAL:
            return R300_ZS_LEQUAL;
        case PIPE_FUNC_GREATER:
            return R300_ZS_GREATER;
        case PIPE_FUNC_NOTEQUAL:
            return R300_ZS_NOTEQUAL;
        case PIPE_FUNC_GEQUAL:
            return R300_ZS_GEQUAL;
        case PIPE_FUNC_ALWAYS:
            return R300_ZS_ALWAYS;
        default:
            debug_printf("r300: Unknown depth/stencil function %d\n",
                zs_func);
            break;
    }
    return 0;
}

static INLINE uint32_t r300_translate_stencil_op(int s_op)
{
    switch (s_op) {
        case PIPE_STENCIL_OP_KEEP:
            return R300_ZS_KEEP;
        case PIPE_STENCIL_OP_ZERO:
            return R300_ZS_ZERO;
        case PIPE_STENCIL_OP_REPLACE:
            return R300_ZS_REPLACE;
        case PIPE_STENCIL_OP_INCR:
            return R300_ZS_INCR;
        case PIPE_STENCIL_OP_DECR:
            return R300_ZS_DECR;
        case PIPE_STENCIL_OP_INCR_WRAP:
            return R300_ZS_INCR_WRAP;
        case PIPE_STENCIL_OP_DECR_WRAP:
            return R300_ZS_DECR_WRAP;
        case PIPE_STENCIL_OP_INVERT:
            return R300_ZS_INVERT;
        default:
            debug_printf("r300: Unknown stencil op %d", s_op);
            break;
    }
    return 0;
}

static INLINE uint32_t r300_translate_alpha_function(int alpha_func)
{
    switch (alpha_func) {
        case PIPE_FUNC_NEVER:
            return R300_FG_ALPHA_FUNC_NEVER;
        case PIPE_FUNC_LESS:
            return R300_FG_ALPHA_FUNC_LESS;
        case PIPE_FUNC_EQUAL:
            return R300_FG_ALPHA_FUNC_EQUAL;
        case PIPE_FUNC_LEQUAL:
            return R300_FG_ALPHA_FUNC_LE;
        case PIPE_FUNC_GREATER:
            return R300_FG_ALPHA_FUNC_GREATER;
        case PIPE_FUNC_NOTEQUAL:
            return R300_FG_ALPHA_FUNC_NOTEQUAL;
        case PIPE_FUNC_GEQUAL:
            return R300_FG_ALPHA_FUNC_GE;
        case PIPE_FUNC_ALWAYS:
            return R300_FG_ALPHA_FUNC_ALWAYS;
        default:
            debug_printf("r300: Unknown alpha function %d", alpha_func);
            break;
    }
    return 0;
}

/* Texture sampler state. */

static INLINE uint32_t r300_translate_wrap(int wrap)
{
    switch (wrap) {
        case PIPE_TEX_WRAP_REPEAT:
            return R300_TX_REPEAT;
        case PIPE_TEX_WRAP_CLAMP:
            return R300_TX_CLAMP;
        case PIPE_TEX_WRAP_CLAMP_TO_EDGE:
            return R300_TX_CLAMP_TO_EDGE;
        case PIPE_TEX_WRAP_CLAMP_TO_BORDER:
            return R300_TX_CLAMP_TO_BORDER;
        case PIPE_TEX_WRAP_MIRROR_REPEAT:
            return R300_TX_REPEAT | R300_TX_MIRRORED;
        case PIPE_TEX_WRAP_MIRROR_CLAMP:
            return R300_TX_CLAMP | R300_TX_MIRRORED;
        case PIPE_TEX_WRAP_MIRROR_CLAMP_TO_EDGE:
            return R300_TX_CLAMP_TO_EDGE | R300_TX_MIRRORED;
        case PIPE_TEX_WRAP_MIRROR_CLAMP_TO_BORDER:
            return R300_TX_CLAMP_TO_EDGE | R300_TX_MIRRORED;
        default:
            debug_printf("r300: Unknown texture wrap %d", wrap);
            return 0;
    }
}

static INLINE uint32_t r300_translate_tex_filters(int min, int mag, int mip)
{
    uint32_t retval = 0;
    switch (min) {
        case PIPE_TEX_FILTER_NEAREST:
            retval |= R300_TX_MIN_FILTER_NEAREST;
            break;
        case PIPE_TEX_FILTER_LINEAR:
            retval |= R300_TX_MIN_FILTER_LINEAR;
            break;
        case PIPE_TEX_FILTER_ANISO:
            retval |= R300_TX_MIN_FILTER_ANISO;
            break;
        default:
            debug_printf("r300: Unknown texture filter %d\n", min);
            break;
    }
    switch (mag) {
        case PIPE_TEX_FILTER_NEAREST:
            retval |= R300_TX_MAG_FILTER_NEAREST;
            break;
        case PIPE_TEX_FILTER_LINEAR:
            retval |= R300_TX_MAG_FILTER_LINEAR;
            break;
        case PIPE_TEX_FILTER_ANISO:
            retval |= R300_TX_MAG_FILTER_ANISO;
            break;
        default:
            debug_printf("r300: Unknown texture filter %d\n", mag);
            break;
    }
    switch (mip) {
        case PIPE_TEX_MIPFILTER_NONE:
            retval |= R300_TX_MIN_FILTER_MIP_NONE;
            break;
        case PIPE_TEX_MIPFILTER_NEAREST:
            retval |= R300_TX_MIN_FILTER_MIP_NEAREST;
            break;
        case PIPE_TEX_MIPFILTER_LINEAR:
            retval |= R300_TX_MIN_FILTER_MIP_LINEAR;
            break;
        default:
            debug_printf("r300: Unknown texture filter %d\n", mip);
            break;
    }

    return retval;
}

static INLINE uint32_t r300_anisotropy(float max_aniso)
{
    if (max_aniso >= 16.0f) {
        return R300_TX_MAX_ANISO_16_TO_1;
    } else if (max_aniso >= 8.0f) {
        return R300_TX_MAX_ANISO_8_TO_1;
    } else if (max_aniso >= 4.0f) {
        return R300_TX_MAX_ANISO_4_TO_1;
    } else if (max_aniso >= 2.0f) {
        return R300_TX_MAX_ANISO_2_TO_1;
    } else {
        return R300_TX_MAX_ANISO_1_TO_1;
    }
}

/* Buffer formats. */

static INLINE uint32_t r300_translate_colorformat(enum pipe_format format)
{
    switch (format) {
        /* 8-bit buffers */
        case PIPE_FORMAT_I8_UNORM:
            return R300_COLOR_FORMAT_I8;
        /* 16-bit buffers */
        case PIPE_FORMAT_R5G6B5_UNORM:
            return R300_COLOR_FORMAT_RGB565;
        case PIPE_FORMAT_A1R5G5B5_UNORM:
            return R300_COLOR_FORMAT_ARGB1555;
        case PIPE_FORMAT_A4R4G4B4_UNORM:
            return R300_COLOR_FORMAT_ARGB4444;
        /* 32-bit buffers */
        case PIPE_FORMAT_A8R8G8B8_UNORM:
        case PIPE_FORMAT_Z24S8_UNORM:
            return R300_COLOR_FORMAT_ARGB8888;
        /* XXX Not in pipe_format
        case PIPE_FORMAT_A32R32G32B32:
            return R300_COLOR_FORMAT_ARGB32323232;
        case PIPE_FORMAT_A16R16G16B16:
            return R300_COLOR_FORMAT_ARGB16161616;
        case PIPE_FORMAT_A10R10G10B10_UNORM:
            return R500_COLOR_FORMAT_ARGB10101010;
        case PIPE_FORMAT_A2R10G10B10_UNORM:
            return R500_COLOR_FORMAT_ARGB2101010;
        case PIPE_FORMAT_I10_UNORM:
            return R500_COLOR_FORMAT_I10; */
        default:
            debug_printf("r300: Implementation error: "
                "Got unsupported color format %s in %s\n",
                pf_name(format), __FUNCTION__);
            break;
    }
    return 0;
}

static INLINE uint32_t r300_translate_zsformat(enum pipe_format format)
{
    switch (format) {
        /* 16-bit depth, no stencil */
        case PIPE_FORMAT_Z16_UNORM:
            return R300_DEPTHFORMAT_16BIT_INT_Z;
        /* 24-bit depth, 8-bit stencil */
        case PIPE_FORMAT_Z24S8_UNORM:
            return R300_DEPTHFORMAT_24BIT_INT_Z_8BIT_STENCIL;
        default:
            debug_printf("r300: Implementation error: "
                "Got unsupported ZS format %s in %s\n",
                pf_name(format), __FUNCTION__);
            break;
    }
    return 0;
}

/* Translate pipe_format into US_OUT_FMT.
 * Note that formats are stored from C3 to C0. */
static INLINE uint32_t r300_translate_out_fmt(enum pipe_format format)
{
    switch (format) {
        case PIPE_FORMAT_A8R8G8B8_UNORM:
        case PIPE_FORMAT_Z24S8_UNORM:
            return R300_US_OUT_FMT_C4_8 |
                R300_C0_SEL_B | R300_C1_SEL_G |
                R300_C2_SEL_R | R300_C3_SEL_A;
        default:
            debug_printf("r300: Implementation error: "
                "Got unsupported output format %s in %s\n",
                pf_name(format), __FUNCTION__);
            return R300_US_OUT_FMT_UNUSED;
    }
    return 0;
}

/* Non-CSO state. (For now.) */

static INLINE uint32_t r300_translate_gb_pipes(int pipe_count)
{
    switch (pipe_count) {
        case 1:
            return R300_GB_TILE_PIPE_COUNT_RV300;
            break;
        case 2:
            return R300_GB_TILE_PIPE_COUNT_R300;
            break;
        case 3:
            return R300_GB_TILE_PIPE_COUNT_R420_3P;
            break;
        case 4:
            return R300_GB_TILE_PIPE_COUNT_R420;
            break;
    }
    return 0;
}

static INLINE uint32_t translate_vertex_data_type(int type) {
    switch (type) {
        case EMIT_1F:
        case EMIT_1F_PSIZE:
            return R300_DATA_TYPE_FLOAT_1;
            break;
        case EMIT_2F:
            return R300_DATA_TYPE_FLOAT_2;
            break;
        case EMIT_3F:
            return R300_DATA_TYPE_FLOAT_3;
            break;
        case EMIT_4F:
            return R300_DATA_TYPE_FLOAT_4;
            break;
        case EMIT_4UB:
            return R300_DATA_TYPE_BYTE;
            break;
        default:
            debug_printf("r300: Implementation error: "
                    "Bad vertex data type!\n");
            assert(0);
            break;
    }

    return 0;
}

#endif /* R300_STATE_INLINES_H */
