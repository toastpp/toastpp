#!/usr/bin/env python
##########################################################################
# 
# Copyright 2009 VMware, Inc.
# Copyright 2008 Tungsten Graphics, Inc., Cedar Park, Texas.
# All Rights Reserved.
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sub license, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice (including the
# next paragraph) shall be included in all copies or substantial portions
# of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT.
# IN NO EVENT SHALL VMWARE AND/OR ITS SUPPLIERS BE LIABLE FOR
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
##########################################################################


from gallium import *
from base import *


def lods(*dims):
    size = max(dims)
    lods = 0
    while size:
        lods += 1
        size >>= 1
    return lods


def minify(dims, level = 1):
    return [max(dim>>level, 1) for dim in dims]


def tex_coords(texture, face, level, zslice):
    st = [ 
        [0.0, 0.0], 
        [1.0, 0.0], 
        [1.0, 1.0], 
        [0.0, 1.0],
    ] 
    
    if texture.target == PIPE_TEXTURE_2D:
        return [[s, t, 0.0] for s, t in st]
    elif texture.target == PIPE_TEXTURE_3D:
        depth = texture.get_depth(level)
        if depth > 1:
            r = float(zslice)/float(depth - 1)
        else:
            r = 0.0
        return [[s, t, r] for s, t in st]
    elif texture.target == PIPE_TEXTURE_CUBE:
        result = []
        for s, t in st:
            # See http://developer.nvidia.com/object/cube_map_ogl_tutorial.html
            sc = 2.0*s - 1.0
            tc = 2.0*t - 1.0
            if face == PIPE_TEX_FACE_POS_X:
                rx = 1.0
                ry = -tc
                rz = -sc
            if face == PIPE_TEX_FACE_NEG_X:
                rx = -1.0
                ry = -tc
                rz = sc
            if face == PIPE_TEX_FACE_POS_Y:
                rx = sc
                ry = 1.0
                rz = tc
            if face == PIPE_TEX_FACE_NEG_Y:
                rx = sc
                ry = -1.0
                rz = -tc
            if face == PIPE_TEX_FACE_POS_Z:
                rx = sc
                ry = -tc
                rz = 1.0
            if face == PIPE_TEX_FACE_NEG_Z:
                rx = -sc
                ry = -tc
                rz = -1.0
            result.append([rx, ry, rz])
        return result

def is_pot(n):
    return n & (n - 1) == 0
      
                
class TextureColorSampleTest(TestCase):
    
    tags = (
        'target',
        'format',
        'width',
        'height',
        'depth',
        'last_level',
        'face',
        'level',
        'zslice',
    )

    def test(self):
        dev = self.dev
        
        target = self.target
        format = self.format
        width = self.width
        height = self.height
        depth = self.depth
        last_level = self.last_level
        face = self.face
        level = self.level
        zslice = self.zslice
        
        tex_usage = PIPE_TEXTURE_USAGE_SAMPLER
        geom_flags = 0
        if width != height:
            geom_flags |= PIPE_TEXTURE_GEOM_NON_SQUARE
        if not is_pot(width) or not is_pot(height) or not is_pot(depth):
            geom_flags |= PIPE_TEXTURE_GEOM_NON_POWER_OF_TWO
        
        if not dev.is_format_supported(format, target, tex_usage, geom_flags):
            raise TestSkip
        
        ctx = self.dev.context_create()
    
        # disabled blending/masking
        blend = Blend()
        blend.rgb_src_factor = PIPE_BLENDFACTOR_ONE
        blend.alpha_src_factor = PIPE_BLENDFACTOR_ONE
        blend.rgb_dst_factor = PIPE_BLENDFACTOR_ZERO
        blend.alpha_dst_factor = PIPE_BLENDFACTOR_ZERO
        blend.colormask = PIPE_MASK_RGBA
        ctx.set_blend(blend)
    
        # no-op depth/stencil/alpha
        depth_stencil_alpha = DepthStencilAlpha()
        ctx.set_depth_stencil_alpha(depth_stencil_alpha)
    
        # rasterizer
        rasterizer = Rasterizer()
        rasterizer.front_winding = PIPE_WINDING_CW
        rasterizer.cull_mode = PIPE_WINDING_NONE
        rasterizer.bypass_vs_clip_and_viewport = 1
        ctx.set_rasterizer(rasterizer)
    
        # samplers
        sampler = Sampler()
        sampler.wrap_s = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.wrap_t = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.wrap_r = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.min_mip_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.min_img_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.mag_img_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.normalized_coords = 1
        sampler.min_lod = 0
        sampler.max_lod = PIPE_MAX_TEXTURE_LEVELS - 1
        ctx.set_sampler(0, sampler)
    
        #  texture 
        texture = dev.texture_create(
            target = target,
            format = format, 
            width = width, 
            height = height,
            depth = depth, 
            last_level = last_level,
            tex_usage = tex_usage,
        )
        
        expected_rgba = FloatArray(height*width*4) 
        texture.get_surface(
            face = face,
            level = level,
            zslice = zslice,
        ).sample_rgba(expected_rgba)
        
        ctx.set_sampler_texture(0, texture)

        #  framebuffer 
        cbuf_tex = dev.texture_create(
            PIPE_FORMAT_A8R8G8B8_UNORM, 
            width, 
            height,
            tex_usage = PIPE_TEXTURE_USAGE_RENDER_TARGET,
        )

        cbuf = cbuf_tex.get_surface()
        fb = Framebuffer()
        fb.width = width
        fb.height = height
        fb.nr_cbufs = 1
        fb.set_cbuf(0, cbuf)
        ctx.set_framebuffer(fb)
        rgba = FloatArray(4);
        rgba[0] = 0.5
        rgba[1] = 0.5
        rgba[2] = 0.5
        rgba[3] = 0.5
        ctx.clear(PIPE_CLEAR_COLOR, rgba, 0.0, 0)
        del fb
    
        # vertex shader
        vs = Shader('''
            VERT1.1
            DCL IN[0], POSITION, CONSTANT
            DCL IN[1], GENERIC, CONSTANT
            DCL OUT[0], POSITION, CONSTANT
            DCL OUT[1], GENERIC, CONSTANT
            0:MOV OUT[0], IN[0]
            1:MOV OUT[1], IN[1]
            2:END
        ''')
        #vs.dump()
        ctx.set_vertex_shader(vs)
    
        # fragment shader
        op = {
            PIPE_TEXTURE_1D: "1D", 
            PIPE_TEXTURE_2D: "2D", 
            PIPE_TEXTURE_3D: "3D", 
            PIPE_TEXTURE_CUBE: "CUBE",
        }[target]
        fs = Shader('''
            FRAG1.1
            DCL IN[0], GENERIC[0], LINEAR
            DCL OUT[0], COLOR, CONSTANT
            DCL SAMP[0], CONSTANT
            0:TEX OUT[0], IN[0], SAMP[0], %s
            1:END
        ''' % op)
        #fs.dump()
        ctx.set_fragment_shader(fs)

        nverts = 4
        nattrs = 2
        verts = FloatArray(nverts * nattrs * 4)
    
        x = 0
        y = 0
        w, h = minify((width, height), level)
    
        pos = [
            [x, y],
            [x+w, y],
            [x+w, y+h],
            [x, y+h],
        ]
    
        tex = tex_coords(texture, face, level, zslice)
    
        for i in range(0, 4):
            j = 8*i
            verts[j + 0] = pos[i][0] # x
            verts[j + 1] = pos[i][1] # y
            verts[j + 2] = 0.0 # z
            verts[j + 3] = 1.0 # w
            verts[j + 4] = tex[i][0] # s
            verts[j + 5] = tex[i][1] # r
            verts[j + 6] = tex[i][2] # q
            verts[j + 7] = 1.0
    
        ctx.draw_vertices(PIPE_PRIM_TRIANGLE_FAN,
                          nverts, 
                          nattrs, 
                          verts)
    
        ctx.flush()
    
        cbuf = cbuf_tex.get_surface()
        
        self.assert_rgba(cbuf, x, y, w, h, expected_rgba, 4.0/256, 0.85)
        

class TextureDepthSampleTest(TestCase):
    
    tags = (
        'target',
        'format',
        'width',
        'height',
        'depth',
        'last_level',
        'face',
        'level',
        'zslice',
    )

    def test(self):
        dev = self.dev
        
        target = self.target
        format = self.format
        width = self.width
        height = self.height
        depth = self.depth
        last_level = self.last_level
        face = self.face
        level = self.level
        zslice = self.zslice
        
        tex_usage = PIPE_TEXTURE_USAGE_SAMPLER
        geom_flags = 0
        if width != height:
            geom_flags |= PIPE_TEXTURE_GEOM_NON_SQUARE
        if not is_pot(width) or not is_pot(height) or not is_pot(depth):
            geom_flags |= PIPE_TEXTURE_GEOM_NON_POWER_OF_TWO
        
        if not dev.is_format_supported(format, target, tex_usage, geom_flags):
            raise TestSkip
        
        ctx = self.dev.context_create()
    
        # disabled blending/masking
        blend = Blend()
        blend.rgb_src_factor = PIPE_BLENDFACTOR_ONE
        blend.alpha_src_factor = PIPE_BLENDFACTOR_ONE
        blend.rgb_dst_factor = PIPE_BLENDFACTOR_ZERO
        blend.alpha_dst_factor = PIPE_BLENDFACTOR_ZERO
        blend.colormask = PIPE_MASK_RGBA
        ctx.set_blend(blend)
    
        # depth/stencil/alpha
        depth_stencil_alpha = DepthStencilAlpha()
        depth_stencil_alpha.depth.enabled = 1
        depth_stencil_alpha.depth.writemask = 1
        depth_stencil_alpha.depth.func = PIPE_FUNC_LESS
        ctx.set_depth_stencil_alpha(depth_stencil_alpha)
    
        # rasterizer
        rasterizer = Rasterizer()
        rasterizer.front_winding = PIPE_WINDING_CW
        rasterizer.cull_mode = PIPE_WINDING_NONE
        rasterizer.bypass_vs_clip_and_viewport = 1
        ctx.set_rasterizer(rasterizer)
    
        # samplers
        sampler = Sampler()
        sampler.wrap_s = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.wrap_t = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.wrap_r = PIPE_TEX_WRAP_CLAMP_TO_EDGE
        sampler.min_mip_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.min_img_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.mag_img_filter = PIPE_TEX_MIPFILTER_NEAREST
        sampler.normalized_coords = 1
        sampler.min_lod = 0
        sampler.max_lod = PIPE_MAX_TEXTURE_LEVELS - 1
        ctx.set_sampler(0, sampler)
    
        #  texture 
        texture = dev.texture_create(
            target = target,
            format = format, 
            width = width, 
            height = height,
            depth = depth, 
            last_level = last_level,
            tex_usage = tex_usage,
        )
        
        expected_rgba = FloatArray(height*width*4) 
        texture.get_surface(
            face = face,
            level = level,
            zslice = zslice,
        ).sample_rgba(expected_rgba)
        
        ctx.set_sampler_texture(0, texture)

        #  framebuffer 
        cbuf_tex = dev.texture_create(
            PIPE_FORMAT_A8R8G8B8_UNORM, 
            width, 
            height,
            tex_usage = PIPE_TEXTURE_USAGE_RENDER_TARGET,
        )

        zsbuf_tex = dev.texture_create(
            PIPE_FORMAT_Z24X8_UNORM, 
            width, 
            height,
            tex_usage = PIPE_TEXTURE_USAGE_RENDER_TARGET,
        )

        cbuf = cbuf_tex.get_surface()
        zsbuf = zsbuf_tex.get_surface()
        fb = Framebuffer()
        fb.width = width
        fb.height = height
        fb.nr_cbufs = 1
        fb.set_cbuf(0, cbuf)
        fb.set_zsbuf(zsbuf)
        ctx.set_framebuffer(fb)
        rgba = FloatArray(4);
        rgba[0] = 0.5
        rgba[1] = 0.5
        rgba[2] = 0.5
        rgba[3] = 0.5
        ctx.clear(PIPE_CLEAR_DEPTHSTENCIL, rgba, 1.0, 0)
        del fb
    
        # vertex shader
        vs = Shader('''
            VERT1.1
            DCL IN[0], POSITION, CONSTANT
            DCL IN[1], GENERIC, CONSTANT
            DCL OUT[0], POSITION, CONSTANT
            DCL OUT[1], GENERIC, CONSTANT
            0:MOV OUT[0], IN[0]
            1:MOV OUT[1], IN[1]
            2:END
        ''')
        #vs.dump()
        ctx.set_vertex_shader(vs)
    
        # fragment shader
        op = {
            PIPE_TEXTURE_1D: "1D", 
            PIPE_TEXTURE_2D: "2D", 
            PIPE_TEXTURE_3D: "3D", 
            PIPE_TEXTURE_CUBE: "CUBE",
        }[target]
        fs = Shader('''
            FRAG1.1
            DCL IN[0], GENERIC[0], LINEAR
            DCL SAMP[0], CONSTANT
            DCL OUT[0].z, POSITION
            0:TEX OUT[0].z, IN[0], SAMP[0], %s
            1:END
        ''' % op)
        #fs.dump()
        ctx.set_fragment_shader(fs)

        nverts = 4
        nattrs = 2
        verts = FloatArray(nverts * nattrs * 4)
    
        x = 0
        y = 0
        w, h = minify((width, height), level)
    
        pos = [
            [x, y],
            [x+w, y],
            [x+w, y+h],
            [x, y+h],
        ]
    
        tex = tex_coords(texture, face, level, zslice)
    
        for i in range(0, 4):
            j = 8*i
            verts[j + 0] = pos[i][0] # x
            verts[j + 1] = pos[i][1] # y
            verts[j + 2] = 0.0 # z
            verts[j + 3] = 1.0 # w
            verts[j + 4] = tex[i][0] # s
            verts[j + 5] = tex[i][1] # r
            verts[j + 6] = tex[i][2] # q
            verts[j + 7] = 1.0
    
        ctx.draw_vertices(PIPE_PRIM_TRIANGLE_FAN,
                          nverts, 
                          nattrs, 
                          verts)
    
        ctx.flush()
    
        zsbuf = zsbuf_tex.get_surface()
        
        self.assert_rgba(zsbuf, x, y, w, h, expected_rgba, 4.0/256, 0.85)
        



def main():
    dev = Device()
    suite = TestSuite()
    
    targets = [
        PIPE_TEXTURE_2D,
        PIPE_TEXTURE_CUBE,
        PIPE_TEXTURE_3D,
    ]
    
    color_formats = [
        PIPE_FORMAT_A8R8G8B8_UNORM,
        PIPE_FORMAT_X8R8G8B8_UNORM,
        #PIPE_FORMAT_A8R8G8B8_SRGB,
        PIPE_FORMAT_R5G6B5_UNORM,
        PIPE_FORMAT_A1R5G5B5_UNORM,
        PIPE_FORMAT_A4R4G4B4_UNORM,
        PIPE_FORMAT_A8_UNORM,
        PIPE_FORMAT_L8_UNORM,
        PIPE_FORMAT_YCBCR,
        PIPE_FORMAT_DXT1_RGB,
        #PIPE_FORMAT_DXT1_RGBA,
        #PIPE_FORMAT_DXT3_RGBA,
        #PIPE_FORMAT_DXT5_RGBA,
    ]
    
    depth_formats = [
        PIPE_FORMAT_Z32_UNORM,
        PIPE_FORMAT_Z24S8_UNORM,
        PIPE_FORMAT_Z24X8_UNORM,
        PIPE_FORMAT_Z16_UNORM,
    ]
    
    sizes = [64, 32, 16, 8, 4, 2, 1]
    #sizes = [1020, 508, 252, 62, 30, 14, 6, 3]
    #sizes = [64]
    #sizes = [63]
    
    faces = [
        PIPE_TEX_FACE_POS_X,
        PIPE_TEX_FACE_NEG_X,
        PIPE_TEX_FACE_POS_Y,
        PIPE_TEX_FACE_NEG_Y, 
        PIPE_TEX_FACE_POS_Z, 
        PIPE_TEX_FACE_NEG_Z,
    ]

    for format in color_formats:
        for target in targets:
            for size in sizes:
                if target == PIPE_TEXTURE_3D:
                    depth = size
                else:
                    depth = 1
                for face in faces:
                    if target != PIPE_TEXTURE_CUBE and face:
                        continue
                    levels = lods(size)
                    for last_level in range(levels):
                        for level in range(0, last_level + 1):
                            zslice = 0
                            while zslice < depth >> level:
                                test = TextureColorSampleTest(
                                    dev = dev,
                                    target = target,
                                    format = format, 
                                    width = size,
                                    height = size,
                                    depth = depth,
                                    last_level = last_level,
                                    face = face,
                                    level = level,
                                    zslice = zslice,
                                )
                                suite.add_test(test)
                                zslice = (zslice + 1)*2 - 1
    for format in depth_formats:
        target = PIPE_TEXTURE_2D
        depth = 1
        face = 0
        last_level = 0
        level = 0
        zslice = 0
        for size in sizes:
            test = TextureDepthSampleTest(
                dev = dev,
                target = target,
                format = format, 
                width = size,
                height = size,
                depth = depth,
                last_level = last_level,
                face = face,
                level = level,
                zslice = zslice,
            )
            suite.add_test(test)
    suite.run()


if __name__ == '__main__':
    main()
