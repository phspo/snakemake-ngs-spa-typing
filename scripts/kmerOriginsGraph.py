# Input: A kmer, A .json file encoding origns of kmers, A reference genome

import cairo


filename = 'test.svg'
WIDTH, HEIGHT = 512,512

surface = cairo.SVGSurface(filename, WIDTH, HEIGHT)
cr = cairo.Context(surface)

cr.set_source_rgb(0.3, 0.2, 0.5)  # Solid color
cr.rectangle(12,13,14,15)
cr.fill()

cr.save()
