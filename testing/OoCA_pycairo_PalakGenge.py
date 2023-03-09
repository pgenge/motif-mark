#!/usr/bin/env python

import cairo
import math

#I tried to make a faux rothko-esque "painting"

#setting width and height parameters
width, height = 400, 400

#make drawing surface and save as pdf
surface = cairo.PDFSurface("OoCA_PalakGenge.pdf",width, height)
context = cairo.Context(surface)

#drawing line
context.set_line_width(3)
context.move_to(25,130)        #(x,y)
context.line_to(375,130)
context.stroke()

#drawing and filling rectangles
context.rectangle(25,25,350,100)
context.set_source_rgb(0.27, 0.37, 0.19)
context.fill()
context.rectangle(25,150,350,200)
context.set_source_rgb(0.85, 0.0, 0.3)
context.fill()

#writing out the picture
surface.finish()

