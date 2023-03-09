#!/usr/bin/env python


import argparse
from bioinfo import oneline_fasta, iupac_notation_dict
import math
import gzip
import re 
import cairo
import itertools

#Argparse for input and output files
def get_args():
    #optional arguments
    parser = argparse.ArgumentParser(description = 'File Inputs for Motif Marking of Fasta')
    parser.add_argument('-f', '--fasta_file', help = 'Specify input .fasta file')
    parser.add_argument('-m', '--motifs_file', help = 'Specify input motifs text file')
    # parser.add_argument('-o', '--output_file', help = 'Specify output image (png) file as same name as fasta file name') #naming output based on input so not needed
    return parser.parse_args()
args = get_args()


#open the files, get fasta file name for final file naming and load iupac notation
fasta = oneline_fasta(args.fasta_file) #make seq in fasta file one line
fa = open('onelinefasta.fa', 'r') #open the file
fasta_name = args.fasta_file.split(".")[0] #get the file name to name output file
motifs = open(args.motifs_file, 'r') #open known motifs text file
iupac = iupac_notation_dict #this is a dictionary of the IUPAC notation for unknown/ambiguous base calls

#fix ambiguous sequences in motifs file and add them to list to be referenced later
motif_list = [] #make empty list to add motifs in motif text file
unambig_motifs = [] #make empty list to add ambiguous motifs that have been replaced with the IUPAC notation

for line in motifs:
  line = line.strip('\n').upper() #take each line in motifs file and make it upper case
  motif_seq = str(line) #make the line a string and each string is a motif sequence
  #print(motif_seq)
  motif_list.append(motif_seq) #add all motifs including ambiguous to motif list
  #print(motif_list)
#take each motif in motif list and replace the ambiguous motifs to have IUPAC notation so that it can be searched for in the fasta seq
for motif_seq in motif_list: 
  for key, value in iupac.items(): #get the iupac motif single base and the bases of what it could be if ambiguous
      motif_seq = motif_seq.replace(key, value) 
  unambig_motifs.append(motif_seq) #add the replaced ambigous motifs to a list
  
#print(unambig_motifs)

class gene:
   """
   This is the gene class. A gene is a proteing coding DNA sequence that can be found in a fasta file
   """
   def __init__(self, gene_start, gene_stop, gene_pos):
      self.start = gene_start #gene start position
      self.stop = gene_stop #gene stop (can be the length of gene)
      self.pos = gene_pos #y position on surface

      #methods: draw using pycairo, set line thickness, colors, positioning and write out
   def draw_gene(self):
    """draws gene on cairo surface"""
    context.set_line_width(1)
    context.set_source_rgba(0,0,0) 
    context.move_to(self.start + margin, self.pos)
    context.line_to(self.stop + margin, self.pos)
    context.stroke()

class exon:
   """
   This is an exon. An exon is a feature of a DNA sequence that contains the information for coding a protein.
   """
   def __init__(self, exon_start, exon_stop, exon_pos):
      self.start = exon_start #exon start pos
      self.stop = exon_stop #exon stop position (can be length of exon)
      self.pos = exon_pos #y position on surface

    #methods: draw using pycairo, set line thickness, colors, positioning and write out
   def draw_exon(self):
      """draws exon on gene"""
      context.set_line_width(10)
      context.set_source_rgb(0,0,0)
      context.move_to(self.start + margin, self.pos)
      context.line_to(self.stop + margin, self.pos)
      context.stroke()

#set colors for motifs to be used when calling motif method draw_motif
yellow = [0.8, 0.7, 0, 0.8]
green = [0.4, 0.7, 0.2, 0.8]
blue = [0.4, 0.6, 1, 0.82]
pink = [0.9, 0.5, 0.6, 0.82]
orange = [1, 0.5, 0.1, 0.82] #maximum of 5 motifs
#context.set_source_rgb(0,0,0) black
colors_list = [yellow, green, blue, pink, orange] #make a list of colors to iterate through later
#print(colors_list)

class motif:
   """
   This is a Motif which is a region of DNA that has a short 
   pattern of nucleotides. Motifs serve specific biological functions such as 
   being the binding site for a regulatory protein etc.
   """
   def __init__(self, motif_start, motif_stop, motif_pos):
      self.start = motif_start
      self.stop = motif_stop
      self.pos = motif_pos

   #methods: draw using pycairo, set line thickness, colors, positioning and write out   
   def draw_motif(self, colors_list):
      """draws motif seq on gene"""
      context.set_line_width(9)
      context.set_source_rgba(colors_list[0], colors_list[1], colors_list[2], colors_list[3])
      context.move_to(self.start + margin, self.pos)
      context.line_to(self.stop + margin, self.pos)
      context.stroke()

#Making surface and Context for cairo output
#where pycairo should start for drawing on page
margin = 50 #set pycairo margins for where to start drawing and to align
#set surface and context for pycairo to draw on
width = 1300 + margin #hard coded width :(
height = 800 #hardcoded heigh of surface :(
pngname = (f'{fasta_name}.png') #make the output file name
surface = cairo.PDFSurface("Figure_1.pdf", width, height) #make a pdf surface first then use this to write to png
context = cairo.Context(surface) #give cairo the surface to draw on and call it context
context.set_source_rgb(1,1,1) #paint background white
context.paint()
context.set_source_rgb(0,0,0) #set font color to black
context.select_font_face("Calibri", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
context.set_font_size(15)
context.move_to(15, 25)
context.show_text(f'Motif Marking on Genes from Fasta file: {fasta_name}.fa')
context.save()

#create a key
keyaligner = 1300-300 #making an aligner so that text lines up on the key
keyspaces = 50  #use for y space between lines
context.move_to(keyaligner, margin)
context.show_text('Key:')

#draw gene key
context.set_font_size(13)
context.move_to(keyaligner, margin + keyspaces)
context.show_text('Gene/Intron')
gene_line = gene(keyaligner + 60, keyaligner + 90, margin + 45)
gene_line.draw_gene()

#draw exon key
context.move_to(keyaligner, margin + keyspaces + 30)
context.show_text('Exon/Cassette')
exon_line = exon(keyaligner + 60, keyaligner + 90, margin + keyspaces + 25)
exon_line.draw_exon()

#draw motif key and color code the motif key's text to corresponding motif on gene
motifkey_num = 0
for motif_seq, color in zip(motif_list, colors_list):
  motif_color = colors_list[motifkey_num]
  context.move_to(keyaligner, margin + keyspaces + 65)
  #print(motif)
  context.select_font_face("Calibri", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
  context.set_source_rgba(motif_color[0], motif_color[1], motif_color[2], motif_color[3])
  context.show_text(motif_seq)
  keyspaces += 30
  motifkey_num += 1

# write the code to actually draw things out using classes and methods
seq_id_ypos = 60 #setting y position for gene name
seq_ypos = 80 #setting y position for gene/intro

for lines in fa:
  #line = fa.readline().strip('\n') #don't have to do this/will mess up the line calling
  lines = lines.strip('\n') #strip new line/get each line
  if lines == '':
     break #finish when there are no lines left
  else:
    seq = lines
    #this chunk is to write out header
    if seq.startswith('>'): #grab header/gene name
        context.set_source_rgb(0,0,0)
        context.move_to(margin, seq_id_ypos)
        context.set_font_size(12.5)
        gene_name = seq.strip('>') #remove > from header
        #print(gene_name)
        context.show_text(f'{gene_name}') #write out gene name given the position
        seq_id_ypos += 70 #create space for gene line and this will become next header position
    else:
        #draw gene
        gene_seq = gene(0, len(seq), seq_ypos) #create gene object and pass start, stop based on length of gene, and y position)
        gene_seq.draw_gene() #draw gene using class method
        #find the exon
        exon_seq = re.finditer('[A-Z]+', seq) #find all exons in gene using regex
        for start_stop in exon_seq: 
            exon_start_stop = start_stop.span() #spanning the sequence of the exon will give us the start and stop positions of this exon
        exon_draw = exon(exon_start_stop[0], exon_start_stop[1], seq_ypos) #make exon object given the start, stop from span and the ypos which should be same as the gene we're currently on
        exon_draw.draw_exon()
        motif_num = 0    #have to initialize the count of motifs so that it can be iterated through
        #find and draw motifs
        for entry,color in zip(unambig_motifs, colors_list): #this accesses the elements inside both lists
            motif_start_stop = re.finditer(str(entry), lines) #find all motifs in each gene line of the fasta
            for start_stop in motif_start_stop: 
              motif_color = colors_list[motif_num] #this iterates through colors list using the motif number and assigns the a corresponding color in the list to the motif that it's on
              motif_start_stop = start_stop.span() #this finds the start and stop positions of the motif
              #print(motif_color)
              motif_seq = motif(motif_start_stop[0], motif_start_stop[1], seq_ypos) #make motif class and pass start, stop, and yposition (same as gene)
              motif_seq.draw_motif(motif_color) #draw motif using motif methods
            motif_num += 1 #increment motif so that it can iterate through all motifs and assign colors accordingly
        seq_ypos += 70 #increment yposition of motif drawing so that it aligns with the genes each time


surface.write_to_png (pngname) #make png file using the file name of the input fasta
surface.finish() #have to do this to make pycairo write out
fa.close() #close fasta
motifs.close() #close motif text file