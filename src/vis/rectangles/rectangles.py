#!/bin/python3

import sys
import time
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# from datetime import datetime

# constants
# matplotlib calculates in inch
# 1 mm := 0.0393701 inch
mm = 0.0393701
eos = "#eos"

# parameters
ism = sys.stdin
ism_pts = None
delay = 0
colour_edge = (0,0,0,0.05)
colour_fill = (0,0,0,0.05)
colour_pts = (0,0,0,1)
fill = False
size_i = 100
marker_i = '+'
lw_i = 1.5
zorder_i = 100
grid_buckets = np.sqrt(64)
mk_image_path = ''
mk_image_ppi = 300
silent = False
domain = [0,0,1,1]
gridlines = [7, 7]

# figure size / [mm]
w = 210.0/2 - 35
h = w

def read_next_pointset(stream, delimiter = ' '):
  pts = []
  eof = True
  for ln in stream:
    # line ln constraints \n (end-of-line feed); need to remove it
    ln = ln.rstrip()
    # read until end-of-set (aka #eos; by custom convention) obtained
    if (ln == eos):
      # definitely a valid set (the empty set)
      eof = False
      break
    # skip comments
    if (len(ln) == 0 or ln[0] == '#'):
      continue
    # receiving numerical input
    eof = False    
    # interpret point coordinates as array
    pt = ln.split(delimiter)
    # append point pt to pointset pts
    pts.append(pt)
  return (np.array(pts).astype(np.float), eof)

def draw_boundary(ax):
  ax.hlines([domain[1], domain[3]], domain[0], domain[2], ls='-', color='#000', lw=0.75, zorder=3)
  ax.vlines([domain[0], domain[2]], domain[1], domain[3], ls='-', color='#000', lw=0.75, zorder=3)

def draw_grid(ax):
  l = np.linspace(domain[1], domain[3], gridlines[1])
  ax.hlines(l, domain[0], domain[2], ls='-', color='#aaa', lw=.5, zorder=2)
  l = np.linspace(domain[0], domain[2], gridlines[0])
  ax.vlines(l, domain[1], domain[3], ls='-', color='#aaa', lw=.5, zorder=2)

def draw_axes_ticks(ax):
  ax.set_xticks(np.linspace(domain[0],domain[2], gridlines[0]))
  ax.set_yticks(np.linspace(domain[1],domain[3], gridlines[1]))
  ax.tick_params(labelleft=False, labelbottom=False, left=True, right=True, bottom=True, top=True)

def limit_axes(ax):
  dlim=1.0/32.0
  ax.set_xlim(domain[0]-dlim, domain[2]+dlim)
  ax.set_ylim(domain[1]-dlim, domain[3]+dlim)
  ax.set_aspect(1)

def init_figure(w, h):
  w_ = w
  h_ = h
  t = 0.95
  l = 0.05
  b = 0.05
  r = 1-b

  # figure environment
  fig, ax = plt.subplots(1,1, figsize=(w_*mm, h_*mm))
  plt.subplots_adjust(left=l, bottom=b, right=r, top=t)

  # # boundary, grid, axis ticks
  draw_boundary(ax)
  draw_grid(ax)
  draw_axes_ticks(ax)
  limit_axes(ax)  

  return fig, ax

def draw_rectangles(ax, rects, pts, size=100, marker='+', linewidth=1.5, zorder=100):
  # skip nil drawing
  if (rects.shape[0] == 0):
    return []

  # need to split coordinates
  # storage format: rectangle: (left bottom right top)
  # or: (lower-point upper-point)
  rd0v0 = rects[:,0]
  rd1v0 = rects[:,1]
  rd0v1 = rects[:,2]
  rd1v1 = rects[:,3]

  # draw points
  # - need to remember shapes for future removal (aka frame update)
  r = []
  for i in np.arange(rd0v0.shape[0]):
    re = matplotlib.patches.Rectangle((rd0v0[i],rd1v0[i]), rd0v1[i]-rd0v0[i], rd1v1[i]-rd1v0[i], edgecolor=colour_edge, facecolor=colour_fill, zorder=zorder, fill=fill)
    r.append(ax.add_patch(re))

  if (pts.shape[0] > 0):
    r.append(ax.scatter(pts[:,0], pts[:,1], color=colour_pts, s=size, marker=marker, linewidth=linewidth, zorder=zorder+1))

  return r

def clear_shapes(shapes):
  for i in np.arange(len(shapes)):
    shapes[i].remove()

def show_figure(fig):
  plt.ion()
  plt.show()

def draw_figure(fig):
  # plt.draw()
  # need to flush events (and not use plt.pause(0.001) which would steel window focus)
  fig.canvas.flush_events()

def visualise():
  # internal
  seq_i = 0
  pts_head = []
  pts_i = []
  pts_shape_head = []
  pts_shape_i = []
  ptspt_i = np.array([])
  ism_eof = True
  ism_pts_eof = True

  # initialize figure
  fig, ax = init_figure(w, h)

  # show figure with non-blocking (allows future figure updates)
  show_figure(fig)
  draw_figure(fig)

  # iterate through input of pointset sequence
  while True:
    # retrieve point se
    pts_i, ism_eof = read_next_pointset(ism)

    if (ism_pts != None):
      ptspt_i, ism_pts_eof = read_next_pointset(ism_pts)

    # stop reading on empty point set
    # if (len(pts_i) == 0):
    if (ism_eof):
      break

    # option: delay between sequence elements (visual inspection by humans)
    # (who (human) is faster than a computer?)
    if (delay > 0):
      time.sleep(delay)

    # remove previous pointset, except the head pointset (first one)
    if (seq_i > 1):
      clear_shapes(pts_shape_i)

    # draw pointset
    # - the first pointset gets colour_edge, subsequent one colour_i  
    pts_shape_i = draw_rectangles(ax, pts_i, ptspt_i, size=size_i, marker=marker_i, linewidth=lw_i, zorder=zorder_i)

    # flush drawing commands (for 55 points, might need about 30ms)
    # print(datetime.now().time())
    draw_figure(fig)
    
    # need to remember initial pointset
    if (seq_i == 0):
      pts_head = pts_i
      pts_shape_head = pts_shape_i

    # remember iteration number
    seq_i = seq_i + 1

    # option: store figure
    if (mk_image_path != ''):
      p = mk_image_path.format(i=seq_i)
      plt.savefig(p, dpi=mk_image_ppi)

    # show console output
    if (not silent):
      if (mk_image_path != ''):
        print(f'sequence={seq_i}   mkout={p}     ', end='\r')
      else:
        print(f'sequence={seq_i}   count={pts_i.shape[0]}      ', end='\r')

  if (not silent):
    print()
    print(f'sequence size: {seq_i}')

  # keep figure open (requires non-interactive mode aka ioff)
  # and plt.show enters blocking mode
  plt.ioff()
  plt.show()

  # clear everything
  fig.clear()
  plt.close(fig)


# handle program arguments
parser = argparse.ArgumentParser(description='Visualise a point set sequence')

parser.add_argument('--domain', type=float, nargs=4, default=domain, help='Problem domain in d=2 dimensions, formatted as [min d=0, min d=1, max d=0, max d=1]. Default: [0,0,1,1].')
parser.add_argument('--gridlines', type=int, nargs=2, default=gridlines, help='Number of gridlines within the d=2 dimensional problem domain, included the domain boundary itself, formatted as [num d=0, num d=1]. Default: [7,7].')

parser.add_argument('--i', help='A path to a sequence of rectangles to be visualised')
parser.add_argument('--pts', help='A path to a point set (sequence) to be visualised')
parser.add_argument('--delay', type=float, default=delay, help='Number of seconds to delay between frame updates')
parser.add_argument('--image-path', default='', help='A path to images generated during this sequence, containing {i}. Example: "seq-{i}.png"')
parser.add_argument('--image-ppi', type=float, default='300', help='Resolution of images to be generated, in ppi unit. Example: 300.')
parser.add_argument('--silent', action='store_true', help='Hide logging/debuggin output.')

parser.add_argument('--colour-alpha-edge', type=float, help='Alpha value of drawing each rectangle')
parser.add_argument('--colour-rgba-edge', type=float, nargs=4, help='rgba value of drawing each rectangle')
parser.add_argument('--colour-edge', help='Colour value (C0,...,C9; hex colour, colour name) of drawing each rectangle')

parser.add_argument('--colour-alpha-fill', type=float, help='Alpha value of drawing each rectangle')
parser.add_argument('--colour-rgba-fill', type=float, nargs=4, help='rgba value of drawing each rectangle')
parser.add_argument('--colour-fill', help='Colour value (C0,...,C9; hex colour, colour name) of drawing each rectangle')
parser.add_argument('--fill', type=int, help='Enable to fill rectangles')

parser.add_argument('--colour-alpha-pts', type=float, help='Alpha value of drawing points')
parser.add_argument('--colour-rgba-pts', type=float, nargs=4, help='rgba value of drawing points')
parser.add_argument('--colour-pts', help='Colour value (C0,...,C9; hex colour, colour name) of drawing points')

args = parser.parse_args()
delay = args.delay
mk_image_path = args.image_path
mk_image_ppi = args.image_ppi
domain = args.domain
gridlines = args.gridlines

if (args.colour_rgba_edge != None):
  colour_edge = tuple(args.colour_rgba_edge)
if (args.colour_edge != None):
  colour_edge = matplotlib.colors.to_rgba(args.colour_edge)
if (args.colour_alpha_edge != None):
  colour_edge = (colour_edge[0], colour_edge[1], colour_edge[2], args.colour_alpha_edge)

if (args.colour_rgba_fill != None):
  colour_fill = tuple(args.colour_rgba_fill)
if (args.colour_fill != None):
  colour_fill = matplotlib.colors.to_rgba(args.colour_fill)
if (args.colour_alpha_fill != None):
  colour_fill = (colour_fill[0], colour_fill[1], colour_fill[2], args.colour_alpha_fill)

if (args.colour_rgba_pts != None):
  colour_pts = tuple(args.colour_rgba_pts)
if (args.colour_pts != None):
  colour_pts = matplotlib.colors.to_rgba(args.colour_pts)
if (args.colour_alpha_pts != None):
  colour_pts = (colour_pts[0], colour_pts[1], colour_pts[2], args.colour_alpha_pts)

if (args.fill != None):
  fill = args.fill != 0

silent = args.silent

if (args.i != None):
  ism = open(args.i, 'rt')

if (args.pts != None):
  ism_pts = open(args.pts, 'rt')

# start visualisation of the point set sequence
visualise()