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
graph_layout = False
axes_yscale = 'linear'
axes_xscale = 'linear'
lw_i = 1.5
zorder_i = 100
delay = 0
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

# def draw_boundary(ax):
#   ax.hlines([domain[1], domain[3]], domain[0], domain[2], ls='-', color='#000', lw=0.75, zorder=3)
#   ax.vlines([domain[0], domain[2]], domain[1], domain[3], ls='-', color='#000', lw=0.75, zorder=3)

# def draw_grid(ax):
#   l = np.linspace(domain[1], domain[3], gridlines[1])
#   ax.hlines(l, domain[0], domain[2], ls='-', color='#aaa', lw=.5, zorder=2)
#   l = np.linspace(domain[0], domain[2], gridlines[0])
#   ax.vlines(l, domain[1], domain[3], ls='-', color='#aaa', lw=.5, zorder=2)

# def draw_axes_ticks(ax):
#   ax.set_xticks(np.linspace(domain[0],domain[2], gridlines[0]))
#   ax.set_yticks(np.linspace(domain[1],domain[3], gridlines[1]))
#   ax.tick_params(labelleft=False, labelbottom=False, left=True, right=True, bottom=True, top=True)

# def limit_axes(ax):
#   dlim=1.0/32.0
#   ax.set_xlim(domain[0]-dlim, domain[2]+dlim)
#   ax.set_ylim(domain[1]-dlim, domain[3]+dlim)
#   ax.set_aspect(1)

def init_figure(w, h):
  w_ = w
  h_ = h
  t = 0.95
  l = 0.05
  b = 0.05
  r = 1-b

  # figure environment
  fig, ax = plt.subplots(1,1, figsize=(w_*mm, h_*mm))
  # plt.subplots_adjust(left=l, bottom=b, right=r, top=t)

  # boundary, grid, axis ticks
  # draw_boundary(ax)
  # draw_grid(ax)
  # draw_axes_ticks(ax)
  # limit_axes(ax)  

  ax.set_xscale(axes_xscale)
  ax.set_yscale(axes_yscale)

  return fig, ax

def draw_graphs(ax, graphs, linewidth=1.5, zorder=100):
  # skip nil drawing
  if (graphs.shape[0] == 0):
    return []

  # draw points
  # - need to remember shapes for future removal (aka frame update)
  r = []
  arg = []
  d = np.arange(graphs.shape[1])

  if (graph_layout):
    d = d[1:]
    arg = graphs[:,0]
  else:
    arg = np.arange(graphs.shape[0])

  for i in d:
    val = graphs[:,i]
    r.append(ax.plot(arg, val))

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
  ism_eof = True

  # initialize figure
  fig, ax = init_figure(w, h)

  # show figure with non-blocking (allows future figure updates)
  show_figure(fig)
  draw_figure(fig)

  # iterate through input of pointset sequence
  while True:
    # retrieve point se
    pts_i, ism_eof = read_next_pointset(ism)

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
    pts_shape_i = draw_graphs(ax, pts_i, linewidth=lw_i, zorder=zorder_i)

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
parser.add_argument('--gridlines', type=int, nargs=2, default=gridlines, help='Number of gridlines within the d=2 dimensional problem domain, including the domain boundary itself, formatted as [num d=0, num d=1]. Default: [7,7].')

parser.add_argument('--i', help='A path to a sequence of graphs to be visualised')
parser.add_argument('--graph-layout', action='store_true', help='Whether the graphs contain arguments or not')

parser.add_argument('--delay', type=float, default=delay, help='Number of seconds to delay between frame updates')
parser.add_argument('--image-path', default='', help='A path to images generated during this sequence, containing {i}. Example: "seq-{i}.png"')
parser.add_argument('--image-ppi', type=float, default='300', help='Resolution of images to be generated, in ppi unit. Example: 300.')
parser.add_argument('--silent', action='store_true', help='Hide logging/debuggin output.')

parser.add_argument('--xscale', default=axes_xscale, help='The x-axis scale. Could be: {"linear", "log", "symlog", "logit", ...}')
parser.add_argument('--yscale', default=axes_yscale, help='The y-axis scale. Could be: {"linear", "log", "symlog", "logit", ...}')

args = parser.parse_args()
delay = args.delay
mk_image_path = args.image_path
mk_image_ppi = args.image_ppi
domain = args.domain
gridlines = args.gridlines

axes_xscale = args.xscale
axes_yscale = args.yscale

if (args.graph_layout != None):
  graph_layout = args.graph_layout

silent = args.silent

if (args.i != None):
  ism = open(args.i, 'rt')

# start visualisation of the point set sequence
visualise()