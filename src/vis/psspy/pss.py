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
delay = 0
color_head = 'gray'
color_i = 'black'
size_i = 100
marker_i = '+'
lw_i = 1.5
zorder_i = 100
grid_buckets = np.sqrt(64)
mk_image_path = ''
mk_image_ppi = 300
silent = False

# figure size / [mm]
w = 210.0/2 - 35
h = w

def read_next_pointset(delimiter = ' '):
  pts = []
  for ln in ism:
    # line ln constraints \n (end-of-line feed); need to remove it
    ln = ln.rstrip()
    # read until end-of-set (aka #eos; by custom convention) obtained
    if (ln == eos):
      break
    # skip comments
    if (len(ln) == 0 or ln[0] == '#'):
      continue
    # interpret point coordinates as array
    pt = ln.split(delimiter)
    # append point pt to pointset pts
    pts.append(pt)
  return np.array(pts).astype(np.float)

def draw_boundary(ax):
  ax.hlines([0,1], 0, 1, ls='-', color='#000', lw=0.75, zorder=3)
  ax.vlines([0,1], 0, 1, ls='-', color='#000', lw=0.75, zorder=3)

def draw_grid(ax, grid_buckets):
  l =  np.arange(grid_buckets)[1:] / grid_buckets
  ax.hlines(l, 0, 1, ls='-', color='#aaa', lw=.5, zorder=2)
  ax.vlines(l, 0, 1, ls='-', color='#aaa', lw=.5, zorder=2)   

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
  draw_grid(ax, grid_buckets)
  ax.set_xticks(np.linspace(0,1,5))
  ax.set_yticks(np.linspace(0,1,5))
  ax.tick_params(labelleft=False, labelbottom=False, left=True, right=True, bottom=True, top=True)

  # viewport: domain, aspect ratio
  dlim=1.0/32.0
  ax.set_xlim(-dlim, 1+dlim)
  ax.set_ylim(-dlim, 1+dlim)
  ax.set_aspect(1)

  return fig, ax

def draw_pointset(ax, pts, color='black', size=100, marker='+', linewidth=1.5, zorder=100):
  # need to split coordinates
  p0 = pts[:,0]
  p1 = pts[:,1]

  # draw points
  # - need to remember shapes for future removal (aka frame update)
  pts_shapes = ax.scatter(p0, p1, c=color, s=size, marker=marker, linewidth=linewidth, zorder=zorder)

  return pts_shapes

def show_figure(fig):
  plt.ion()
  plt.show()

def draw_figure(fig):
  # plt.draw()
  # need to flush events (and not use plt.pause(0.001) which would steel window focus)
  fig.canvas.flush_events()

def drawing_style(i):
  if (i == 0):
    color = color_head
  else:
    color = color_i
  return color

def visualise():
  # internal
  seq_i = 0
  pts_head = []
  pts_i = []
  pts_shape_head = []
  pts_shape_i = []

  # initialize figure
  fig, ax = init_figure(w, h)

  # show figure with non-blocking (allows future figure updates)
  show_figure(fig)
  draw_figure(fig)

  # iterate through input of pointset sequence
  while True:
    # retrieve point se
    pts_i = read_next_pointset()

    # stop reading on empty point set
    if (len(pts_i) == 0):
      break

    # option: delay between sequence elements (visual inspection by humans)
    # (who (human) is faster than a computer?)
    if (delay > 0):
      time.sleep(delay)

    # remove previous pointset, except the head pointset (first one)
    if (seq_i > 1):
      pts_shape_i.remove()

    # draw pointset
    # - the first pointset gets color_head, subsequent one color_i  
    color = drawing_style(seq_i)
    pts_shape_i = draw_pointset(ax, pts_i, color=color, size=size_i, marker=marker_i, linewidth=lw_i, zorder=zorder_i)
    
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
        print(f'sequence={seq_i}   mkout={p}', end='\r')
      else:
        print(f'sequence={seq_i}', end='\r')

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

parser.add_argument('--i', default='', help='A path to a point set (sequence) to be visualised')
parser.add_argument('--delay', type=float, default=delay, help='Number of seconds to delay between frame updates')
parser.add_argument('--image-path', default='', help='A path to images generated during this sequence, containing {i}. Example: "seq-{i}.png"')
parser.add_argument('--image-ppi', type=float, default='300', help='Resolution of images to be generated, in ppi unit. Example: 300.')
parser.add_argument('--silent', action='store_true', help='Hide logging/debuggin output.')

args = parser.parse_args()
delay = args.delay
mk_image_path = args.image_path
mk_image_ppi = args.image_ppi
silent = args.silent

if (args.i != ''):
  ism = open(args.i, 'rt')

# start visualisation of the point set sequence
visualise()