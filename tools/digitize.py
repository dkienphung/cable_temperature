"""Digitize single-curve line plots from the rendered paper pages."""
import sys, json
import numpy as np
from PIL import Image

def load(path, box):
    im = np.asarray(Image.open(path).convert('L'), dtype=float)
    h, w = im.shape
    x0, y0, x1, y1 = [int(round(v)) for v in (box[0]*w, box[1]*h, box[2]*w, box[3]*h)]
    return im[y0:y1, x0:x1] < 110   # dark-pixel mask

def frame(mask):
    rows = mask.sum(1); cols = mask.sum(0)
    H, W = mask.shape
    r = np.where(rows > 0.5*W)[0]; c = np.where(cols > 0.5*H)[0]
    top, bot = r.min(), r.max(); left, right = c.min(), c.max()
    return top, bot, left, right

def ticks_x(mask, bot, left, right):
    # tick marks: short dark runs directly below the bottom frame line
    band = mask[bot+3:bot+12, left+3:right-2].sum(0)
    idx = np.where(band >= 4)[0] + left + 3
    groups = np.split(idx, np.where(np.diff(idx) > 2)[0]+1)
    return [g.mean() for g in groups if len(g)]

def ticks_y(mask, top, bot, left):
    band = mask[top+2:bot-2, left-12:left-3].sum(1)
    idx = np.where(band >= 4)[0] + top + 2
    groups = np.split(idx, np.where(np.diff(idx) > 2)[0]+1)
    return [g.mean() for g in groups if len(g)]

def curve(mask, top, bot, left, right):
    inner = mask[top+4:bot-3, left+4:right-3]
    xs, ys = [], []
    for j in range(inner.shape[1]):
        rr = np.where(inner[:, j])[0]
        if len(rr):
            xs.append(j + left + 4); ys.append(np.median(rr) + top + 4)
    return np.array(xs), np.array(ys)

def run(path, box, xvals, yvals):
    m = load(path, box)
    top, bot, left, right = frame(m)
    tx = ticks_x(m, bot, left, right); ty = ticks_y(m, top, bot, left)
    print('frame', top, bot, left, right, 'xticks', np.round(tx,1), 'yticks', np.round(ty,1), file=sys.stderr)
    # map using the frame edges = axis limits (xvals=(xmin,xmax), yvals=(ymin,ymax))
    px = np.polyfit([left, right], xvals, 1)
    py = np.polyfit([bot, top], yvals, 1)
    print('tick check x', np.round(np.polyval(px, tx), 2), 'y', np.round(np.polyval(py, ty), 4), file=sys.stderr)
    xs, ys = curve(m, top, bot, left, right)
    return np.polyval(px, xs), np.polyval(py, ys)

if __name__ == '__main__':
    cfg = json.loads(sys.argv[1])
    X, Y = run(cfg['path'], cfg['box'], np.array(cfg['xt'], float), np.array(cfg['yt'], float))
    grid = np.arange(cfg['g0'], cfg['g1']+1e-9, cfg.get('dg', 2.5))
    out = [(g, float(np.median(Y[np.abs(X-g) < 0.4]))) for g in grid if np.any(np.abs(X-g) < 0.4)]
    print(json.dumps(out))
