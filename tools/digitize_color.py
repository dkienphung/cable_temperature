"""Colour-based digitizer for the multi-curve result figures (Figs. 9, 10b, 11)."""
import sys, json
import numpy as np
from PIL import Image

def frame_and_ticks(rgb):
    dark = rgb.mean(2) < 90
    H, W = dark.shape
    r = np.where(dark.sum(1) > 0.5*W)[0]; c = np.where(dark.sum(0) > 0.5*H)[0]
    top, bot, left, right = r.min(), r.max(), c.min(), c.max()
    def groups(idx):
        g = np.split(idx, np.where(np.diff(idx) > 2)[0]+1)
        return [float(x.mean()) for x in g if len(x)]
    bx = dark[bot+3:bot+14, :].sum(0); tx = groups(np.where(bx >= 5)[0])
    by = dark[:, max(left-14,0):left-3].sum(1); ty = groups(np.where(by >= 5)[0])
    return (top, bot, left, right), tx, ty

def main():
    cfg = json.loads(sys.argv[1])
    rgb = np.asarray(Image.open(cfg['path']).convert('RGB'), dtype=float)
    (top, bot, left, right), tx, ty = frame_and_ticks(rgb)
    if cfg.get('probe'):
        print(json.dumps(dict(frame=[int(top), int(bot), int(left), int(right)], xt=tx, yt=ty)))
        return
    px = np.polyfit(cfg['xpix'], cfg['xval'], 1); py = np.polyfit(cfg['ypix'], cfg['yval'], 1)
    out = {}
    for name, spec in cfg['curves'].items():
        col = np.array(spec['rgb'], float); tol = spec.get('tol', 60)
        dist = np.sqrt(((rgb - col)**2).sum(2))
        mask = dist < tol
        mask[:top+4] = False; mask[bot-3:] = False; mask[:, :left+4] = False; mask[:, right-3:] = False
        for (x0, y0, x1, y1) in cfg.get('exclude', []):
            mask[y0:y1, x0:x1] = False
        xs, ys = [], []
        for j in range(mask.shape[1]):
            rr = np.where(mask[:, j])[0]
            if len(rr) >= spec.get('minpix', 1):
                xs.append(j); ys.append(np.median(rr))
        X = np.polyval(px, np.array(xs)); Y = np.polyval(py, np.array(ys))
        grid = np.arange(cfg['g0'], cfg['g1']+1e-9, cfg['dg'])
        hw = cfg['dg']/2
        pts = []
        for g in grid:
            sel = np.abs(X - g) <= hw
            if sel.sum() >= 2:
                pts.append([float(g), round(float(np.median(Y[sel])), 4)])
        out[name] = pts
    print(json.dumps(out))

main()
