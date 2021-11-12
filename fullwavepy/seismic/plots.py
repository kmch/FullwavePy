"""
Various plots combining objects defined 
in different modules.
"""
from autologging import logged
import matplotlib.pyplot as plt

from arrau.a3d import Arr3d
from fullwavepy.plot.misc import plot_square
from plotea.mpl2d import figure

@logged
class PlotExp:
  """
  A plotting mix-in for PROTEUS(Experiment).
  """
  def plot_acq(self, extent=[[-6e4,6e4],[-12e3,3e4]], acq=True, bt=True, \
    box=True, zoom=True, fig=None, markersize=7, cmap=None):
    if self.all_not_read:
      self.read_all()
    if fig is None:
      figure(20,7)
    if bt:
      self._plot_bathy_topo(extent, cmap)
    if box:
      self._plot_box(extent)
    if acq:
      self._plot_stations(extent, markersize)
    if zoom:
      self._zoom(extent)
  def _plot_bathy_topo(self, extent=None, cmap=None):
    if extent is None:
      bt = self.bt
    else:
      bt = self.bt.extract(extent)
    if cmap is None:
      bt.plot()
    else:
      bt.plot(cmap=cmap)
  def _plot_box(self, extent, **kwargs):
    if not extent is None:
      [[x1, x2], [y1, y2]] = extent
      plot_square(x1, x2, y1, y2, **kwargs)
  def _plot_stations(self, extent=None, markersize=7):
    """
    Plot acquisition with mpl.
    """
    srcs, recs = self.srcs, self.recs
    sio, who, lan = [self.pool[key] for key in ['sio', 'who', 'lan']]
    kw_recs = dict(linestyle='', \
      markersize=markersize, 
            markeredgecolor='k', markeredgewidth=1, markerfacecolor='w')
    plt.plot(srcs.sx, srcs.sy, linestyle='', marker='.', label='shot',
        markersize=1, markerfacecolor='k', markeredgecolor='k')
    plt.plot(sio.gx, sio.gy, marker='o', label='SIO', **kw_recs)
    plt.plot(who.gx, who.gy, marker='^', label='WHOI', **kw_recs)
    plt.plot(lan.gx, lan.gy, marker='s', label='land', **kw_recs)
    shift = 3e2
    for ID, x, y in zip(recs.id, recs.gx, recs.gy):
      xytext = (x+shift, y+shift)
      plt.annotate(s=str(ID), xy=(x,y), xytext=xytext, clip_on=True,
            bbox={'facecolor': 'w', 'edgecolor': 'k', 'alpha': .5, },
            # arrowprops={'arrowstyle': '->', 'color': 'k'}
            )
    plt.gca().set_aspect('equal')
  def _zoom(self, extent):
    if not extent is None:
      [[x1, x2], [y1, y2]] = extent
      plt.xlim(x1,x2)
      plt.ylim(y1,y2)
# OLD Plots for ch08_Kolumbo_volcano/Results.ipynb
def plot_bathy(bt, srcs, recs, sio, who, lan, x1, x2, y1, y2):
  bt.plot(mode='shade', cmap='cmo.deep', label='metres b.s.l.', alpha=.2)
  bt.plot(mode='cr', colors='w')
  plot_acq_geom(srcs, recs, sio, who, lan)
  plt.xlim(x1,x2)
  plt.ylim(y1,y2)
def plot_fit(fit, it, sids):
  assert it > 0
  for sid in sids:
    y = fit[sid][:it]
    x = np.arange(1, len(y)+1)
    
    plt.plot(x, y, marker='o', c='r', alpha=.5)
  plt.xlabel('iteration')
  plt.ylabel('total misfit (%)')
  plt.xlim(-1,161)
  plt.ylim(0,100)
def plot_mod(proj, box, bt, ref, it, axis=1, value=9e3, vmin=-1.5e3):
  a = proj.o.vp.it[it].read()
  vp = Arr3d(a, extent=a.extent).extract(box.extent)
  vp = Arr3d(vp.arr - ref.arr, extent=box.extent)
  clip = 0 # m/s
  vp.arr[vp.arr > clip] = 0

  kws = dict(value=value, unit='m', axis=axis)
  vp.slice(**kws).plot(cmap='cet_fire_r', # 'hot_r'
             label='vp, m/s',
             vmin=vmin, vmax=0,
            )
  vp.slice(**kws).plot(mode='cr', colors='w')
  if axis != 2:
    bt.slice(**kws).plot(color='Grey')
  # vp.slices.list[0].plot_slice_lines()
    plt.gca().invert_yaxis()   
def plot_phase(proj, it, sid, freq, x1, x2, y1, y2):
  self = proj.o.dc.it[it][sid]
  self._get_phase(freq=freq)
  ph_type = 'dif'
  plt.scatter(self.head.sx, self.head.sy, vmin=-np.pi, vmax=+np.pi, cmap='cet_bkr',
           c=self.head['phase %s (%s Hz)' % (ph_type, freq)])
  plt.scatter(self.head.gx[0], self.head.gy[0], s=20**2, 
          marker='*', c='w', edgecolors='k')
  plt.xlim(x1, x2)
  plt.ylim(y1, y2)
  plt.gca().set_aspect('equal')
