"""
"""
from autologging import logged, traced
import matplotlib.pyplot as plt
import numpy as np

from fullwavepy.ndat.arrays import Arr3d # a-b
from fullwavepy.plot.generic import * 
from fullwavepy.project.generic.au import ProjBox
from fullwavepy.project.generic.au import ProjGeometry
from fullwavepy.project.types.basic import ProjSyn, ProjInv

# PROJECT PARAMS
@logged
class Discret(ProjGeometry):
    def __init__(self, dt, ns, dx, dims):
        self.dt = dt
        self.ns = ns
        self.dx = dx
        self.dict = dict(dt=dt, ns=ns, dx=dx, dims=dims)
    def check(self, f_max, v_min=1500, v_max=8000):
        from fwipy.solve.api import check_discret
        self.__log.info('\n\nLOW KERNEL')
        check_discret(self.dx, self.dt, f_max, v_min, v_max, kernel='low')
        self.__log.info('\n\nHIGH KERNEL')
        check_discret(self.dx, self.dt, f_max, v_min, v_max, kernel='high')
class Box(ProjBox):
    def __init__(self, x1, x2, y1, y2, z1, z2):
        self.box = [x1, x2, y1, y2, z1, z2]
        self.extent = [[x1, x2], [y1, y2], [z1, z2]]
        self.x1, self.x2, self.y1, self.y2, self.z1, self.z2 = self.box
        self._init_slices()
    def _init_slices(self):
        self.sl = {}
    def plot(self, figsize=None, label=None, aspect='equal'):
        from fullwavepy.plot.misc import plot_box
        if figsize is not None:
            figure(*figsize)
        kwargs = dict(label=label)
        plt.subplot(131)
        plot_box(self.box[0], self.box[1], 
                 self.box[2], self.box[3], **kwargs)
        plt.gca().set_aspect(aspect)

        plt.subplot(132)
        plot_box(self.box[0], self.box[1], 
                 self.box[4], self.box[5], **kwargs)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(aspect)

        plt.subplot(133)
        plot_box(self.box[2], self.box[3], 
                 self.box[4], self.box[5], **kwargs)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(aspect)
    def plot_zslice(self, figsize=None, label=None, aspect='equal', **kwargs):
        from fullwavepy.plot.misc import plot_box
        if figsize is not None:
            figure(*figsize)
        kwargs = dict(label=label)
        plot_box(self.box[0], self.box[1], 
                 self.box[2], self.box[3], **kwargs)
        plt.gca().set_aspect(aspect)
def box2dims(box, dx):
    x1, x2, y1, y2, z1, z2 = box
    assert x2 >= x1
    assert y2 >= y1
    assert z2 >= z1
    nx1 = int((x2 - x1) / dx) + 1 
    nx2 = int((y2 - y1) / dx) + 1  
    nx3 = int((z2 - z1) / dx) + 1     
    return nx1, nx2, nx3

# WORKFLOWS
def prep_inp_syn(p, tvp, rsg, sp, rnf, cat=0):
    p.i.tvp.create(tvp, cat=cat)
    p.i.rsg.create(rsg, cat=cat)
    p.i.sp.create(**dict(sp), cat=cat)
    p.i.sp.run(cat=cat)
    p.i.rnf.create(**rnf, cat=cat)
    if cat:
        p.i.rnf.cat()
def syn(pname, box, dis, acq, svp, rsg, bnd, fs, fw, exe, eph, ts, **kwargs):
    """
    """
    plot = kwargs.get('plot', False)
    # 1. INIT
    path = path_eph if eph else './'
    p = ProjSyn(pname, box=box.box, **dis.dict, path=path, exe=exe, env={'SLAVES_WAVEFIELDSVTR': fw}, cat=0)
    # 2. PREP
    ebnd = bnd + 10
    rnf = dict(btop=bnd, etop=ebnd, b_abs=bnd, e_abs=ebnd)
    rnf = dict(rnf, btop=0, etop=0) if fs else rnf
    prep_inp_syn(p, svp, rsg, acq.sp, rnf)
    # 3. RUN
    no = 0
    p.i.bash.no[no].prep(cat=0)
    p.o.rm(ls=0)
    p.i.bash.no[no].run(truncate=1000)
    # 4. PLOT
    if plot:
        syn_plot_out(p, ts)
    return p
def shot_snap(p, sid, it, ts, bw=False, **kwargs):
    verbose = kwargs.get('verbose', True)
    srcs = p.i.s.read()
    
    if 'x' in kwargs and 'y' in kwargs and 'z' in kwargs:
        x = kwargs.get('x')
        y = kwargs.get('y')
        z = kwargs.get('z')
    else:
        # round to nearest integer
        x, y, z = [[int(round(j)) for j in i] for i in srcs.li if i.ID == sid][0]
        x = x + p.elef - 1 # -1 to convert grid node to array index
        y = y + p.efro - 1
        z = z + p.etop - 1
    if verbose:
        print('x,y,z', x,y,z)
    a = kwargs.get('clip', 10)

    hw = kwargs.get('hw', 45)
    if verbose:
        print('hw: ', hw)    
    if bw:
        f = p.o.bw.it[it][sid][ts]
    else:
        f = p.o.fw.it[it][sid][ts]
    
    kwgs = dict(vmin=-a, vmax=a, noextent=True)
    kwgs['overwrite_mmp'] = kwargs.get('overwrite_mmp', True)
    kwgs['overwrite'] = kwargs.get('overwrite', True)    
    
    if verbose:
        print('kwgs', kwgs)
    # xlim = (80,130)
    # ylim = (50,0)
    
    plt.figure(figsize=(18,6))
    plt.suptitle(f.fname)
    plt.subplot(131)
    f.read().plot(x=x, **kwgs)
    if hw is not None:
        plt.xlim(y-hw, y+hw)
        plt.ylim(z+hw, z-hw)
    plt.gca().set_aspect('equal')
    plt.subplot(132)
    f.read().plot(y=y,**kwgs)
    if hw is not None:
        plt.xlim(x-hw, x+hw)
        plt.ylim(z+hw, z-hw)
    plt.gca().set_aspect('equal')    
    plt.subplot(133)
    f.read().plot(z=z, **kwgs)
    if hw is not None:
        plt.xlim(x-hw, x+hw)
        plt.ylim(y-hw, y+hw)  
    plt.gca().set_aspect('equal')   
def syn_plot_out(p, ts):
    figure(14,5)
    plt.subplot(221)
    _ = plt.plot(p.i.rsg.read(overwrite_mmp=1, overwrite=1)[0,0]) 
    plt.subplot(222)
    _ = plt.plot(p.o.syn.read(overwrite_mmp=1, overwrite=1)[0,0])
    plt.subplot(223)
    _ = plt.imshow(p.i.tvp.read(overwrite_mmp=1, overwrite=1)[:,0,:].T)
    plt.subplot(224)
    _ = plt.imshow(p.o.fw.it[1][1][ts].read(overwrite_mmp=1, overwrite=1)[:,0,:].T)

# MODELS AND THEIR INTERFACES
def a_minus_b(a, b, clip=500, **kwargs):
    kwargs['cmap'] = kwargs.get('cmap', 'RdBu') 
    Arr3d(a.read() - b.read()).plot(vmin=-clip, vmax=clip, **kwargs)
def extent2absorb(extent, dx, etop, e_abs):
    """
    Extend `extent` to absorbing layers.
    """    
    if etop <= 2:
        etop = 2
    [[x1,x2],[y1,y2],[z1,z2]] = extent
    elef, erig, efro, ebac, ebot = [e_abs] * 5
    x1 = x1 - elef * dx
    y1 = y1 - efro * dx
    z1 = z1 - etop * dx
    x2 = x2 + erig * dx
    y2 = y2 + ebac * dx
    z2 = z2 + ebot * dx
    return [[x1,x2],[y1,y2],[z1,z2]]
def extract_vp_and_fs(exp_svp, exp_bt, box, dx, plot=True):
    assert dx == exp_svp.dx[0]
    assert len(set(exp_svp.dx == 1)) # dx=dy=dz
    vp = exp_svp.copy().carve(box.box)
    node_fs = -box.z1 / dx # z1 was assumed to be final one
    bt = exp_bt.carve(box.box)
    bt.extent = bt.extent[:-1] # collapse from 3d to 2d
    fs = bt.extract_freesurf(add=node_fs, dx=dx) # note the sign
    if plot:
        plot_bt_and_model(bt, vp, k_fs=node_fs-1)
    return vp, fs
def check_compatibility_vp_and_fs(vp, fs, k_fs, vp_min=0):
    def format_plot(val=None):
        plt.gca().set_aspect('equal')
        plt.gca().invert_yaxis()
        try:
            plt.colorbar()
        except:
            print('UFuncTypeError: Cannot plot the colorbar when const (%s)'% val)
    kws = dict(vmin=-1, vmax=1)
    for k in np.arange(k_fs+1,k_fs-9,-1):
        plt.figure(figsize=(18,5))

        plt.subplot(131)
        plt.title('Vp > %s at depth node %s' % (str(vp_min), str(k)))
        a = (vp > vp_min).astype(int)[...,k]
        plt.imshow(a.T, **kws)
        format_plot()
        
        plt.subplot(132)
        plt.title('fs < %s' % str(k))
        b = (fs < k).astype(int)[...,0]
        # kfs = k + 2
        # plt.title('fs < %s' % str(kfs))
        # b = (fs < kfs - 1).astype(int)[...,0] # NOTE test 21-05-20
        plt.imshow(b.T, **kws)
        format_plot()

        plt.subplot(133)
        plt.title('Difference (middle - left)')
        plt.imshow((b - a).T, **kws)
        format_plot(val=a[0])  
def make_vp_compatible_with_fs(vp, fs, vel_air):
    assert vp.shape[:-1] == fs.shape[:-1]
    nvp = np.copy(vp)
    nx, ny, nz = vp.shape
    n = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz-1,-1,-1):
                if fs[i,j] < k and nvp[i,j,k] == 0:
                    new = nvp[i,j,k+1]
                    assert new > 0
                    print('Replacing nvp[%s,%s,%s]=%s with %s' % (i,j,k,nvp[i,j,k],new))
                    nvp[i,j,k] = new
                    n += 1
    print('Updated %s values' % n)   


    # print('min velocity: ', np.min(nvp))
    # nvp_min = np.min(nvp[nvp>0])
    # print('min non-zero velocity: ', nvp_min)
    # dv = 20
    # print('Setting air velocity %s m/s smaller than current min' % dv)
    # nnvp_min = nvp_min - dv
    nnvp_min = vel_air
    nnvp = np.clip(nvp, nnvp_min, None)
    print('new min velocity: ', nnvp_min)

    plt.figure(figsize=(15,5))
    _ = plt.hist(np.ravel(nnvp), bins=50, range=(0,4000))
    plt.xlabel('vp [m/s]')
    plt.ylabel('no. of grid cells')
    plt.grid()

    return nnvp  
def compare_vp_and_fs(vp, fs, k_fs):
    plt.figure(figsize=(14,7))
    plt.subplot(121)
    ax = vp.plot(z=k_fs)
    aspeqt(ax)
    plt.subplot(122)
    ax = fs.plot(z=0, cmap='nipy_spectral')
    # levels 7 doesn't work <= 
    # plt.contour(fs_kameni[...,0].T, levels=[7]) #, extent=fs_kameni.extent[:-1].flatten(), levels=[7], colors='k')
    plt.gca().set_aspect('equal')
def plot_bt_and_model(bt, vp, k_fs):
    def figure():
        plt.figure(figsize=(20,10))
    figure()
    k = dict(aspect='e')
    plt.subplot(121)
    vp.plot(z=k_fs, cmap='magma', **k)
    plt.subplot(122)
    bt.plot(z=0, **k)
    figure()
    plt.subplot(121)
    bt.extract_freesurf(add=0, dx=vp.dx[0]).plot(z=0, cmap=[], center_cmap=1, **k)
    plt.subplot(122)
    bt.extract_seabed(dx=vp.dx[0]).plot(z=0, cmap='magma', center_cmap=0, **k)
    plt.title('extracted seabed')
    print(vp.shape, bt.shape)
def qc_mod3d(mod, z_zoom):
    figure(16,5)
    plt.subplot(121)
    plt.imshow(mod[0].T)
    plt.xlim(-2,2)
    plt.ylim(2,-2)
    plt.subplot(122)
    plt.imshow(mod[0].T)
    plt.ylim(z_zoom+2, z_zoom-2) 
def plot_4zslices(a, z, xlim, ylim):
    def fmt():
        plt.xlim(xlim)
        plt.ylim(ylim)
        aspeqt(plt.gca())
    assert len(z) == 4
    figure(14,12)
    plt.subplot(221)
    a.plot(z=z[0])
    fmt()
    plt.subplot(222)
    a.plot(z=z[1])
    fmt()
    plt.subplot(223)
    a.plot(z=z[2])
    fmt()
    plt.subplot(224)
    a.plot(z=z[3])
    fmt()
def nb_scroller():
  fig = figure(5,20)
  tracker = (vp != 0).scrollall(fig, cmap='viridis')
  return fig.canvas.mpl_connect('scroll_event', tracker.onscroll)

# DATA QC AND PROCESSING
def qc_datafile(datafile, ep, cmap1='Greys', cmap2='hot', \
    **kwargs):
    txlim = kwargs.get('txlim', None)
    tylim = kwargs.get('tylim', None)
    fxlim = kwargs.get('fxlim', None)
    fylim = kwargs.get('fylim', None)
    kwargs['win'] = dict(ep=[ep])
    datafile.read(**kwargs)
    figure(16,8)
    plt.suptitle(datafile.name + ', line ' + str(ep))    
    plt.subplot(121)
    datafile.array.plot(cmap=cmap1, center_cmap=1, **kwargs)
    plt.xlim(txlim)
    plt.ylim(tylim)
    plt.xlabel('trace no.')
    plt.ylabel('sample')
    plt.subplot(122)
    datafile.array.plot(cmap=cmap2, center_cmap=0, spect='ampl', dt=datafile.dt, **kwargs)
    plt.xlim(fxlim)
    plt.ylim(fylim)
    plt.xlabel('trace no.')
    plt.ylabel('frequency [Hz]')
    plt.gca().set_aspect('auto') 
def qc_filt(p, psyn, sid=None, ep=None, overwrite=False, overwrite_mmp=False):
    sids = [s.ID for s in p.i.s.read().li]
    eps = sorted(p.i.obs.read_header()['ep'].unique())
    if sid is None:
        sid = sids[0]    
    if ep is None:
        ep = eps[0]    
    kwargs = dict(win=dict(tracf=[sid], ep=[ep]), norm='max', 
                  overwrite=overwrite, overwrite_mmp=overwrite_mmp)
    plt.figure(figsize=(25,14))
    plt.suptitle('\nQC of data processing. From left to right: raw, filtered, '+
                 'final (filtered & muted), synthetic.',fontsize=25)
    xlabel = 'trace index [-]'
    tylabel = 'time sample [-]'
    fylabel = 'frequency [Hz]'
    
    plt.subplot(2,4,1)
    p.i.obs.raw.plot(**kwargs, cbar=0)
    plt.xlabel(xlabel)
    plt.ylabel(tylabel)
    plt.subplot(2,4,2)
    p.i.obs.fil.plot(**kwargs, cbar=0)
    plt.xlabel(xlabel)
    #     plt.ylabel(tylabel)    
    plt.subplot(2,4,3)
    p.i.obs.plot(**kwargs, cbar=0)
    plt.xlabel(xlabel)
    #     plt.ylabel(tylabel)
    plt.subplot(2,4,4)
    psyn.o.syn.plot(**kwargs)    
    plt.xlabel(xlabel)
    #     plt.ylabel(tylabel)    
    
    #     kwargs = dict(kwargs, norm=None)
    fkwargs = dict(spect='ampl', dt=p.dt, cmap='hot', center_cmap=0)
    ylim = (10,0)
    plt.subplot(2,4,5)
    p.i.obs.raw.plot(**kwargs, **fkwargs, cbar=0)
    plt.ylim(ylim)
    plt.xlabel(xlabel)
    plt.ylabel(fylabel)     
    plt.subplot(2,4,6)
    p.i.obs.fil.plot(**kwargs, **fkwargs, cbar=0)
    plt.ylim(ylim)
    plt.xlabel(xlabel)
    #     plt.ylabel(fylabel)       
    plt.subplot(2,4,7)
    p.i.obs.plot(**kwargs, **fkwargs, cbar=0)    
    plt.ylim(ylim)
    plt.xlabel(xlabel)
    #     plt.ylabel(fylabel)       
    plt.subplot(2,4,8)
    psyn.o.syn.plot(**kwargs, **fkwargs)    
    plt.ylim(ylim)     
    #     plt.ylabel(fylabel)        
    plt.xlabel(xlabel)
def plot_data_and_geom(datafile, md, sid, ep, **kwargs):
    from fullwavepy.plot.plt2d import colorbar
    save = kwargs.get('save', False)
    df = md[md.tracf==sid]
    attr = 'ep'
    val = ep
    figx = kwargs.get('figx', 16)
    figy = kwargs.get('figx', 14)
    txlim = kwargs.get('txlim', None)
    tylim = kwargs.get('tylim', None)
    fxlim = kwargs.get('fxlim', None)
    fylim = kwargs.get('fylim', None)
    kwargs['win'] = dict(tracf=[sid], ep=[ep])
    kwargs['norm'] = 'max'
    datafile.read(**kwargs)

    from matplotlib.gridspec import GridSpec
    gs = GridSpec(2,2, width_ratios=[2,1], height_ratios=[1,2])
    fig = figure(figx,figy)
    fig.add_subplot(gs[0,0])
    plot_acq_geom(df, attr, val, sid, ep)
    fig.add_subplot(gs[1,:])
    datafile.array.plot(cmap='seismic', center_cmap=1, **kwargs)
    plt.xlim(txlim)
    plt.ylim(tylim)
    plt.xlabel('trace no.')
    plt.ylabel('sample')
    fig.add_subplot(gs[0,1])
    try:
      dt = datafile.dt
    except AttributeError:
      dt = md['dt'][0]
      if dt > 1000: # microsec
        dt = dt / 1e6 # sec
    datafile.array.plot(cmap='hot', center_cmap=0, spect='ampl', dt=dt, **kwargs)
    plt.xlim(fxlim)
    plt.ylim(fylim)
    plt.xlabel('trace no.')
    plt.ylabel('frequency [Hz]')
    if save:
        plt.savefig('dataqc_sid%s_ep%s.png' % (sid, ep))
        plt.close()
def set_ticks(datafile, **kwargs):
    #     datafile.read_header()
    #     decim = 10 #kwargs['decim']
    #     hw = 'offset' #kwargs['hw']
    #     locs = np.arange(len(datafile.head))[::decim]
    #     labels = datafile.head[hw][::decim]
    #     rotation = lambda decim : np.clip(90 - 10 * (decim - 1), 0, 90)
    #     _ = plt.xticks(locs, labels, rotation=rotation(decim))
    #     _ = plt.xlabel(hw)
    #     _ = plt.ylabel('sample')
    pass # it's not that easy because we have to window the header too!
def plot_acq_geom(df, attr, val, sid, ep):
    from fullwavepy.plot.plt2d import colorbar
    # plt.title('Station %s, SEGY keyword: %s; highlighted value: %s' % (sid, attr, val))
    plt.title('Station %s, line: %s' % (sid, ep))
    sc = plt.scatter(df['sx'], df['sy'], c=df[attr])
    colorbar(sc, plt.gca())
    ndf = df[df[attr] == val]
    plt.scatter(ndf['sx'], ndf['sy'], c='red')
    plt.scatter(df['gx'], df['gy'], c='magenta') #marker=dict(color=color, size=size))
    aspeqt(plt.gca())
def plot_station_pools():
    # great plot with different pools (WHOI, SIO, land)
    # markersize = 50 
    # plt.figure(figsize=(10,7))
    # # topo.plot(center_cmap=1, extent=topo.extent[:-2])
    # kwargs = dict(marker='^', s=markersize)
    # # rr = r[r.pool=='SIO']
    # # plt.scatter(rr.gx, rr.gy, c='b', **kwargs)
    # # rr = r[r.pool=='WHOI']
    # # plt.scatter(rr.gx, rr.gy, c='r', **kwargs)
    # # rr = r[r.pool=='land']
    # # plt.scatter(rr.gx, rr.gy, c='y', **kwargs)
    # # shift = 1e2
    # # for ID, x, y in zip(r.tracf, r.gx, r.gy):
    # #     plt.annotate(s=str(ID)[1:], xy=(x+shift, y+shift), clip_on=True)
    # plt.scatter(s.sx, s.sy, c=s.ep, cmap='hsv', s=10, vmin=40, vmax=43)
    # plt.colorbar()
    # # p01.pbox.plot(c='b')
    # # plt.xlim(xlim)
    # # plt.ylim(ylim)
    # # plt.gca().set_aspect('equal')
    # plt.grid('--', c='Grey')
    pass
def get_phase(dc, freq):
  try:
    phi = dc.read_header(overwrite=0)['phase dif (%s Hz)' % freq]
  except KeyError:
    dc._get_phase(freq)
    phi = dc.read_header(overwrite=0)['phase dif (%s Hz)' % freq]
  return phi
def phase_hists(p, sid, freq, it2, it1=1, bins=40):
    dc1 = p.o.dc.it[it1][sid]
    dc2 = p.o.dc.it[it2][sid]
    
    ph1 = get_phase(dc1, freq)
    ph2 = get_phase(dc2, freq)
    
    figure(10,6)
    ph1.hist(color='r', label='iteraton %s' % it1, bins=bins)
    ph2.hist(color='k', label='iteraton %s' % it2, bins=bins, alpha=0.7)
    _ = plt.xlabel('phase difference [rad]')
    _ = plt.ylabel('counts')    
    plt.legend()
def find_shoot_dir(lid, md):
    line = md[md['ep']==lid]
    f0, f1 = line['fldr'][:2]
    x0, x1 = line['sx'][:2]
    if f0 > f1:
        raise ValueError 
    direction = 'left_to_right' if x1 > x0 else 'right_to_left'
    return direction
def find_marginal_shot(lid, md, shoot_dir):
    line = md[md['ep']==lid]
    if shoot_dir == 'left_to_right':
        min_or_max = lambda x : min(x)
    elif shoot_dir == 'right_to_left':
        min_or_max = lambda x : max(x)
    else:
        raise ValuError()
    
    sx = min_or_max(line['sx'])
    sy = line[line['sx'] == sx]['sy'].unique()
    if len(sy) > 1:
        raise ValueError('It is supposed to have only one value, rows being different receivers')
    else:
        sy = sy[0]
    
    return sx, sy
def arrow_for_shot_line(lid, md, dx):
    shoot_dir = find_shoot_dir(lid, md)
    sx, sy = find_marginal_shot(lid, md, shoot_dir)
    if shoot_dir == 'left_to_right':
        arrow = [sx-dx, sy, dx, 0]
    else:
        arrow = [sx+dx, sy, dx, 0]
    return arrow
def plot_out_data(p, it, sid, lid, freq, interleave,\
     phase=True, save=False, overwrite=1, **kwargs):
    dc = p.o.dumpcomp
    # print(dc)
    p.i.rnf.read_blocks(new_block_activator='freq') # NEEDED FOR THIS RUNFILE
    f = dc.it[it][sid]
    #     print('High-cut freq: %s Hz' % freq)
    f.split(overwrite=0)
    
    if interleave:
        plt.figure(figsize=(15,5))
        plt.suptitle('Interleaved syn and obs data (10xsyn, 10xobs, 10xsyn,...)')
        
        args = [f.obs.lid[lid]]
        kws = dict(norm='max', overwrite=overwrite)
        def set_ticks(c):
            chunk_size = 10
            ti = np.arange(len(f.syn.lid[lid]))[::chunk_size] - .5        
            plt.gca().set_xticks(ti)
            plt.gca().grid(axis='x', c=c, linestyle='-.', linewidth=2)            
        
        plt.subplot(121)
        f.syn.lid[lid].compare(*args, **kws, noextent=1)
        set_ticks('k')


        plt.subplot(122)
        f.syn.lid[lid].compare(*args, **kws, spect='ampl', dt=p.dt, cmap='hot', center_cmap=0)    
        plt.ylim(10,0)
        set_ticks('Grey')    
        if save:
            plt.savefig('ileave_p%s_it%s_sid%s_lid%s_freq%s.png' % (p.name, it, sid, lid, freq))
    if phase:        
        md = p.i.obs.read_header(overwrite=0)
        dx = kwargs.get('arrow_dx', 10*p.dx)
        if lid is None:
            f.plot_phase(freq, overwrite=overwrite, **kwargs)
        else:
            arrow_width = kwargs.get('arrow_width', p.dx)
            arrow_color = kwargs.get('arrow_color', 'k')
            arrow = arrow_for_shot_line(lid, md, dx=dx)
            f.plot_phase(freq, overwrite=overwrite, arrow=arrow,\
                arrow_width=arrow_width, arrow_color=arrow_color, **kwargs)

# VARIA
def add_spherical_anomaly(m, center=(100,235,50), radius=20, anom=-0.4, smooth=True):
    """
    center : tuple of 3 floats
      In nodes.
    radius : float
      In nodes.
    anom : float (+/-)
      Percent of model to add.
    """
    from fullwavepy.seismic.models import sphere
    from fullwavepy.ndat.arrays import Arr3d
    from scipy.ndimage import gaussian_filter as gauss_filt
    
    shape = m.read().shape
    sph = Arr3d(sphere(shape, center=center, radius=radius))
    if smooth:
        # c.shellen empirically:
        A = 100
        sigma = 10
        ##
        sphf = gauss_filt(A*sph, sigma) 
    else:
      sphf = sph
    sphf = sphf / np.max(sphf) * anom 
    sphf = Arr3d(sphf)
    
    x, y, z = center
    sphf.plot(x=x, y=y, z=z)
    
    tvp = m.array + sphf * m.array
    m.prep(Arr3d(tvp, extent=m.extent))
    return m
def tintegrate_wvlt(rsg, plot=False):
    """
    """
    from scipy.integrate import cumtrapz 
    from fullwavepy.ndat.arrays import Arr1d
    
    rs = Arr1d(rsg.read()[0,0])
    rs /= max(rs)
    RS = cumtrapz(rs)
    # APPARENTLY cumtrapz chops one time-sample that we will re-add manually
    # print(RS.shape)
    RS.resize(rs.shape, refcheck=False) # refcheck=True => can't do it
    # (help(RS.resize) missing entries are filled with zeros
    # DIFFERENT from np.resize !!) 
    RS /= max(RS)
    rsg.prep(wavelet=RS)
    if plot:
        plt.plot(rs)
        plt.plot(RS)
    #         qc_integration(rs, RS, rsg.read().dx[-1])
    return rsg
def correct_depth_WRONG(p, ioapi='sgy'): # use SP's addtodepth instead
    from fwipy.shell.specific.linux import run_linux
    assert ioapi == 'sgy'
    for suffix in ['-Sources', '-Receivers']:
        ext = '.geo'
        fname = p.i.path + p.name + suffix + ext 
        interfix = '_not_corrected'
        fname_nc = strip(fname) + interfix + ext
        command = 'cp ' + fname + ' ' + fname_nc
        run_linux(command)
        #print(fname, fname_nc)

        c = read_txt(fname_nc)
        header = c[0]
        header_str = ''
        for word in header:
            header_str += word + ' ' # FIXME: BETTER FORMATTING
        data = c[1: ]
        f = open(fname, 'w')
        f.write(header_str + '\n')

        for line in data:
            source_id = line[0]
            x = line[1]
            y = line[2]
            z = line[3]    
            if suffix == '-Sources':
    #             z = str(p.zsea * p.dx - float(z))
                z = str(p.zsea * p.dx + float(z))
            elif suffix == '-Receivers':
                z = str(p.zsea * p.dx + float(z))
            else:
                raise ValueError('suffix: ' + str(suffix))

            f.write(source_id + ' ' + x  + ' ' + y  + ' ' + z  + '\n')
        f.close()
    p.i.s.cat()
