import numpy as np
import os, sys
import struct
import fns_plot_data as Fplt
import datetime as dtm


def fn_bfile_info(bfile):
    # routine to get all output fields from binary files:

    # get info like dimensions and variable names from .b file
    bid    = open(bfile, 'r')
    lines = bid.readlines()
    bid.close()

    int_list = ['nx', 'ny', 'Nrecs', 'Norder'] # these are integers not floats
    binfo     = {'recnos': {}}  # dictionary with info about fields in corresponding .a file

    do_vlist = 0
    count     = 0
    for lin in lines:
        ls = lin.split()
        if ls != []:
            # skip blank lines

            # if not at "Record number and name"
            if not(ls[0]=='Record' and ls[1]=='number'):
                val = ls[0]
                key = ls[1]
                if not do_vlist:
                    if key in int_list:
                        val = int(val)
                    elif ("T" in val) and ("Z" in val):
                        val = dtm.datetime.strptime(val, "%Y%m%dT%H%M%SZ")
                    else:
                        val = float(val)
                    binfo.update({key:val})
                else:
                    binfo['recnos'].update({key:count})
                    count += 1
            else:
                # have got down to list of variables
                do_vlist = 1

    if binfo['Nrecs']!=len(binfo['recnos']):
        raise ValueError('Inconsistent number of records in file: '+bfile)

    return binfo


def get_array(fid, recno, nx, ny, fmt_size=4, order='F'):
    # routine to get the array from the .a (binary) file
    # * fmt_size = size in bytes of each entry)
    #    > default = 4 (real*4/single precision)
    recs      = nx*ny
    rec_size = recs*fmt_size
    fid.seek(recno*rec_size, 0) # relative to start of file
    #
    if fmt_size==4:
        fmt_py = 'f' # python string for single
    else:
        fmt_py = 'd' # python string for double

    data = fid.read(rec_size) # no of bytes to read
    fld  = struct.unpack(recs*fmt_py, data)
    fld  = np.array(fld)
    fld  = fld.reshape((nx, ny), order=order)

    return fld


def check_names(vname, variables, stop=True):

    if vname in variables:
        return vname

    lists = []

    # ice conc alt names
    List = ['ficem', 'fice', 'ice_conc', 'icec', 'cice', 'area', \
                        'concentration', 'sea_ice_concentration']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    # ice thick alt names
    List = ['hicem', 'hice', 'ice_thick', 'icetk', 'iceh', \
                'sea_ice_thickness', 'thickness']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    # floe size alt names
    List = ['dfloe', 'dmax', 'Dfloe', 'Dmax']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    # wave stress: x component
    List = ['taux', 'tau_x', 'taux_waves']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    # wave stress: y component
    List = ['tauy', 'tau_y', 'tauy_waves']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    # swh
    List = ['Hs', 'hs', 'swh', 'significant_wave_height']
    if vname in List:
        for v in variables:
            if v in List:
                return v

    if stop:
        print('Failed to get '+vname)
        print('\nAvailable variables:')
        for v in variables:
            print(v)

        print('')
        raise ValueError(vname+' not in variable list')
    else:
        return ''


def fn_read_general_binary(afile, vlist=None):
    # routine to get output fields from binary files:

    # get dimensions and variable names from .b file
    bfile = afile[:-2]+'.b'
    binfo = fn_bfile_info(bfile)

    if vlist is None:
        # get all fields in binary files
        recnos = binfo['recnos']
    else:
        # check vbls are in binary files
        recnos = {}
        for vbl in vlist:
            vname = check_names(vbl, list(binfo['recnos']), stop=True)
            recnos.update({vbl:binfo['recnos'][vname]})

    nv = binfo['Nrecs']
    nx = binfo['nx']
    ny = binfo['ny']
    if binfo['Norder']==1:
        order = 'fortran'
    else:
        order = 'C'

    # can now read data from .a file
    sz         = os.path.getsize(afile)
    fmt_size = int(sz/float(nv*nx*ny)) # 4 for single,  8 for double
    aid        = open(afile, 'rb')

    out    = {}
    for vbl in recnos:
        out.update({vbl: get_array(aid, recnos[vbl], nx, ny, order=order, fmt_size=fmt_size)})

    aid.close()

    # outputs
    return out, binfo


class var_info:
    # simple object
    # scalars: obj.value, obj.unit
    # numpy arrays: obj.data, obj.unit, obj.length
    def __init__(self, value, unit=None):

        # set units
        if unit is not '':
            self.unit = unit
        else:
            self.unit = None

        if hasattr(value, 'ndim'):
            # an array
            self.data    = value
            self.length = len(value)
            self.min     = np.min(value)
            self.max     = np.max(value)
            self.first  = value[0]
            self.last    = value[-1]
        else:
            # a number
            self.value = value

        return


def fn_check_grid(outdir):
    # routine to get grid and other parameters
    # from binary files

    afile = outdir+'/wim_grid.a'
    if not os.path.exists(afile):
        return dict()

    grid_prams, info = fn_read_general_binary(afile)

    # extra info
    nx, ny = grid_prams['X'].shape
    grid_prams.update({'nx':nx})
    grid_prams.update({'ny':ny})
    #
    dx = np.mean(grid_prams['scvx'])
    dy = np.mean(grid_prams['scuy'])
    grid_prams.update({'dx':dx})
    grid_prams.update({'dy':dy})

    return grid_prams

def fn_check_init(outdir):
    # routine to get initial fields from binary files:
    lst    = os.listdir(outdir)
    afile = None
    for f in lst:
        if 'wim_init' in f and '.a' in f:
            afile = outdir+'/'+f
            break

    ## ice fields
    ice_fields = fn_read_general_binary(afile, vlist=['icec', 'iceh', 'dfloe'])[0]

    # ice mask
    ice_fields.update({'ICE_MASK':0*ice_fields['icec']})
    ice_fields['ICE_MASK'][ice_fields['icec']>0.05] = 1.0


    ## wave fields
    wave_fields = fn_read_general_binary(afile, vlist=['Hs', 'Tp', 'mwd'])[0]

    # wave mask
    wave_fields.update({'WAVE_MASK':0*wave_fields['Hs']})
    wave_fields['WAVE_MASK'][wave_fields['Hs']>0.0] = 1.0


    return ice_fields, wave_fields


def fn_check_out_bin(outdir, **kwargs):
    # routine to get output fields from binary files:
    lst    = os.listdir(outdir)
    afile = None
    for f in lst:
        if 'wim_out' in f and '.a' in f:
            afile = outdir+'/'+f
            break

    return fn_read_general_binary(afile, vlist=['dfloe', 'taux', 'tauy', 'Hs', 'Tp', 'mwd'])[0]

def fn_check_out_arr(out_arrays):
    # routine to convert out_arrays to dictionary
    out_fields = {}

    keys = ['dfloe', 'taux', 'tauy', 'Hs', 'Tp', 'mwd']
    for n, key in enumerate(keys):
        out_fields.update({key:out_arrays[:, :, n]})

    # outputs
    return out_fields


def fn_check_prog(outdir, cts):
    # routine to get progress fields from binary files:
    # cts is a string eg '010' or '0010' corresponding to the time step
    if type(cts)==type(0):
        # convert from int to str of correct length
        fils  = os.listdir(outdir+'/binaries/prog/')
        # cts0  = fils[0].strip('wim_prog')[:-2]
        n = 0
        while '.swp' in fils[n] or '.DS_Store' in fils[n]:
             n = n+1
        cts0 = fils[n].strip('wim_prog')[:-2]
        fmt  = '%'+str(len(cts0))+'.'+str(len(cts0))+'d'
        cts  = fmt %(cts)
        print(cts0, cts)

    afile = outdir+'/binaries/prog/wim_prog'+cts+'.a'
    return fn_read_general_binary(afile, vlist=['dfloe', 'taux', 'tauy', 'Hs', 'Tp', 'mwd'])[0]

class file_list:

    def __init__(self, directory, pattern,
            grid_prams = None, grid_dir = None):
        self.dirname   = directory
        self.pattern   = pattern
        self.times     = []
        self.files     = []
        self.Nfiles    = 0
        self.variables = []

        # check if directory exists before continuing
        if not os.path.exists(self.dirname):
            return

        self.grid_prams = grid_prams
        if grid_prams is None:
            if grid_dir is not None:
                gp = fn_check_grid(grid_dir)
                if len(gp)>0:
                    self.grid_prams = gp
            if self.grid_prams is None:
                gp = fn_check_grid(self.dirname)
                if len(gp)>0:
                    self.grid_prams = gp

        # if directory exists,  sort the files
        all_files = os.listdir(self.dirname)

        # find the .a files
        alist    = []
        steplist = []
        for pf in all_files:
            if (pattern in pf) and ('.a'==pf[-2:]):
                alist.append(pf)
                stepno = pf.strip(pattern).strip('.a') # step no eg 001 or date yyyymmddThhmmssZ
                if ('T' in stepno) and ('Z' in stepno):
                    # date
                    steplist.append(
                            dtm.datetime.strptime(stepno, '%Y%m%dT%H%M%SZ'))
                else:
                    # integer step no
                    steplist.append(int(stepno))

        # sort according to steplist:
        slist = sorted([(e, i) for i, e in enumerate(steplist)])
        self.times  = [e for e, i in slist]
        self.files  = [alist[i] for e, i in slist]
        self.Nfiles = len(alist)

        if self.Nfiles>0:
            bfile = self.files[0].replace('.a', '.b')
            info  = fn_bfile_info(self.dirname+'/'+bfile)
            self.variables = list(info['recnos'])

        return

    def plot_steps(self, grid_prams, figdir3, zlims=None, vlist=None, **kwargs):
        pdir = self.dirname

        if not os.path.exists(figdir3):
            os.mkdir(figdir3)

        # check variables to stop crash
        if vlist is not None:
            Vlist = []
            for vbl in vlist:
                if vbl in self.variables:
                    Vlist.append(vbl)
            if len(Vlist)==0:
                print("Variables in 'vlist' argument not present in files")
                print("Variables that are available:")
                for vbl in self.variables:
                    print("    "+vbl)
                return
        else:
            Vlist = None


        # determine the plotting limits
        if zlims is None:
            # set colorbar axes automatically
            zdef  = [1.e30, -1.e30]
            zlims = {'icec' :1*zdef, 
                         'iceh' :1*zdef, 
                         'Dmax' :1*zdef, 
                         'tau_x':1*zdef, 
                         'tau_y':1*zdef, 
                         'Hs'    :1*zdef, 
                         'Tp'    :1*zdef, 
                         'mwd'  :1*zdef}

            # determine the plotting limits
            # by checking the arrays
            alist = self.files
            tlist = []
            for pf in alist:
                afile = pdir+'/'+pf
                #
                fields = fn_read_general_binary(afile, vlist=Vlist, **kwargs)[0]
                for key in fields:
                    key2 = check_names(key, list(zlims), stop=False)
                    if key2!="":
                        zmin = fields[key].min()
                        zmax = fields[key].max()
                        if zmin<zlims[key2][0]:
                            zlims[key2][0] = zmin
                        if zmax>zlims[key2][1]:
                            zlims[key2][1] = zmax
                #
                tlist.append('_'+pf[4:-2])

        for vbl in fields:
            key2 = check_names(key, list(zlims), stop=False)
            if key2!="":
                print("range in "+vbl, zlims[key2][0], zlims[key2][1])


        # do the plots
        for vbl in fields:
            print('\n')
            print('Plotting results for '+vbl)
            print('\n')

            figdir3B = figdir3+'/'+vbl
            if not os.path.exists(figdir3B):
                os.mkdir(figdir3B)

            for i, pf in enumerate(alist):
                print(self.times[i])
                afile  = pdir+'/'+pf
                fields = fn_read_general_binary(afile)[0]
                Fplt.fn_plot_gen(grid_prams, fields, figdir=figdir3B, \
                        zlims_in=zlims, text=tlist[i], vlist=[vbl])
        return


    def plot_step(self, grid_prams, time_index=0, vlist=None,
            figdir=None, sep_fig_dirs=True, **kwargs):

        pf    = self.files[time_index]
        afile = self.dirname+'/'+pf
        tstr  = pf[4:-2]

        DO_SAVE = (figdir is not None)
        if DO_SAVE:
            if not os.path.exists(figdir):
                os.mkdir(figdir)

        # do the plots
        fields = fn_read_general_binary(afile, vlist=vlist)[0]

        if not sep_fig_dirs:
            Fplt.fn_plot_gen(grid_prams, fields, figdir=figdir,
                 text=tstr, **kwargs)
        else:
            for vbl in fields:
                subdir = figdir
                if DO_SAVE:
                    subdir = figdir+'/'+vbl
                    if not os.path.exists(subdir):
                        os.mkdir(subdir)
                else:
                    subdir = None
                Fplt.fn_plot_gen(grid_prams, {vbl: fields[vbl]},
                        figdir=subdir, text=tstr, **kwargs)
        return

class wim_results:

    def __init__(self, outdir='.', grid_dir=None):

        self.rootdir = outdir
        self.bindir  = outdir+'/binaries'
        if not os.path.exists(self.bindir):
            raise ValueError(self.bindir+ ' does not exist')


        # check for grid files
        self.grid_dir = None
        if grid_dir is not None:
            lst = os.listdir(grid_dir)
            for f in lst:
                if 'wim_grid' in f and '.a' in f:
                    self.grid_dir = grid_dir
                    break
            
            if self.grid_dir is None:
                print('wim_grid.[a, b] not in '+grid_dir)
                raise ValueError('input "grid_dir" does not contain grid files')
            else:
                self.grid_dir = grid_dir

        else:
            dirs = [self.bindir, outdir+'/../grid']
            i = 0

            while (self.grid_dir is None) and (i<len(dirs)):
                dir_i = dirs[i]
                lst   = os.listdir(dir_i)
                for f in lst:
                    if 'wim_grid' in f and '.a' in f:
                        self.grid_dir = dir_i
                        break
                i += 1

            if self.grid_dir is None:
                raise ValueError('wim_grid.[a, b] not found')

        # check for initial conditions
        binlist = os.listdir(self.bindir)
        self.init_list = file_list(self.bindir, 'wim_init')

        if self.init_list.Nfiles==0:
            # try again in binaries/init
            # print("try again in binaries/init")
            self.init_list = file_list(self.bindir+'/init', 'wim_init')

        if self.init_list.Nfiles>0:
            self.start_time = self.init_list.times[0]


        # check for final conditions
        binlist         = os.listdir(self.bindir)
        self.out_list = file_list(self.bindir, 'wim_out')

        if self.out_list.Nfiles==0:
            # try again in binaries/out
            # print("try again in binaries/final")
            self.out_list = file_list(self.bindir+'/final', 'wim_out')

        if self.out_list.Nfiles>0:
            self.finish_time = self.out_list.times[-1]


        # Check for progress files
        self.prog_dir  = self.bindir+'/prog'
        self.prog_list = file_list(self.prog_dir, 'wim_prog')


        # Check for incident wave files
        self.inc_dir  = self.bindir+'/incwaves'
        self.inc_list = file_list(self.inc_dir, 'wim_inc')


        # final checks
        if (self.init_list.Nfiles==0) and\
            (self.out_list.Nfiles==0)  and\
            (self.prog_list.Nfiles==0):
            raise ValueError('No results in ' + outdir)

        # where to save figures
        self.figdir = outdir+'/figs'

        return


    def get_grid(self):
        return fn_check_grid(self.grid_dir)


    def get_fields(self, field_type, time_index=0, **kwargs):

        if field_type=="initial":
            file_list = self.init_list
            short     = "init"
        elif field_type=="progress":
            file_list = self.prog_list
            short     = "prog"
        elif field_type=="final":
            file_list = self.out_list
            short     = "out"
        elif field_type=="inc":
            file_list = self.inc_list
            short     = "inc"
        else:
            raise ValueError("unknown field_type: "+field_type)

        if file_list.Nfiles==0:
            Start = field_type[0].upper + field_type[1:]
            raise ValueError(Start+' files not outputted: wim_'+short+'*.[a, b]')

        afile = file_list.dirname+'/'+file_list.files[time_index]
        return fn_read_general_binary(afile, **kwargs)


    def initial_fields(self, **kwargs):
        return self.get_fields("initial", **kwargs)


    def out_fields(self, **kwargs):
        return self.get_fields("final", **kwargs)


    def progfields(self, **kwargs):
        return self.get_fields("progress", **kwargs)


    def incfields(self, **kwargs):
        return self.get_fields("inc", **kwargs)


    def plot(self, time_index=None, show=False, field_type="initial", vlist=None, **kwargs):

        if field_type=="initial":
            file_list = self.init_list
            short     = "init"
            outdir    = "init"
        elif field_type=="progress":
            file_list = self.prog_list
            short     = "prog"
            outdir    = "prog"
        elif field_type=="final":
            file_list = self.out_list
            short     = "out"
            outdir    = "final"
        elif field_type=="inc":
            file_list = self.inc_list
            short     = "inc"
            outdir    = "incwaves"
        else:
            raise ValueError("unknown field_type: "+field_type)


        if file_list.Nfiles==0:
            print('No '+field_type+' files wim*'+short+'*.[a, b] in '+
                    file_list.dirname)
            print('Not plotting')
            return


        print('\nPLOTTING '+field_type.upper()+' FILES...\n')


        figdir3 = self.figdir
        if not os.path.exists(figdir3):
            os.mkdir(figdir3)
        figdir3 = self.figdir+'/'+outdir
        if not os.path.exists(figdir3):
            os.mkdir(figdir3)


        grid_prams = self.get_grid()
        if time_index is None:
            # plot all
            file_list.plot_steps(grid_prams, figdir3=figdir3, vlist=vlist, **kwargs)
            print('\nPlots in '+figdir3+'\n')
        else:
            if show:
                figdir3 = None
                # plot step & show
                # - otherwise plot step & save fig
            file_list.plot_step(grid_prams, time_index=time_index, 
                    figdir3=figdir3, vlist=vlist, **kwargs)

        return


    def plot_grid(self, show=False, **kwargs):

        print('\nPLOTTING LAND MASK...\n')

        if not show:
            hold   = False
            figdir = self.figdir
            if not os.path.exists(figdir):
                os.mkdir(figdir)
            figdir += '/grid'
            if not os.path.exists(figdir):
                os.mkdir(figdir)
        else:
            figdir = None
            hold   = True

        grid_prams  = self.get_grid()
        Fplt.fn_plot_gen(grid_prams, grid_prams, vlist=['LANDMASK'], \
                figdir=figdir, hold=hold)

        return


    def plot_initial(self, **kwargs):
        self.plot(field_type="initial", **kwargs)
        return


    def plot_prog(self, **kwargs):
        self.plot(field_type="progress", **kwargs)
        return


    def plot_final(self, **kwargs):
        self.plot(field_type="final", **kwargs)
        return


    def plot_incwaves(self, **kwargs):
        self.plot(field_type="inc", **kwargs)
        return
