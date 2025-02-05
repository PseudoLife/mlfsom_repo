from math import *
import os, time, math, tempfile
from os.path import join
import subprocess
import multiprocessing as mp
from glob import glob
import pandas as pd

mlfsom_path='~/Desktop/MLFSOM'
tmp_path = join(tempfile.gettempdir(), os.environ.get('USER'))
#tmp_path = '/tmp/peter/'  # tmp path of the operating sys


def GetExperimentList(N_grid, start_mos, k_mos, k_cell, k_bfactor, frames, \
    osc=0.01, exposure=0.5, xtal_size=200, beam_fwhm_x=100, beam_fwhm_y=100):
    """
    Used by the fn: SpacialDependentCrystal
    Returns a queue of params of exps, weights on the diffraction pattern,
    and frame number, for mlfsom to calculate the spacial and dose dependent diffraction
    patterns for a non-homogenous crystal model given by:
    
    N_grid by N_grid (even number) grid of cells that are centered on the beam
    start_mos: starting value of mosaicity
    k_mos: increase in mosaicity given by mos_0 + k_mos*dose
    k_cell: fractional increase in unit cell dim given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames
    osc: angle in degrees xtl rotates for each angle slice (0.01 is good for stills)
    exposure: exposure time per frame in sec (0.5 sec is mlfsom default)
    xtal_size: full crystal size in um
    beam_fwhm_x, beam_fwhm_y: in um  

    Note: dose = frame-no*decay-factor where decay-factor is exp(-(x**2+y**2)/(2*sigma**2))
    Returns tuple of list of params for exps to be calc by mlfsom:
    (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    and list of weights to each set of exp on each frame, frame_weights:
    (frame,weight,exp_ID)
    """
    # create list of unique cells and their weights in the crystal
    cell_list=[]    # form is (x,y,weight)
    if beam_fwhm_x == beam_fwhm_y:
        # less cells are needed due to symmetry
        for i in range(0,N_grid/2):
            for j in range(0,i+1):
                if i == j:
                    weight = 4
                else:
                    weight = 8
                cell_list.append((i+.5,j+.5,weight))
    else:
        # more cells are needed due to less symmetry
        for i in range(0,N_grid/2):
            for j in range(0,N_grid/2):
                weight = 4
                cell_list.append((i+.5,j+.5,weight))
    
    sub_xtal_size = xtal_size/float(N_grid)  # sub-crystal size - input for mlfsom
    sub_beam_size = 2.0 * sub_xtal_size  # beam size - input for mlfsom, must be greater than sub_xtal_size
    beam_sigma_x = beam_fwhm_x/2.355
    beam_sigma_y = beam_fwhm_y/2.355

    # calculate dose in each cell for each frame and params for each cell
    experiment_list = []
    # tuple (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    frame_weights = []
    # tuple (frame,weight,experiment_ID)
    for frame in range(0,frames):
        for cell in cell_list:
            ID = len(experiment_list)
            frame_weights.append((frame,cell[2],ID))  
            x_coord = cell[0]*sub_xtal_size
            y_coord = cell[1]*sub_xtal_size
            # decay factor is fractional decay in the int. of the incident beam relative to the max. int
            decay_factor = math.exp(-(x_coord/float(beam_sigma_x))**2/2.0 -(y_coord/float(beam_sigma_y))**2/2.0)
            dose = frame*decay_factor
            sub_exposure = exposure*decay_factor 
            experiment_list.append(\
                (ID, round(start_mos+k_mos*dose,3), round(k_bfactor*dose,2),\
                 round(k_cell*dose,4), osc, sub_exposure, sub_xtal_size, sub_beam_size))
    
    # experiment_list: (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    # frame_wights: (frame,weight,experiment_ID)
    return (frame_weights,experiment_list)


def GetHomogenousExperimentList(start_mos, k_mos, k_cell, k_bfactor, frames, osc=0.01, exposure=0.5):
    """
    Used by fn: HomogenousCrystal
    Returns a queue of params of exps for mlfsom to calc dose dependent diff.
    for a homogenously illuminated crystal model given by:
    
    start_mos: starting value of mos
    k_mos: increase in mosaicity given by mos_0 + k_mos*dose
    k_cell: fractional increase in unit cell dim. given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames
    osc: angle in degrees xtl rotates for each angle slice (0.01 is good for stills)
    exposure: exposure time per frame in sec

    Note: the dose in each frame is given by frame*1
    Returns list of params for exps to be calc by mlfsom:
    (ID,mos,bfactor_inc,cell_inc,osc,exposure,xtal_size,beam_size)
    """
    experiment_list=[]  # tuple(ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    xtal_size = 77.8  # some exp. constants for homogenous beam case
    beam_size = 100
    for frame in range(0,frames):
        dose = frame
        ID = len(experiment_list)
        experiment_list.append(\
            (ID, round(start_mos+k_mos*dose,3), round(k_bfactor*dose,2),\
             round(k_cell*dose,4), osc, exposure, xtal_size, beam_size))
    # experiment_list: (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    return experiment_list


def RunExperiments(prefix,experiment_list,resolution=2.5,solvent_B=35,threads=4):
    """
    Used by fn: SpacialDependentCrystal, HomogenousCrystal, RunFromExpQueue
    Runs each exp. and generates the diff pattern defined in experiment_list
    Saves imgs to mlfsom_path and necessary temp files to mlfsom_path/[prefix]
    experiment_list: tuple (ID,mos,bfactor_inc,cell_inc,osc)
    resolution: mlfsom default 2.5 A. initial resolution of the diffraction
    solvent_B: mlfsom default 35. does this change during the run at all??? needs more inspection
    threads: number of images that are run at once
    Images are saved as {prefix}###_001.img, where ### is the ID
    e.g. prefix = mseq and ID = 1, img => mseq1_001.img
    """
    time_initial = time.time()
    prev_cwd = os.getcwd()
    os.chdir(os.path.expanduser(mlfsom_path))
    os.system('cp 1H87.pdb temp.pdb')
    os.system('rm ' + join(tmp_path,'*')) # clear previous temp files
    
    # set up threading and run exps
    p = mp.Pool(threads)
    experiments_and_prefix = list(zip(experiment_list,[prefix]*len(experiment_list)))
    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('mkdir --parents ' + output_folder)  # create sub-fol to save output files
    for i in range(0,len(experiments_and_prefix),threads):
        # generate input files
        for exp in experiments_and_prefix[i:i+threads]:
            ID = exp[0][0]
            bfactor_inc = exp[0][2]
            cell_inc = exp[0][3]
            print("generating inputs for id="+str(ID))
            subprocess.call(['./change_pdb_param.com', 'temp.pdb', 'input'+str(ID)+'.pdb', \
                'add_cell='+str(cell_inc), 'add_bfactor='+str(bfactor_inc)])           
            subprocess.call(['./ano_sfall.com', 'energy=12660', \
                'input'+str(ID)+'.pdb', 'resolution='+str(resolution), 'solvent_B='+str(solvent_B)])
            subprocess.call(['mv', 'ideal_ano.mtz', 'input'+str(ID)+'.mtz'])
        p.map(RunExperiment,experiments_and_prefix[i:i+threads])
        
        # clear temp files and move results
        print('clearing temp files and moving results')
        os.system('mv ' + join(tmp_path,'mlfsom*.XYI ') + output_folder)    
        os.system('mv '+ join(tmp_path,'mlfsom*predin.txt ') + output_folder)        
        os.system('rm '+ join(mlfsom_path,'fit2d_*'))
        os.system('rm ' + join(tmp_path,'*'))

        # Function below might not be doing its job!!!
        RenameTempFiles(experiments_and_prefix[i:i+threads], output_folder)
    time_final = time.time()
    print 'TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)
    os.chdir(prev_cwd)


def RunExperiment(experiment_and_prefix):
    """
    Used by fn: RunExperiments
    Runs mlfsom with params defined by exp and outprefix
    experiment_and_prefix: two tuple
    first is exp params => (ID,mos,bfactor_inc,cell_inc,osc)
    second is prefix
    Images are saved as {prefix}###_001.img, where ### is the ID
    e.g. prefix = mseq and ID = 1, img => mseq1_001.img
    """
    # (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    experiment = experiment_and_prefix[0]
    prefix = experiment_and_prefix[1]
    ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size = experiment
    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('mkdir --parents ' + output_folder)  # create sub-fol to save output files
    # IS THE PRINT COMMAND BELOW REALLY NECESSARY???
    print ( ' '.join(['------- RUNNING ID='+str(ID), 'mos='+str(mos), 'bfactor_inc='+str(bfactor_inc),\
        'cell_inc='+str(cell_inc), 'osc='+str(osc), 'sub_exposure='+str(sub_exposure), \
        'sub_xtal_size='+str(sub_xtal_size), 'sub_beam_size='+str(sub_beam_size)]) )
    subprocess.call([\
        './mlfsom.com', join(output_folder,prefix+'_'+str(ID)+'_001.img'), 'input'+str(ID)+'.mtz',\
        'frames=1', 'id='+str(ID), 'mosaic='+str(mos), 'osc='+str(osc), 'exposure='+str(sub_exposure), \
        'xtal_size='+str(sub_xtal_size), 'beam_size='+str(sub_beam_size)])
    print('---- ID: '+str(ID)+' DONE ----')


def SpacialDependentCrystal(prefix,N_grid,start_mos,k_mos,k_cell,k_bfactor,\
    frames,resolution=2.5,solvent_B=35,osc=.01,exposure=0.5, \
    xtal_size=200, beam_fwhm_x=100, beam_fwhm_y=100,threads=4):
    """
    Generates files defining the params of exps, weights on the diff. pattern, and frame number
    for mlfsom to calc spacial and dose dep. diff. patterns for a non-homogenous xtl
    e.g. SpacialDependentCrystal(5,0.05,0.07,0.0153,2.228,Gaussian(10/14),14,1,"sp_seq1_")
    
    prefix: name of the exp
    N_grid by N_grid (even number) grid of cells that are centered on the beam
    start_mos: starting value of mos
    k_mos: increase in mos given by mos_0 + k_mos*dose
    k_cell: fractional increase in unit cell dim. given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames    
    resolution: res in A
    solvent_B: solvent B factor - how does this affect the results - needs more inspection 
    osc: angle in deg. xtl rotates for each angle slice (0.01 is good for stills)
    exposure: exposure time per frame in sec 
    xtal_size: total size of the crystal in um
    beam_fwhm_x, beam_fwhm_y: in um
    threads: number of images that are run at once
    
    Files generated:
    exp-{prefix}-desc.txt - lists the params of the exp
    exp-{prefix}-queue.txt - list of exps to be run, intended to be used in RunFromExpQueue
    exp-{prefix}-list.txt - back up list of exps
    Returns tuple of list of params for exps to be calc. by mlfsom given:
    (ID,mos,bfactor_inc,cell_inc,frames,osc)
    """
    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('mkdir --parents ' + output_folder)  # create sub-fol to save output files

    frame_weights,experiment_list = GetExperimentList(N_grid, start_mos, k_mos, k_cell, \
        k_bfactor, frames, osc, exposure, xtal_size, beam_fwhm_x, beam_fwhm_y)
    WriteDescription(prefix, N_grid, start_mos, k_mos, k_cell, k_bfactor, frames, resolution, \
    solvent_B, osc, exposure, xtal_size, beam_fwhm_x, beam_fwhm_y, threads, frame_weights)
    WriteExpQueueAndList(prefix,experiment_list)
    RunExperiments(prefix,experiment_list,resolution,solvent_B,threads)
    # Move files
    os.system('mv ' + join(mlfsom_path,'input*.pdb ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'input*.mtz ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'exp-'+prefix+'*.txt ') + output_folder)
    return experiment_list


def HomogenousCrystal(prefix,start_mos,k_mos,k_cell,k_bfactor,frames,\
    resolution=2.5,solvent_B=35,osc=.01,exposure=0.5,threads=4):
    """
    Generates files defining the params of exps for mlfsom to calc
    dose dep. diff. patterns for a homogenous xtl given by:
    prefix: name of the experiment
    start_mos: starting value of moasicty
    k_mos: increase in mosaicity given by mos_0 + k_mos*dosef
    k_cells: fractional increase in unit cell dim. given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames
    resolution: res in A
    solvent_B: solvent B factor - how does this affect the results - needs more inspection 
    osc: angle in deg. xtl rotates for each angle slice (0.01 is good for stills)
    exposure: exposure time per frame in sec
    threads: number of images that are run at once
    
    Returns list of params for exps to be calculated by mlffsom:
    (ID,mos,bfactor_inc,cell_inc,osc)
    Files generated:
    exp-{prefix}-desc.txt - lists the params of the exp
    exp-{prefix}-queue.txt - list of exps to be run, intended to be used in RunFromExpQueue
    exp-{prefix}-list.txt - back up list of exps
    Returns tuple of list of params for exps to be calc by mlfsom:
    (ID,mos,bfactor_inc,cell_inc,frames,osc)
    """
    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('mkdir --parents ' + output_folder)  # create sub-fol to save output files

    # some exp. constants for homogenous beam case
    N_grid = None
    xtal_size = 77.8
    beam_fwhm_x = 100
    beam_fwhm_y = 100
    frame_weights = []
    experiment_list = GetHomogenousExperimentList(\
        start_mos, k_mos, k_cell, k_bfactor, frames, osc, exposure)
    WriteDescription(prefix, N_grid, start_mos, k_mos, k_cell, k_bfactor, frames, resolution, \
    solvent_B, osc, exposure, xtal_size, beam_fwhm_x, beam_fwhm_y, threads, frame_weights)
    WriteExpQueueAndList(prefix,experiment_list)
    RunExperiments(prefix,experiment_list,resolution,solvent_B,threads)
    # Move files
    os.system('mv ' + join(mlfsom_path,'input*.pdb ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'input*.mtz ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'exp-'+prefix+'*.txt ') + output_folder)
    return experiment_list


def WriteExpQueueAndList(prefix,experiment_list):
    """
    Used by fn: HomogenousCrystal, SpacialDependentCrystal
    Writes list of exps to exp-{prefix}-queue.txt and exp-{prefix}-list.txt
    experiment_list: list of tuples defing exp params:
    (ID,mos,bfactor_inc,cell_inc,osc,sub_exposure,sub_xtal_size,sub_beam_size)
    prefix: name of exp
    """
    f=open("exp-"+prefix+"-queue.txt","w")
    flist=open("exp-"+prefix+"-list.txt","w")
    for exp in experiment_list:
        f.write(str(exp)+"\n")
        flist.write(str(exp)+"\n")
    f.close()
    flist.close()


def WriteDescription(\
    prefix, N_grid, start_mos, k_mos, k_cell, k_bfactor, frames, resolution, \
    solvent_B, osc, exposure, xtal_size, beam_fwhm_x, beam_fwhm_y, threads, frame_weights):
    """
    Used by fn: HomogenousCrystal, SpacialDependentCrystal
    Writes exp. params to a file e.g. exp-prefix-desc.txt
    """
    f=open("exp-"+prefix+"-desc.txt","w")
    f.write("prefix: "+prefix+"\n")
    f.write("N_grid: "+str(N_grid)+"\n")
    f.write("frames: "+str(frames)+"\n")
    f.write("start_mos: "+str(start_mos)+"\n")
    f.write("k_mos: "+str(k_mos)+"\n")
    f.write("k_cell: "+str(k_cell)+"\n")
    f.write("k_bfactor: "+str(k_bfactor)+"\n")
    f.write("resolution: "+str(resolution)+"\n")
    f.write("solvent_B: "+str(solvent_B)+"\n")
    f.write("osc: "+str(osc)+"\n")
    f.write("exposure: "+str(exposure)+"\n")
    f.write("xtal_size: "+str(xtal_size)+"\n")
    f.write("beam_fwhm_x: "+str(beam_fwhm_x)+"\n")
    f.write("beam_fwhm_y: "+str(beam_fwhm_y)+"\n")
    f.write("threads: "+str(threads)+"\n")
    f.write("exp_list_weights:"+"\n")
    f.write("(frame,weight,experiment_ID)"+"\n")
    
    for i in frame_weights:
        f.write(str(i)+"\n")
    f.close()


# This does not work correctly
def RenameTempFiles(experiments_and_prefix,folder):
    """
    NEEDS TO BE MODIFIED COMPLETELY, IT DID NOT WORK ANYWAYS!!!
    Used by fn: RunExperiments
    Renames temp txt files (peak intensities) generated by mlfsom to
    more meaningful file names with exp. params 
    """
    files=glob(folder+"/mlfsom*predin.txt")
    for exp_and_pre in experiments_and_prefix:
        prefix=exp_and_pre[1]
        exp=exp_and_pre[0]
        ID=exp[0]
        mos=exp[1]
        for file in files:
            with open(file) as f:
                lines=f.readlines()
            for l in lines:
                if "ID" in l:
                    fID=l.split()[-1]
                if "MOSAIC" in l:
                    fmos=l.split()[-1]
            if fmos==mos:
                print(ID,mos,fmos,file)


def CleanMessAfterStuck(prefix):
    """
    Standalone fn
    When HomogenousCrystal or SpacialDependentCrystal fn gets stuck, use this fn to
    1) move temp *.XYI and *predin.txt files and input*.pdb and mtz files to appropriate fol
    2) remove unnecessary files e.g. fit2d_001.in etc. and other tmp files  
    """
    prev_cwd = os.getcwd()
    os.chdir(os.path.expanduser(mlfsom_path))
    output_folder = join(mlfsom_path,'data_'+prefix)

    os.system('mv ' + join(tmp_path,'mlfsom*.XYI ') + output_folder)
    os.system('mv '+ join(tmp_path,'mlfsom*predin.txt ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'input*.pdb ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'input*.mtz ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'exp-'+prefix+'*.txt ') + output_folder)     
    os.system('rm '+ join(mlfsom_path,'fit2d_*'))
    os.system('rm ' + join(tmp_path,'*'))
    os.chdir(prev_cwd)


def RunFromQueue(exp_desc_file,exp_queue_file,threads=None):
    """
    Standalone fn
    When certain frames are not generated after using HomogenousCrystal or SpacialDependentCrystal...
    use this fn to rerun those exp.
    
    Takes desc file and queue file listing exps
    Then runs the specified exps just like HomogenousCrystal or SpacialDependentCrystal
    queue file can be created by modifying original exp-XXX-queue.txt file
    exp_queue_file: lists frames to be rerun
    threads: number of frames run at once - uses #threads in the desc file if not specified
    """
    # Read desc file in a dataframe, get the necessary params 
    df_desc = pd.read_csv(exp_desc_file,sep=': ',header=None,index_col=0,engine='python') ####### might need to set sep=': '
    prefix = df_desc.loc['prefix',1]
    resolution = df_desc.loc['resolution',1]
    solvent_B = df_desc.loc['solvent_B',1]
    if threads == None:
        threads = int(df_desc.loc['threads',1])

    # Read queue file, create experiment_list
    with open(exp_queue_file) as myfile:
        lines = myfile.readlines()
    experiment_list = [eval(line.strip()) for line in lines]

    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('mkdir --parents ' + output_folder)  # create sub-fol to save output files
    RunExperiments(prefix,experiment_list,resolution,solvent_B,threads)
    # Move files
    os.system('mv ' + join(mlfsom_path,'input*.pdb ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'input*.mtz ') + output_folder)
    os.system('mv ' + join(mlfsom_path,'exp-'+prefix+'*.txt ') + output_folder)
