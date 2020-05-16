import os, subprocess, time, math, tempfile, shutil
from os.path import join
from glob import glob
#import multiprocessing as mp
import pathos.multiprocessing as mp
import pandas as pd

mlfsom_path=os.path.expanduser('~/Desktop/MLFSOM')
tmp_path = join(tempfile.gettempdir(), os.environ.get('USER'))


def GetGaussianExperimentList(N_grid, start_mos, k_mos, k_cell, k_bfactor, frames, \
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


def GetHomogenousExperimentList(stills, start_mos, k_mos, k_cell, k_bfactor, frames, phi_start, osc, distance, flux):
    """
    Used by fn: HomogenousCrystal
    Returns a queue of params of exps for mlfsom to calc dose dependent diff.
    for a homogenously illuminated crystal model given by:
    
    stills: True for stills, False for continuous data collection e.g. full data set
    start_mos: starting value of mos
    k_mos: increase in mosaicity given by mos_0 + k_mos*dose
    k_cell: fractional increase in unit cell dim. given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames
    osc: angle in degrees xtl rotates for each angle slice (0.01 is good for stills)
    distance: in mm e.g. 250 is good for lysozyme
    flux: in ph/s

    Note: the dose in each frame is given by frame*1
    Returns list of exp. params to be calc by mlfsom:
    (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    """
    experiment_list=[]  # tuple(ID,mos,bfactor_inc,cell_inc,osc,sub_flux,sub_xtal_size,sub_beam_size)
    xtal_size = 77.8  # some exp. constants for homogenous beam case
    beam_size = 100
    for frame in range(0,frames):
        dose = frame
        ID = len(experiment_list) + 1    # ID of the first exp is set to 1
        if stills == True:
            phi = phi_start
        else:
            phi = phi_start + frame*osc
        experiment_list.append(\
            (ID, round(start_mos+k_mos*dose,5), round(k_bfactor*dose,3),\
             round(k_cell*dose,6), osc, round(phi,2), round(distance,4), "%.2e" %flux, xtal_size, beam_size))
    # experiment_list: (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    return experiment_list


def RunExperiments(pdb_file,prefix,experiment_list,resolution,solvent_B,threads):
    """
    Used by fn: SpacialDependentCrystal, HomogenousCrystal, RunFromExpQueue
    Runs each exp. and generates the diff pattern defined in experiment_list
    Saves imgs to mlfsom_path and necessary temp files to mlfsom_path/[prefix]

    pdb_file: input pdb file e.g. '1H87.pdb'
    prefix: name of the exp
    experiment_list: tuple (ID,mos,bfactor_inc,cell_inc,osc)
    resolution: mlfsom default 2.5 A. initial resolution of the diffraction
    solvent_B: mlfsom default 35. does this change during the run at all??? needs more inspection
    threads: number of images that are run at once
    Images are saved as {prefix}###_001.img, where ### is the ID
    e.g. prefix = mseq and ID = 1, img => mseq1_001.img
    """
    time_initial = time.time()
    prev_cwd = os.getcwd()
    os.system('rm -rf ' + join(tmp_path,'*')) # clear previous temp files
    
    # set up threading and run exps
    pool = mp.Pool(threads)
    
    experiments_and_prefix = list(zip(experiment_list,[prefix]*len(experiment_list)))
    output_folder = join(mlfsom_path,'data_'+prefix)
    os.system('rm -rf ' + output_folder)  # delete the same named folder if exists
    os.system('mkdir --parents ' + output_folder)  # create fol to save output files
    for i in range(0,len(experiments_and_prefix),threads):
        # generate input files
        for exp in experiments_and_prefix[i:i+threads]:
            ID = exp[0][0]

            output_subfolder = join(output_folder,'exp_'+str(ID))
            os.system('mkdir ' + output_subfolder)
            os.chdir(output_subfolder)

            subprocess.call(['cp', join(mlfsom_path,pdb_file), 'refined.pdb'])
            os.system( 'cp %s %s' %(join(mlfsom_path,'required/*'),output_subfolder) )

            bfactor_inc = exp[0][2]
            cell_inc = exp[0][3]
            print("\nGenerating inputs for ID="+str(ID))
            
            # Generate inputID.pdb for each exp using refined.pdb which is just the copy of original input pdb
            subprocess.call([ join(mlfsom_path,'change_pdb_param.com'), 'refined.pdb', 'input'+str(ID)+'.pdb', \
                'add_cell='+str(cell_inc), 'add_bfactor='+str(bfactor_inc) ])        
            
            # Generate inputID.mtz from inputID.pdb by adding solvent contribution  
            subprocess.call([ join(mlfsom_path,'ano_sfall.com'), 'input'+str(ID)+'.pdb', \
                'resolution='+str(resolution), 'solvent_B='+str(solvent_B) ])
            subprocess.call(['mv', 'ideal_ano.mtz', 'input'+str(ID)+'.mtz'])
    
    pool.map(RunExperiment,experiments_and_prefix)
    pool.close()
    pool.join()
    os.chdir(prev_cwd)
    print('Clearing temp files and moving results\n')

    tmp_sub_fols = [x for x in os.listdir(tmp_path) if x.startswith('exp_')]
    for exp_fol in tmp_sub_fols:
        os.system('mv ' + join(tmp_path,exp_fol,'mlfsom*.XYI ') + output_folder)
        os.system('mv '+ join(tmp_path,exp_fol,'mlfsom*predin.txt ') + output_folder)

    data_sub_fols = [x for x in os.listdir(output_folder) if x.startswith('exp_')]
    for exp_fol in data_sub_fols:
        os.system('mv ' + join(output_folder,exp_fol,'*.img ') + output_folder)
        os.system('mv ' + join(output_folder,exp_fol,'input*.pdb ') + output_folder)
        os.system('mv ' + join(output_folder,exp_fol,'input*.mtz ') + output_folder)
        os.system('rm -rf '+ join(output_folder,exp_fol))

    os.system('rm -rf ' + join(tmp_path,'*'))
    time_final = time.time()
    print 'TOTAL TIME ELAPSED: %i min' %int((time_final-time_initial)/60)


def RunExperiment(experiment_and_prefix):
    """
    Used by fn: RunExperiments
    Runs mlfsom with params defined by exp and outprefix
    experiment_and_prefix: two tuple
    first is exp params => (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    second is prefix
    Images are saved as {prefix}###_001.img, where ### is the ID
    """
    # (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    experiment = experiment_and_prefix[0]
    prefix = experiment_and_prefix[1]
    ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size = experiment
    
    output_folder = join(mlfsom_path,'data_'+prefix)
    output_subfolder = join(output_folder,'exp_'+str(ID))
    os.chdir(output_subfolder)

    # Create sub folders under tmp dir to place temp files of each exp in separate fol
    tmp_exp_path = join(tmp_path,'exp_'+str(ID))
    os.system('mkdir '+tmp_exp_path)
    os.environ["CCP4_SCR"] = tmp_exp_path
       
    print('------------------ RUNNING ID: '+str(ID)+' ------------------')
    subprocess.call([\
        './mlfsom.com', \
        prefix+'_'+str(ID)+'_001.img', \
        'input'+str(ID)+'.mtz',\
        'frames=1', \
        'id='+str(ID), \
        'mosaic='+str(mos), \
        'osc='+str(osc), \
        'phi='+str(phi), \
        'distance='+str(distance), \
        'flux='+str(sub_flux), \
        'xtal_size='+str(sub_xtal_size), \
        'beam_size='+str(sub_beam_size)])
    print('------------------ ID: '+str(ID)+' DONE ------------------')


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

    frame_weights,experiment_list = GetGaussianExperimentList(N_grid, start_mos, k_mos, k_cell, \
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


def HomogenousCrystal(pdb_file,prefix,stills,start_mos,k_mos,k_cell,k_bfactor,frames,\
    resolution=2.5,solvent_B=35,phi_start=342,osc=.01,distance=250,flux=8.4e10,threads=4):
    """
    Generates files defining the params of exps for mlfsom to calc
    dose dep. diff. patterns for a homogenous xtl given by:
    
    pdb_file: input pdb file e.g. 1H87.pdb
    prefix: name of the experiment
    stills: True for stills, False for continuous  data collection e.g. full data set
    start_mos: starting value of moasicty
    k_mos: increase in mosaicity given by mos_0 + k_mos*dosef
    k_cells: fractional increase in unit cell dim. given by cell_dimenion_0*(1+k_cell*dose)
    k_bfactor: increase in b-factors given by B_0 + k_bfactor*dose
    frames: number of frames
    resolution: res in A
    solvent_B: solvent B factor - how does this affect the results - needs more inspection 
    osc: angle in deg. xtl rotates for each angle slice (0.01 is good for stills)
    flux: in ph/s
    threads: number of images that are run at once
    
    Files generated:
    exp-{prefix}-desc.txt - lists the params of the exp
    exp-{prefix}-queue.txt - list of exps to be run, intended to be used in RunFromQueue
    exp-{prefix}-list.txt - back up list of exps
    Returns tuple of list of params for exps to be calc by mlfsom:
    (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    """
    # exp. constants for homogenous beam case
    N_grid = None
    xtal_size = 77.8
    beam_fwhm_x = 100
    beam_fwhm_y = 100
    frame_weights = []
    experiment_list = GetHomogenousExperimentList(\
        stills, start_mos, k_mos, k_cell, k_bfactor, frames, phi_start, osc, distance, flux)
    RunExperiments(pdb_file,prefix,experiment_list,resolution,solvent_B,threads)
    WriteDescription(pdb_file, prefix, stills, N_grid, start_mos, k_mos, \
        k_cell, k_bfactor, frames, resolution, solvent_B, phi_start, osc, distance, flux, \
        xtal_size, beam_fwhm_x, beam_fwhm_y, threads, frame_weights)
    WriteExpQueueAndList(prefix,experiment_list)
    return experiment_list


def WriteExpQueueAndList(prefix,experiment_list):
    """
    Used by fn: HomogenousCrystal, SpacialDependentCrystal
    Writes list of exps to exp-{prefix}-queue.txt and exp-{prefix}-list.txt
    experiment_list: list of tuples defing exp params:
    (ID,mos,bfactor_inc,cell_inc,osc,phi,distance,sub_flux,sub_xtal_size,sub_beam_size)
    prefix: name of exp
    """
    output_folder = join(mlfsom_path,'data_'+prefix)
    f_queue=open(join(output_folder,"exp-"+prefix+"-queue.txt"),"w")
    f_list=open(join(output_folder,"exp-"+prefix+"-list.txt"),"w")
    for exp in experiment_list:
        f_queue.write(str(exp)+"\n")
        f_list.write(str(exp)+"\n")
    f_queue.close()
    f_list.close()


def WriteDescription(\
    pdb_file, prefix, stills, N_grid, start_mos, k_mos, k_cell, k_bfactor, \
    frames, resolution, solvent_B, phi_start, osc, distance, flux, \
    xtal_size, beam_fwhm_x, beam_fwhm_y, threads, frame_weights):
    """
    Used by fn: HomogenousCrystal, SpacialDependentCrystal
    Writes exp. params to a file e.g. exp-prefix-desc.txt
    """
    output_folder = join(mlfsom_path,'data_'+prefix)
    f = open(join(output_folder,"exp-"+prefix+"-desc.txt"),"w")
    f.write("pdb_file: "+pdb_file+"\n")
    f.write("prefix: "+prefix+"\n")
    f.write("stills: "+str(stills)+"\n")
    f.write("N_grid: "+str(N_grid)+"\n")
    f.write("frames: "+str(frames)+"\n")
    f.write("start_mos: "+str(start_mos)+"\n")
    f.write("k_mos: "+str(k_mos)+"\n")
    f.write("k_cell: "+str(k_cell)+"\n")
    f.write("k_bfactor: "+str(k_bfactor)+"\n")
    f.write("resolution: "+str(resolution)+"\n")
    f.write("solvent_B: "+str(solvent_B)+"\n")
    f.write("phi_start: "+str(phi_start)+"\n")
    f.write("osc: "+str(osc)+"\n")
    f.write("distance: "+str(distance)+"\n")
    f.write("flux: "+"%.2e" %flux + "\n")
    f.write("xtal_size: "+str(xtal_size)+"\n")
    f.write("beam_fwhm_x: "+str(beam_fwhm_x)+"\n")
    f.write("beam_fwhm_y: "+str(beam_fwhm_y)+"\n")
    f.write("threads: "+str(threads)+"\n")
    f.write("exp_list_weights:"+"\n")
    f.write("(frame,weight,experiment_ID)"+"\n")
    
    for i in frame_weights:
        f.write(str(i)+"\n")
    f.close()
